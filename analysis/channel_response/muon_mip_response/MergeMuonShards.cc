// Merge disjoint channel workers; recompute medians from full-precision values.
#include "MuonMIPCommon.hh"
#include "MuonMIPPlotStyle.hh"
#include "../adc_energy/FastPng.hh"
#include <TFileMerger.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TParameter.h>
#include <iomanip>
#include <iostream>
#include <tuple>
using namespace mmr;
namespace {
struct Table {std::string header;std::vector<std::string> rows;std::map<std::string,size_t> columns;};
Table readTable(const fs::path& path){std::ifstream f(path);require(bool(f),"Missing worker table: "+path.string());Table t;require(bool(std::getline(f,t.header)),"Empty worker table: "+path.string());auto h=split(t.header);for(size_t i=0;i<h.size();++i)t.columns.emplace(h[i],i);for(std::string line;std::getline(f,line);){if(line.empty())continue;require(split(line).size()==h.size(),"Incomplete worker table row: "+path.string());t.rows.push_back(line);}return t;}
std::ofstream output(const fs::path& path){std::ofstream f(path);require(bool(f),"Cannot write "+path.string());f<<std::setprecision(12);return f;}
std::string value(const Table&t,const std::string&row,const std::string&column){auto fields=split(row);require(t.columns.count(column),"Missing table column "+column);return fields.at(t.columns.at(column));}
Table joinTables(const std::vector<fs::path>& parts,const fs::path& out,const std::string& name,bool optional=false){
 Table all;for(const auto&part:parts){auto path=part/name;if(!fs::exists(path)){require(optional,"Missing "+path.string());continue;}auto t=readTable(path);if(all.header.empty()){all.header=t.header;all.columns=t.columns;}require(all.header==t.header,"Worker table headers differ: "+name);all.rows.insert(all.rows.end(),t.rows.begin(),t.rows.end());}
 if(all.header.empty())return all;
 std::stable_sort(all.rows.begin(),all.rows.end(),[](const std::string&a,const std::string&b){return std::stoll(a.substr(0,a.find('\t')))<std::stoll(b.substr(0,b.find('\t')));});
 auto f=output(out/name);f<<all.header<<'\n';for(const auto&row:all.rows)f<<row<<'\n';return all;
}
double median(std::vector<double> a){if(a.empty())return NAN;std::sort(a.begin(),a.end());return a.size()%2?a[a.size()/2]:.5*(a[a.size()/2-1]+a[a.size()/2]);}
struct Record{int sipm=0,fallback=0,valid=0;double normalized=NAN;};
} // namespace
int main(int argc,char**argv)try{
 require(argc==4,"Usage: MergeMuonShards parts_directory config.tsv output_directory");Config config(argv[2]);fs::path partsBase=argv[1],out=argv[3];gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetCanvasPreferGL(false);
 std::vector<std::pair<int,fs::path>> indexed;for(const auto&entry:fs::directory_iterator(partsBase)){if(!entry.is_directory())continue;std::string name=entry.path().filename();if(name.rfind("worker_",0)==0)indexed.emplace_back(std::stoi(name.substr(7)),entry.path());}std::sort(indexed.begin(),indexed.end());require(!indexed.empty(),"No completed fit workers");std::vector<fs::path> parts;for(const auto&v:indexed)parts.push_back(v.second);
 require(!fs::exists(out/"spectra.root"),"Merged ROOT output already exists");
 std::vector<Table> statuses;for(const auto&part:parts){statuses.push_back(readTable(part/"analysis_status.tsv"));require(statuses.back().rows.size()==1,"Worker status is incomplete");}
 for(const auto&name:{"fits.tsv","comparison.tsv","efficiency.tsv","calibration_closure.tsv","plateau_stability.tsv"})joinTables(parts,out,name);
 joinTables(parts,out,"dependency_fits.tsv",!config.integer("dependency-fits",1));
 joinTables(parts,out,"performance.tsv");
 // Group medians must use the original doubles, never rounded public TSVs.
 std::map<std::tuple<int,int,int>,Record> records;
 double commonTarget=NAN;
 for(const auto&part:parts){auto file=openFile((part/"spectra.root").string());if(!std::isfinite(commonTarget)){auto*p=dynamic_cast<TParameter<double>*>(file->Get("geometry_common_target"));require(p,"Missing common matching target");commonTarget=p->GetVal();}
  auto*t=tree(*file,"WorkerSummary");int cid=0,sample=0,matched=0,sipm=0,fallback=0,valid=0;double normalized=0;
  bind(t,"CellID",&cid);bind(t,"sample",&sample);bind(t,"matched",&matched);bind(t,"sipm_um",&sipm);bind(t,"mip_fallback",&fallback);bind(t,"valid",&valid);bind(t,"normalized",&normalized);
  for(Long64_t i=0;i<t->GetEntries();++i){read(t,i);require(sample>=0&&sample<2&&matched>=0&&matched<2,"Invalid worker summary");require(records.emplace(std::make_tuple(cid,sample,matched),Record{sipm,fallback,valid,normalized}).second,"A channel is owned by multiple workers");}t->ResetBranchAddresses();
 }
 std::map<std::tuple<int,int,bool,int>,std::vector<double>> groups;
 for(const auto&[key,r]:records)if(r.valid)groups[{std::get<1>(key),r.sipm,bool(r.fallback),std::get<2>(key)}].push_back(r.normalized);
 auto group=output(out/"group_summary.tsv");group<<"sample\tSiPM_um\tMIP_adoption\tgeometry_condition\tvalid_channel_count\tmedian_temperature_MPV_over_MIP\n";
 for(const auto&[key,a]:groups)group<<(std::get<0>(key)?"mc":"data")<<'\t'<<std::get<1>(key)<<'\t'<<(std::get<2>(key)?"fallback":"individual")<<'\t'<<(std::get<3>(key)?"matched":"unmatched")<<'\t'<<a.size()<<'\t'<<median(a)<<'\n';
 std::map<std::tuple<int,bool,bool,int>,std::array<std::vector<double>,3>> pairedGroups;
 for(const auto&[key,d]:records){if(std::get<1>(key)!=0||!d.valid)continue;auto i=records.find({std::get<0>(key),1,std::get<2>(key)});require(i!=records.end(),"Missing partner channel summary");const auto&m=i->second;if(!m.valid)continue;auto& a=pairedGroups[{d.sipm,bool(d.fallback),bool(m.fallback),std::get<2>(key)}];a[0].push_back(d.normalized);a[1].push_back(m.normalized);a[2].push_back(d.normalized/m.normalized);}
 auto paired=output(out/"common_group_summary.tsv");paired<<"SiPM_um\tdata_MIP_adoption\tmc_MIP_adoption\tgeometry_condition\tcommon_valid_channel_count\tmedian_data_temperature_MPV_over_MIP\tmedian_mc_temperature_MPV_over_MIP\tmedian_channel_data_over_mc_ratio\n";
 for(const auto&[key,a]:pairedGroups)paired<<std::get<0>(key)<<'\t'<<(std::get<1>(key)?"fallback":"individual")<<'\t'<<(std::get<2>(key)?"fallback":"individual")<<'\t'<<(std::get<3>(key)?"matched":"unmatched")<<'\t'<<a[0].size()<<'\t'<<median(a[0])<<'\t'<<median(a[1])<<'\t'<<median(a[2])<<'\n';
 std::map<std::tuple<int,int,int,std::string>,long long> fitCounts;
 const std::vector<std::string> stages={"pedestal_adc","path_adc","temperature_adc","saved_hg_mip","detected_truth_mev","linear_prediction_adc"};
 for(const auto&part:parts){
  auto counts=readTable(part/"fit_status_counts.tsv");
  for(const auto&row:counts.rows){
   auto stage=value(counts,row,"stage");auto i=std::find(stages.begin(),stages.end(),stage);
   require(i!=stages.end(),"Unknown fit stage in worker summary");
   fitCounts[{value(counts,row,"sample")=="mc",value(counts,row,"geometry_condition")=="matched",int(i-stages.begin()),value(counts,row,"status")}]+=std::stoll(value(counts,row,"channels"));
  }
 }
 auto fc=output(out/"fit_status_counts.tsv");fc<<"sample\tgeometry_condition\tstage\tstatus\tchannels\n";for(const auto&[key,n]:fitCounts)fc<<(std::get<0>(key)?"mc":"data")<<'\t'<<(std::get<1>(key)?"matched":"unmatched")<<'\t'<<stages.at(std::get<2>(key))<<'\t'<<std::get<3>(key)<<'\t'<<n<<'\n';
 // Every channel object has exactly one owner. Only partial global maps add.
 {TFileMerger merger(false,false);merger.SetPrintLevel(0);require(merger.OutputFile((out/"spectra.root").c_str(),"CREATE"),"Cannot create merged ROOT");for(const auto&part:parts)require(merger.AddFile((part/"spectra.root").c_str(),false),"Cannot add worker ROOT");require(merger.Merge(),"ROOT worker merge failed");merger.CloseOutputFile();}
 {TFile merged((out/"spectra.root").c_str(),"UPDATE");require(!merged.IsZombie(),"Cannot update merged ROOT metadata");TParameter<double>("geometry_common_target",commonTarget).Write("geometry_common_target",TObject::kOverwrite);merged.Close();}
 // Parent-global counters appear in every worker; copy those exactly once.
 const auto& first=statuses.front();std::map<std::string,long long> counters;
 const std::set<std::string> global={"geometry_common_bins","crossings_read","generated_truth_valid","dependency_fits_enabled"};
 for(const auto&[column,index]:first.columns){if(column=="diagnostic_status")continue;long long total=0,initial=std::stoll(split(first.rows[0]).at(index));for(const auto&t:statuses){require(t.header==first.header,"Worker status headers differ");long long x=std::stoll(value(t,t.rows[0],column));if(global.count(column))require(x==initial,"Parent global counters differ across workers: "+column);else total+=x;}counters[column]=global.count(column)?initial:total;}
 auto status=output(out/"analysis_status.tsv");status<<first.header<<'\n';auto header=split(first.header);for(size_t i=0;i<header.size();++i){if(i)status<<'\t';if(header[i]=="diagnostic_status")status<<(counters["saved_HG_mismatches"]||counters["selected_energy_mismatches"]?"saved_calibration_mismatch":counters["common_matched_valid_channels"]?"paired_fits_available":"no_common_valid_matched_fit");else status<<counters.at(header[i]);}status<<'\n';
 for(const auto&part:parts){auto base=part/"comparison";if(!fs::exists(base))continue;for(const auto&entry:fs::recursive_directory_iterator(base)){if(!entry.is_regular_file()||entry.path().extension()!=".png"||entry.path().filename()=="muon_mip_channel_map.png")continue;auto dest=out/fs::relative(entry.path(),part);fs::create_directories(dest.parent_path());require(!fs::exists(dest),"Duplicate channel PNG: "+dest.string());fs::copy_file(entry.path(),dest);}}
 // Regenerate the single global map from the combined channel values.
 {auto file=openFile((out/"spectra.root").string());TCanvas canvas("muon_response","",mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.SetBatch(true);canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(2,2);for(int s=0;s<2;++s)for(int m=0;m<2;++m){canvas.cd(1+2*s+m);std::string sample=s?"mc":"data",condition=m?"matched":"unmatched";auto*h=dynamic_cast<TH2D*>(file->Get((sample+"_"+condition+"_mpv_map").c_str()));require(h,"Missing merged MPV map");h->SetMinimum(.5);h->SetMaximum(1.5);h->SetStats(false);h->SetTitle((sample+", "+condition+";Local channel (36 chip + channel);Layer;ADC MPV / MIP").c_str());h->Draw("COLZ");}fs::create_directories(out/"comparison");fastPng(canvas,(out/"comparison/muon_mip_channel_map.png").c_str());canvas.Clear();}
 require(bool(group)&&bool(paired)&&bool(fc)&&bool(status),"Final worker table write failed");std::cout<<"Merged "<<parts.size()<<" channel workers into "<<out<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<"MergeMuonShards: "<<e.what()<<'\n';return 1;}
