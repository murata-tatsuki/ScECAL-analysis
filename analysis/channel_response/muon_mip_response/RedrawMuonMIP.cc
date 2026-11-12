// Repaint existing muon PNGs using saved ROOT objects and fit metadata only.
#include "MuonMIPCommon.hh"
#include "MuonMIPPlotStyle.hh"
#include "../adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TList.h>
#include <TROOT.h>
#include <TStyle.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <tuple>

namespace redraw {
using mmr::require;
namespace fs=std::filesystem;
using FitKey=std::tuple<int,int,std::string,std::string>;
using DependencyKey=std::tuple<int,int,std::string,int>;
struct FitEntry {mmr_plot::FitInfo fit;double plateauDistance=INFINITY;};
struct Bin {int index;double low,high,mean;};
enum class Kind {Response,MPVProfiles,MeanProfiles,ChannelMap};
struct Item {int cell;std::string condition;Kind kind;fs::path png;};
struct Stamp {std::uintmax_t bytes;fs::file_time_type time;};
Stamp stamp(const fs::path& p){return {fs::file_size(p),fs::last_write_time(p)};}
bool same(const Stamp&a,const Stamp&b){return a.bytes==b.bytes&&a.time==b.time;}
struct Table {
 std::ifstream input;std::map<std::string,size_t> column;
 explicit Table(const fs::path& path):input(path){require(bool(input),"Cannot read "+path.string());std::string line;require(bool(std::getline(input,line)),"Empty table "+path.string());auto h=mmr::split(line);for(size_t i=0;i<h.size();++i)require(column.emplace(h[i],i).second,"Duplicate table column");}
 std::string value(const std::vector<std::string>& row,const std::string& key)const{auto it=column.find(key);require(it!=column.end()&&it->second<row.size(),"Missing table column "+key);return row[it->second];}
 bool read(std::vector<std::string>& row){std::string line;if(!std::getline(input,line)){require(input.eof()&&!input.bad(),"Failed reading saved fit table");return false;}row=mmr::split(line);require(row.size()==column.size(),"Incomplete table row");return true;}
};
int sample(const std::string& name){require(name=="data"||name=="mc","Invalid sample in fit metadata");return name=="mc";}
std::set<int> channels(const std::string& specification){std::set<int> out;if(specification=="all")return out;for(const auto& v:mmr::split(specification,',')){size_t n=0;int id=std::stoi(v,&n);require(n==v.size()&&mmr::connected(id),"Invalid physical channel "+v);out.insert(id);}require(!out.empty(),"Empty channel selection");return out;}
std::vector<Item> discover(const fs::path& base,const std::set<int>& selected){
 std::vector<Item> out;const fs::path comparison=base/"comparison";require(fs::is_directory(comparison),"Missing comparison directory");
 const std::regex response("^cell([0-9]+)_response[.]png$"),profiles("^cell([0-9]+)_profiles[.]png$"),means("^cell([0-9]+)_mean_profiles[.]png$");
 for(const auto& f:fs::recursive_directory_iterator(comparison)){
  if(!f.is_regular_file())continue;
  std::smatch m;const std::string name=f.path().filename().string();Kind kind;
  if(std::regex_match(name,m,response))kind=Kind::Response;
  else if(std::regex_match(name,m,profiles))kind=Kind::MPVProfiles;
  else if(std::regex_match(name,m,means))kind=Kind::MeanProfiles;
  else continue;
  int id=std::stoi(m[1].str());require(mmr::connected(id),"Invalid CellID in PNG name");if(!selected.empty()&&!selected.count(id))continue;
  auto relative=fs::relative(f.path(),comparison);auto p=relative.begin();require(p!=relative.end(),"Invalid comparison PNG path");std::string condition=p->string();
  if(condition!="matched"&&condition!="unmatched")continue;
  out.push_back({id,condition,kind,f.path()});
 }
 const fs::path map=comparison/"muon_mip_channel_map.png";
 if(selected.empty()&&fs::is_regular_file(map))out.push_back({-1,"",Kind::ChannelMap,map});
 std::sort(out.begin(),out.end(),[](const Item&a,const Item&b){return std::tie(a.cell,a.condition,a.kind,a.png)<std::tie(b.cell,b.condition,b.kind,b.png);});
 return out;
}
std::map<FitKey,FitEntry> loadFits(const fs::path& path,double plateau,const std::set<int>& wanted){
 Table table(path);std::vector<std::string> row;std::map<FitKey,FitEntry> out;
 while(table.read(row)){
  int id=std::stoi(table.value(row,"CellID"));if(!wanted.count(id))continue;
  const std::string condition=table.value(row,"geometry_condition");double distance=0;
  if(condition!="generated_independent"){
   double sigma=std::stod(table.value(row,"plateau_sigmas"));distance=std::abs(sigma-plateau);
   if(!std::isfinite(sigma)||distance>5e-10*std::max(1.,std::abs(plateau)))continue;
  }
  FitEntry entry;entry.plateauDistance=distance;auto& f=entry.fit;
  f.status=table.value(row,"status");f.valid=f.status=="ok";
  f.mpv=std::stod(table.value(row,"Landau_MPV"));f.mpvError=std::stod(table.value(row,"MPV_error"));f.landau=std::stod(table.value(row,"Landau_width"));f.gaussian=std::stod(table.value(row,"Gaussian_sigma"));f.chi2ndf=std::stod(table.value(row,"chi2_ndf"));
  const FitKey key{id,sample(table.value(row,"sample")),condition,table.value(row,"stage")};auto it=out.find(key);
  // The configured row precedes plateau variants. The closest numeric plateau
  // wins; exact ties retain that first row, including 12-digit rounding ties.
  if(it==out.end()||distance<it->second.plateauDistance)out[key]=entry;
 }
 return out;
}
std::map<DependencyKey,std::vector<Bin>> loadBins(const fs::path& path,const std::set<int>& wanted){
 Table table(path);std::vector<std::string> row;std::map<DependencyKey,std::vector<Bin>> out;
 const std::array<std::string,4> axes={"longitudinal","transverse","temperature","run"};
 while(table.read(row)){
  int id=std::stoi(table.value(row,"CellID"));if(!wanted.count(id)||table.value(row,"status")!="ok")continue;
  std::string axis=table.value(row,"axis");auto found=std::find(axes.begin(),axes.end(),axis);require(found!=axes.end(),"Unknown dependency axis");
  int s=sample(table.value(row,"sample")),a=std::distance(axes.begin(),found);std::string c=table.value(row,"geometry_condition");
  Bin b{std::stoi(table.value(row,"bin_index")),std::stod(table.value(row,"bin_low")),std::stod(table.value(row,"bin_high")),std::stod(table.value(row,"mean_axis_value"))};
  require(std::isfinite(b.low)&&std::isfinite(b.high)&&std::isfinite(b.mean)&&b.high>b.low,"Invalid dependency bin metadata");out[{id,s,c,a}].push_back(b);
 }
 for(auto& [key,bins]:out){std::sort(bins.begin(),bins.end(),[](const Bin&a,const Bin&b){return a.index<b.index;});for(size_t i=1;i<bins.size();++i)require(bins[i-1].index!=bins[i].index,"Duplicate valid dependency bin");}
 return out;
}
template<class T>std::unique_ptr<T> owned(TFile& file,const std::string& path){
 TObject* object=file.Get(path.c_str());if(!object)return {};
 auto* result=dynamic_cast<T*>(object);require(result,"Unexpected ROOT object type: "+path);
 if(auto* h=dynamic_cast<TH1*>(result))h->SetDirectory(nullptr);
 else{auto slash=path.rfind('/');auto* dir=slash==std::string::npos?static_cast<TDirectory*>(&file):file.GetDirectory(path.substr(0,slash).c_str());if(dir)dir->GetList()->Remove(result);}
 if(auto* f=dynamic_cast<TF1*>(result))f->AddToGlobalList(false);
 return std::unique_ptr<T>(result);
}
mmr_plot::FitInfo info(const std::map<FitKey,FitEntry>& fits,int id,int s,const std::string& condition,const std::string& stage,const TF1* curve){
 auto found=fits.find({id,s,condition,stage});mmr_plot::FitInfo out;
 if(found!=fits.end())out=found->second.fit;
 if(out.valid){require(curve&&curve->GetNpar()>=4,"Valid fit metadata has no saved TF1 for cell "+std::to_string(id)+"/"+stage);out.landau=curve->GetParameter(0);out.mpv=curve->GetParameter(1);out.gaussian=curve->GetParameter(3);}
 return out;
}
void response(TFile& root,TCanvas& canvas,const Item& item,const std::string& energy,const std::map<FitKey,FitEntry>& fits,const fs::path& png){
 const std::array<std::string,6> stages={"pedestal_adc","path_adc","temperature_adc","saved_hg_mip","detected_truth_mev","linear_prediction_adc"};
 const std::array<std::string,4> titles={"Pedestal subtracted HG","HG, 2 mm path normalization","HG, path and 20 C normalization","Saved HG energy, path normalized"};
 const std::array<std::string,4> units={"ADC","ADC","ADC","MIP"};
 std::array<std::array<std::unique_ptr<TH1D>,6>,2> hist;
 std::array<std::array<std::unique_ptr<TF1>,6>,2> curve;
 std::array<std::array<mmr_plot::FitInfo,6>,2> fitsInfo;
 for(int s=0;s<2;++s)for(int st=0;st<(s?6:4);++st){
  const std::string sample=s?"mc":"data",name=sample+"_"+item.condition+"_"+std::to_string(item.cell)+"_"+stages[st];
  const std::string path=sample+"/"+item.condition+"/"+energy+"/"+name;
  hist[s][st]=owned<TH1D>(root,path);
  if(st<4){curve[s][st]=owned<TF1>(root,path+"_langau");fitsInfo[s][st]=info(fits,item.cell,s,item.condition,stages[st],curve[s][st].get());}
 }
 const std::string generatedPath="mc/generated/"+energy+"/generated_truth_"+std::to_string(item.cell);
 auto generated=owned<TH1D>(root,generatedPath);
 canvas.Clear();canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(3,2);std::vector<std::unique_ptr<TH1D>> drawn;
 for(int st=0;st<4;++st)mmr_plot::drawResponsePanel(canvas,st+1,hist[0][st].get(),hist[1][st].get(),"Cell "+std::to_string(item.cell)+" | "+titles[st],units[st],drawn,&fitsInfo[0][st],&fitsInfo[1][st],curve[0][st].get(),curve[1][st].get());
 mmr_plot::drawResponsePanel(canvas,5,generated.get(),hist[1][4].get(),"MC generated vs reconstructed detected truth","MeV",drawn,nullptr,nullptr,nullptr,nullptr,"Generated truth","Detected crossings");
 canvas.cd(5);TLatex annotation;annotation.SetNDC();annotation.SetTextSize(.030);annotation.DrawLatex(.13,.92,"Independent truth uses primary strip chord");
 mmr_plot::drawResponsePanel(canvas,6,hist[1][5].get(),hist[1][2].get(),"MC linear prediction vs final HG","ADC",drawn,nullptr,nullptr,nullptr,nullptr,"Linear prediction","Final HG (20 C)");
 canvas.cd(6);annotation.DrawLatex(.13,.92,"Same reconstructed crossing cohort");
 fastPng(canvas,png.c_str());canvas.Clear();drawn.clear();
}
std::unique_ptr<TGraphAsymmErrors> displayGraph(TFile& root,int id,int s,const std::string& condition,const std::string& energy,int axis,const std::map<DependencyKey,std::vector<Bin>>& bins){
 const std::array<std::string,4> names={"longitudinal","transverse","temperature","run"};
 const std::string name=(s?"mc":"data")+std::string("_cell")+std::to_string(id)+"_"+condition+"_"+names[axis]+"_mpv";
 auto source=owned<TGraph>(root,"comparison/"+condition+"/"+energy+"/dependencies/"+name);
 auto found=bins.find({id,s,condition,axis});const size_t expected=found==bins.end()?0:found->second.size();
 if(!source){require(expected==0,"Valid dependency fits have no saved ROOT graph for cell "+std::to_string(id));return {};}
 require(static_cast<size_t>(source->GetN())==expected,"Saved ROOT graph and valid dependency table point counts differ for cell "+std::to_string(id));
 const auto* asymmetric=dynamic_cast<const TGraphAsymmErrors*>(source.get());const auto* symmetric=dynamic_cast<const TGraphErrors*>(source.get());
 require(asymmetric||symmetric,"Dependency graph is neither TGraphErrors nor TGraphAsymmErrors");
 auto out=std::make_unique<TGraphAsymmErrors>(source->GetN());out->SetName((name+"_display").c_str());out->SetTitle(source->GetTitle());
 for(int i=0;i<source->GetN();++i){
  double x=0,y=0;source->GetPoint(i,x,y);const auto& bin=found->second[i];double tolerance=1e-7+1e-10*std::max(1.,std::abs(x));
  require(std::isfinite(x)&&std::isfinite(y)&&std::abs(x-bin.mean)<=tolerance,"Saved ROOT point does not match dependency bin mean for cell "+std::to_string(id));
  double left=axis==3?0:x-bin.low,right=axis==3?0:bin.high-x;
  require(left>=-tolerance&&right>=-tolerance,"Saved graph point lies outside its dependency bin");left=std::max(0.,left);right=std::max(0.,right);
  double lowY=asymmetric?asymmetric->GetErrorYlow(i):symmetric->GetErrorY(i),highY=asymmetric?asymmetric->GetErrorYhigh(i):symmetric->GetErrorY(i);
  require(std::isfinite(lowY)&&std::isfinite(highY)&&lowY>=0&&highY>=0,"Invalid saved graph Y error");out->SetPoint(i,x,y);out->SetPointError(i,left,right,lowY,highY);
 }
 return out;
}
void profiles(TFile& root,TCanvas& canvas,const Item& item,const std::string& energy,const std::map<DependencyKey,std::vector<Bin>>& bins,const fs::path& png){
 std::array<std::array<std::unique_ptr<TGraphAsymmErrors>,2>,4> ownedGraphs;
 std::array<std::array<const TGraphAsymmErrors*,2>,4> pointers{};
 for(int axis=0;axis<4;++axis)for(int s=0;s<2;++s){ownedGraphs[axis][s]=displayGraph(root,item.cell,s,item.condition,energy,axis,bins);pointers[axis][s]=ownedGraphs[axis][s].get();}
 mmr_plot::drawMpvProfiles(canvas,item.cell,item.condition=="matched",pointers,png.string());canvas.Clear();
}
void meanProfiles(TFile& root,TCanvas& canvas,const Item& item,const std::string& energy,const fs::path& png){
 const std::array<std::string,4> axes={"Strip longitudinal position [mm]","Strip transverse position [mm]","Temperature [C]","Run"};
 std::array<std::array<std::unique_ptr<TProfile>,2>,4> profiles;
 canvas.Clear();canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(2,2);
 for(int axis=0;axis<4;++axis){
  canvas.cd(axis+1);TProfile* first=nullptr;
  for(int s=0;s<2;++s){
   const std::string sample=s?"mc":"data";
   const std::string name=sample+"_"+std::to_string(item.cell)+"_"+(item.condition=="matched"?"1":"0")+"_profile"+std::to_string(axis);
   profiles[axis][s]=owned<TProfile>(root,sample+"/"+item.condition+"/"+energy+"/"+name);
   auto* profile=profiles[axis][s].get();if(!profile||!profile->GetEntries())continue;
   profile->SetTitle(("Channel HG response;"+axes[axis]+";Path/temperature corrected ADC / current MIP").c_str());
   profile->SetStats(false);profile->SetLineColor(s?kRed+1:kBlack);profile->SetMarkerColor(s?kRed+1:kBlack);profile->SetMarkerStyle(s?24:20);
   profile->Draw(first?"SAME":"");if(!first)first=profile;
  }
  TLegend legend(.73,.73,.91,.87);
  for(int s=0;s<2;++s)if(profiles[axis][s]&&profiles[axis][s]->GetEntries())legend.AddEntry(profiles[axis][s].get(),s?"MC":"data","lp");
  legend.DrawClone();
 }
 fastPng(canvas,png.c_str());canvas.Clear();
}
void channelMap(TFile& root,TCanvas& canvas,const fs::path& png){
 std::array<std::array<std::unique_ptr<TH2D>,2>,2> maps;
 canvas.Clear();canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(2,2);
 for(int s=0;s<2;++s)for(int matched=0;matched<2;++matched){
  canvas.cd(1+2*s+matched);const std::string sample=s?"mc":"data",condition=matched?"matched":"unmatched";
  maps[s][matched]=owned<TH2D>(root,sample+"_"+condition+"_mpv_map");
  require(bool(maps[s][matched]),"Missing saved channel map: "+sample+"_"+condition+"_mpv_map");
  auto* h=maps[s][matched].get();h->SetMinimum(.5);h->SetMaximum(1.5);h->SetStats(false);
  h->SetTitle((sample+", "+condition+";Local channel (36 chip + channel);Layer;ADC MPV / MIP").c_str());h->Draw("COLZ");
 }
 fastPng(canvas,png.c_str());canvas.Clear();
}
void wait(pid_t child,int& status){while(true){pid_t r=waitpid(child,&status,0);if(r==child)return;if(r<0&&errno==EINTR)continue;throw std::runtime_error("Cannot wait for PNG redraw worker");}}
} // namespace redraw
int main(int argc,char** argv)try{
 using namespace redraw;
 require(argc>=2,"Usage: RedrawMuonMIP result_directory [--workers N] [--channels all|CellID,...]");
 const fs::path base=fs::canonical(argv[1]);int workers=4;std::string selection="all";
 for(int i=2;i<argc;++i){std::string option=argv[i];require(i+1<argc,"Missing redraw option value");if(option=="--workers")workers=std::stoi(argv[++i]);else if(option=="--channels")selection=argv[++i];else throw std::runtime_error("Unknown redraw option "+option);}
 require(workers>=1&&workers<=64,"Redraw workers must be between 1 and 64");
 const auto selected=channels(selection);auto items=discover(base,selected);
 if(items.empty()){std::cout<<"No existing muon PNG matches the requested channels"<<std::endl;return 0;}
 mmr::Config config((base/"inputs/configuration.tsv").string());const std::string energy=std::to_string(config.integer("energy",100))+"GeV";
 std::set<int> wanted;bool hasProfiles=false,hasResponse=false;
 for(const auto& item:items){
  if(item.kind!=Kind::ChannelMap)wanted.insert(item.cell);
  hasProfiles=hasProfiles||item.kind==Kind::MPVProfiles;
  hasResponse=hasResponse||item.kind==Kind::Response;
 }
 const std::array<fs::path,3> sources={base/"spectra.root",base/"fits.tsv",base/"dependency_fits.tsv"};std::map<fs::path,Stamp> original;
 original.emplace(sources[0],stamp(sources[0]));if(hasResponse)original.emplace(sources[1],stamp(sources[1]));if(hasProfiles)original.emplace(sources[2],stamp(sources[2]));
 const auto fits=hasResponse?loadFits(sources[1],config.number("plateau-sigmas",3),wanted):std::map<FitKey,FitEntry>{};
 const auto bins=hasProfiles?loadBins(sources[2],wanted):std::map<DependencyKey,std::vector<Bin>>{};
 for(const auto& entry:original)require(same(entry.second,stamp(entry.first)),"Scientific source changed while reading redraw metadata: "+entry.first.string());
 workers=std::min<int>(workers,std::max<size_t>(1,wanted.size()));std::map<int,int> owner;int channelIndex=0;for(int id:wanted)owner[id]=channelIndex++%workers;
 const fs::path binary=fs::read_symlink("/proc/self/exe").parent_path(),project=binary.parent_path().parent_path().parent_path();
 const fs::path scratchRoot=project/".scratch/muon_mip_response_redraw";fs::create_directories(scratchRoot);
 const std::string templatePath=(scratchRoot/("redraw_"+std::to_string(getpid())+"_XXXXXX")).string();std::vector<char> pattern(templatePath.begin(),templatePath.end());pattern.push_back('\0');
 char* created=mkdtemp(pattern.data());require(created,"Cannot create unique PNG redraw scratch directory");const fs::path scratch=created;
 std::cout<<"Repainting "<<items.size()<<" existing PNGs for "<<wanted.size()<<" channels with "<<workers<<" workers; scratch="<<scratch<<std::endl;
 auto render=[&](int worker){
  // ROOT input is opened after fork; no file handle or canvas is shared.
  gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetCanvasPreferGL(false);
  auto input=mmr::openFile(sources[0].string(),config.integer("allow-recovered",0)!=0);
  TCanvas canvas(("muon_redraw_"+std::to_string(worker)).c_str(),"",mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.SetBatch(true);canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);
  const fs::path directory=scratch/("worker_"+std::to_string(worker));fs::create_directories(directory);
  long long done=0;
  for(const auto& item:items){
   if(item.kind==Kind::ChannelMap?worker!=0:owner.at(item.cell)!=worker)continue;
   require(fs::is_regular_file(item.png),"Existing PNG disappeared during redraw: "+item.png.string());
   const std::array<std::string,4> names={"response.png","profiles.png","mean_profiles.png","channel_map.png"};
   const fs::path temporary=directory/names[static_cast<size_t>(item.kind)];
   switch(item.kind){
    case Kind::Response:response(*input,canvas,item,energy,fits,temporary);break;
    case Kind::MPVProfiles:profiles(*input,canvas,item,energy,bins,temporary);break;
    case Kind::MeanProfiles:meanProfiles(*input,canvas,item,energy,temporary);break;
    case Kind::ChannelMap:channelMap(*input,canvas,temporary);break;
   }
   require(fs::is_regular_file(temporary)&&fs::file_size(temporary)>0,"PNG drawing produced no output");
   // Same project filesystem: rename replaces only this worker's assigned PNG.
   fs::rename(temporary,item.png);++done;
   if(done%100==0)std::cout<<"Redraw worker "<<worker<<": "<<done<<" PNGs"<<std::endl;
  }
  input->Close();std::cout<<"Redraw worker "<<worker<<" complete: "<<done<<" PNGs"<<std::endl;
 };
 if(workers==1){try{render(0);}catch(const std::exception& e){throw std::runtime_error(std::string(e.what())+"; temporary redraw files preserved in "+scratch.string());}}
 else{
  std::vector<pid_t> children;
  for(int worker=0;worker<workers;++worker){std::cout.flush();std::cerr.flush();pid_t pid=fork();if(pid<0){for(pid_t child:children){int status=0;wait(child,status);}throw std::runtime_error("Cannot fork PNG redraw worker; temporary files preserved in "+scratch.string());}
   if(pid==0){try{render(worker);std::cout.flush();std::cerr.flush();_exit(0);}catch(const std::exception& e){std::cerr<<"Redraw worker "<<worker<<": "<<e.what()<<"; scratch="<<scratch<<std::endl;_exit(1);}}
   children.push_back(pid);
  }
  bool failed=false;for(pid_t child:children){int status=0;wait(child,status);failed=failed||!WIFEXITED(status)||WEXITSTATUS(status)!=0;}
  require(!failed,"PNG redraw worker failed; temporary files preserved in "+scratch.string());
 }
 for(const auto& entry:original)require(same(entry.second,stamp(entry.first)),"Scientific source changed during PNG redraw: "+entry.first.string()+"; scratch="+scratch.string());
 fs::remove_all(scratch);std::error_code cleanup;fs::remove(scratchRoot,cleanup);
 std::cout<<"PNG redraw complete; saved ROOT and fit tables were read-only"<<std::endl;return 0;
}catch(const std::exception& e){std::cerr<<"RedrawMuonMIP: "<<e.what()<<std::endl;return 1;}
