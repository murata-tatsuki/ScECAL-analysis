// Event-level diagnostics on the fixed selected cohort saved by IntegratedReview.
// No raw/calibration input is read and no calibration or event selection is changed.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include "../../analysis/channel_response/adc_energy/FastPng.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace fs=std::filesystem;
namespace {
constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
void need(bool ok,const std::string& msg){if(!ok)throw std::runtime_error(msg);}
std::vector<std::string> fields(const std::string& line){std::vector<std::string> result;std::istringstream input(line);std::string item;while(std::getline(input,item,'\t'))result.push_back(item);return result;}
struct Event {
 int sample=0,file=0,run=0;long long entry=0;bool cog=false,tail=false,window=false;
 double x=0,y=0,e=0,central=0,n=0,ncentral=0,noT=0,back=0,edge=0,depth=0,temperature=NaN,lg=NaN;
 std::array<double,30> layers{};bool hasLayers=false;
 double perHit()const{return 1000*e/n;}
 double centerFraction()const{return central/e;}
 double tempProxy()const{return std::log(e/noT);}
};
using Samples=std::array<std::vector<Event>,2>;
Samples readEvents(const fs::path& path){
 std::ifstream input(path);need(bool(input),"Cannot read "+path.string());std::string line;need(bool(std::getline(input,line)),"Empty events table");
 const auto names=fields(line);std::map<std::string,size_t> index;for(size_t i=0;i<names.size();++i)index.emplace(names[i],i);
 for(const std::string key:{"sample","file","entry","run","cog20","tail_rich","energy_window","x_mm","y_mm","energy_GeV","central_GeV","hits","central_hits","no_temperature_GeV","back5_fraction","edge_fraction","depth_layer"})need(index.count(key),"Missing column "+key);
 Samples result;std::set<std::tuple<int,int,long long>> identities;
 while(std::getline(input,line)){
  auto f=fields(line);need(f.size()==names.size(),"Malformed events row "+path.string());
  auto text=[&](const std::string& name)->const std::string&{return f.at(index.at(name));};
  auto number=[&](const std::string& name){return std::stod(text(name));};
  Event e;e.sample=text("sample")=="mc"?1:0;need(text("sample")=="data"||text("sample")=="mc","Invalid sample");
  e.file=std::stoi(text("file"));e.entry=std::stoll(text("entry"));e.run=std::stoi(text("run"));
  need(identities.emplace(e.sample,e.file,e.entry).second,"Duplicate sample/file/entry");
  e.cog=std::stoi(text("cog20"));e.tail=std::stoi(text("tail_rich"));e.window=std::stoi(text("energy_window"));
  e.x=number("x_mm");e.y=number("y_mm");e.e=number("energy_GeV");e.central=number("central_GeV");e.n=number("hits");e.ncentral=number("central_hits");
  e.noT=number("no_temperature_GeV");e.back=number("back5_fraction");e.edge=number("edge_fraction");e.depth=number("depth_layer");
  if(index.count("mean_temperature_C"))e.temperature=number("mean_temperature_C");
  if(index.count("temperature_C"))e.temperature=number("temperature_C");
  if(index.count("lg_energy_GeV"))e.lg=number("lg_energy_GeV");
  e.hasLayers=index.count("layer0_GeV")>0;
  if(e.hasLayers)for(int l=0;l<30;++l)e.layers[l]=number("layer"+std::to_string(l)+"_GeV");
  result[e.sample].push_back(e);
 }
 return result;
}
struct Stats {
 double sw=0,sw2=0,sx=0,sxx=0;long long count=0;
 void add(double x,double w=1){if(!std::isfinite(x)||!(w>0))return;++count;sw+=w;sw2+=w*w;sx+=w*x;sxx+=w*x*x;}
 double mean()const{return sw?sx/sw:NaN;}
 double variance()const{return sw?std::max(0.,sxx/sw-mean()*mean()):NaN;}
 double rms()const{return std::sqrt(variance());}
 double resolution()const{return rms()/mean();}
 double neff()const{return sw2?sw*sw/sw2:0;}
};
double quantile(std::vector<std::pair<double,double>> v,double p){
 v.erase(std::remove_if(v.begin(),v.end(),[](auto q){return !std::isfinite(q.first)||!(q.second>0);}),v.end());
 if(v.empty())return NaN;std::sort(v.begin(),v.end());double total=0;for(auto q:v)total+=q.second;
 double running=0;for(size_t i=0;i<v.size();++i){double before=running;running+=v[i].second;if(running>=p*total){if(i==0)return v[i].first;double f=std::clamp((p*total-before)/v[i].second,0.,1.);return v[i-1].first+f*(v[i].first-v[i-1].first);}}
 return v.back().first;
}
struct Distribution {
 Stats st;double q16=NaN,q50=NaN,q84=NaN,q95=NaN,q99=NaN,core=NaN;
 Distribution(const std::vector<double>& v,const std::vector<double>& w={}){
  std::vector<std::pair<double,double>> q;for(size_t i=0;i<v.size();++i){double weight=w.empty()?1:w.at(i);st.add(v[i],weight);q.emplace_back(v[i],weight);}
  q16=quantile(q,.16);q50=quantile(q,.5);q84=quantile(q,.84);q95=quantile(q,.95);q99=quantile(q,.99);core=.5*(q84-q16)/q50;
 }
};
double correlation(const std::vector<double>& x,const std::vector<double>& y){
 need(x.size()==y.size(),"Correlation size mismatch");Stats a,b;double cross=0;long long n=0;for(size_t i=0;i<x.size();++i)if(std::isfinite(x[i])&&std::isfinite(y[i])){a.add(x[i]);b.add(y[i]);cross+=x[i]*y[i];++n;}
 return a.rms()>0&&b.rms()>0?(cross/n-a.mean()*b.mean())/(a.rms()*b.rms()):NaN;
}
std::vector<double> values(const std::vector<Event>& e,const std::string& field){
 std::vector<double> v;v.reserve(e.size());for(const auto& q:e){
  double x=field=="energy"?q.e:field=="hits"?q.n:field=="central_hits"?q.ncentral:field=="peripheral_hits"?q.n-q.ncentral:field=="energy_per_hit_MeV"?q.perHit():field=="central_fraction"?q.centerFraction():field=="back5_fraction"?q.back:field=="edge_fraction"?q.edge:field=="depth_layer"?q.depth:field=="temperature_C"?q.temperature:field=="lg_fraction"?q.lg/q.e:field=="temperature_correction_log"?q.tempProxy():field=="x_mm"?q.x:q.y;
  v.push_back(x);
 }return v;
}
std::string sampleName(int s){return s?"mc":"data";}
std::ofstream table(const fs::path& path){std::ofstream out(path);need(bool(out),"Cannot write "+path.string());out<<std::setprecision(12);return out;}
uint64_t mix(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
struct FoldSet {std::vector<int> fold;std::string grouping;int count=5,groups=0;};
FoldSet folds(const std::vector<Event>& ev,bool grouped){
 FoldSet f;std::set<int> runs,files;for(auto e:ev){runs.insert(e.run);files.insert(e.file);}
 std::map<int,int> assigned;
 if(grouped){const auto& keys=runs.size()>=3?runs:files;f.grouping=runs.size()>=3?"run":"input_file";f.groups=keys.size();f.count=std::min(5,f.groups);need(f.count>=2,"Too few independent run/file groups");int i=0;for(int k:keys)assigned[k]=i++%f.count;}
 else{f.grouping="event_hash";f.groups=ev.size();}
 for(const auto& e:ev)f.fold.push_back(grouped?assigned.at(f.grouping=="run"?e.run:e.file):mix(uint64_t(e.entry)^(uint64_t(e.file+1)*0x9e3779b97f4a7c15ULL))%5);
 return f;
}
std::vector<double> features(const Event& e,const std::string& model){
 std::vector<double> f={1};
 if(model.find("position")!=std::string::npos){double x=e.x/20,y=e.y/20;for(double v:{x,y,x*x,y*y,x*y})f.push_back(v);}
 if(model.find("nhit")!=std::string::npos){double n=std::log(e.n/300.);f.push_back(n);f.push_back(n*n);}
 if(model.find("leak")!=std::string::npos){for(double v:{e.back,e.edge,e.depth/30,e.centerFraction()})f.push_back(v);}
 if(model=="temperature_proxy")f.push_back(e.tempProxy());
 if(model=="temperature"){double t=(e.temperature-20)/5;f.push_back(t);f.push_back(t*t);}
 return f;
}
std::vector<double> solve(std::vector<std::vector<double>> a,std::vector<double> b){
 int n=b.size();for(int i=0;i<n;++i){int pivot=i;for(int j=i+1;j<n;++j)if(std::abs(a[j][i])>std::abs(a[pivot][i]))pivot=j;need(std::abs(a[pivot][i])>1e-14,"Singular CV model");std::swap(a[i],a[pivot]);std::swap(b[i],b[pivot]);double d=a[i][i];for(int k=i;k<n;++k)a[i][k]/=d;b[i]/=d;for(int j=0;j<n;++j)if(j!=i){double f=a[j][i];for(int k=i;k<n;++k)a[j][k]-=f*a[i][k];b[j]-=f*b[i];}}
 return b;
}
struct CV {std::vector<double> energy;std::vector<std::vector<double>> coefficients;double maxFactor=0;};
CV corrected(const std::vector<Event>& ev,const FoldSet& partition,const std::string& model){
 const int p=features(ev.front(),model).size();CV result;result.energy.resize(ev.size());
 for(int fold=0;fold<partition.count;++fold){
  std::vector<std::vector<double>> a(p,std::vector<double>(p));std::vector<double>b(p);double trainLog=0;long long n=0;
  for(size_t i=0;i<ev.size();++i)if(partition.fold[i]!=fold){auto f=features(ev[i],model);double y=std::log(ev[i].e);trainLog+=y;++n;for(int j=0;j<p;++j){b[j]+=f[j]*y;for(int k=0;k<p;++k)a[j][k]+=f[j]*f[k];}}
  need(n>20*p,"Too few training events");trainLog/=n;
  for(int j=1;j<p;++j)a[j][j]+=1e-7*n;
  auto coefficient=solve(a,b);result.coefficients.push_back(coefficient);
  for(size_t i=0;i<ev.size();++i)if(partition.fold[i]==fold){auto f=features(ev[i],model);double prediction=std::inner_product(f.begin(),f.end(),coefficient.begin(),0.);double factor=std::exp(trainLog-prediction);need(std::isfinite(factor)&&factor>0,"Invalid CV correction");result.maxFactor=std::max(result.maxFactor,std::max(factor,1/factor));result.energy[i]=ev[i].e*factor;}
 }
 return result;
}
using Key=std::array<int,3>;
Key binKey(const Event&e,int beam,const std::string& mode){
 const double nwidth=beam==5?10:beam==20?20:40;
 return {mode.find("position")!=std::string::npos?int(std::floor(e.x/5)):0,mode.find("position")!=std::string::npos?int(std::floor(e.y/5)):0,mode.find("nhit")!=std::string::npos?int(std::floor(e.n/nwidth)):0};
}
struct Reweighted {std::array<std::vector<double>,2> weights;std::array<long long,2> covered{};std::array<double,2> neff{};};
Reweighted commonWeights(const Samples& ev,const std::array<FoldSet,2>& partitions,int beam,const std::string& mode){
 Reweighted result;for(int s=0;s<2;++s)result.weights[s].resize(ev[s].size());
 const int count=std::min(partitions[0].count,partitions[1].count);
 // Use the same fold modulus on both samples when independent group counts differ.
 for(int fold=0;fold<count;++fold){
  std::map<Key,std::array<long long,2>> train;
  for(int s=0;s<2;++s)for(size_t i=0;i<ev[s].size();++i)if(partitions[s].fold[i]%count!=fold)++train[binKey(ev[s][i],beam,mode)][s];
  for(int s=0;s<2;++s)for(size_t i=0;i<ev[s].size();++i)if(partitions[s].fold[i]%count==fold){
   auto it=train.find(binKey(ev[s][i],beam,mode));if(it==train.end())continue;const auto&n=it->second;if(n[0]<5||n[1]<5)continue;
   result.weights[s][i]=mode.rfind("mc_to_data_",0)==0?(s?double(n[0])/n[1]:1.):double(std::min(n[0],n[1]))/n[s];++result.covered[s];
  }
 }
 for(int s=0;s<2;++s){Stats st;for(double w:result.weights[s])st.add(0,w);result.neff[s]=st.neff();}return result;
}
void text(double y,const std::string& label,int color=kBlack){TLatex t;t.SetNDC();t.SetTextSize(.033);t.SetTextColor(color);t.DrawLatex(.13,y,label.c_str());}
std::string number(double x,int p=3){std::ostringstream out;out<<std::fixed<<std::setprecision(p)<<x;return out.str();}
std::pair<double,double> limits(const std::array<std::vector<double>,2>& v,bool positive=false){
 std::vector<std::pair<double,double>> pooled;for(const auto& sample:v)for(double x:sample)if(std::isfinite(x))pooled.emplace_back(x,1);
 double lo=quantile(pooled,.001),hi=quantile(pooled,.999);need(std::isfinite(lo)&&std::isfinite(hi),"Invalid plot range");
 double pad=std::max(1e-4,.08*(hi-lo));return {positive?std::max(0.,lo-pad):lo-pad,hi+pad};
}
void overlay(TCanvas& c,int pad,const std::array<std::vector<double>,2>& v,const std::string& title,const std::string& axis,TFile& root,const std::string& name,const std::array<std::vector<double>,2>* weights=nullptr,bool positive=false,bool integerBins=false){
 auto* plotPad=c.cd(pad);plotPad->SetLeftMargin(.14);plotPad->SetBottomMargin(.12);plotPad->SetRightMargin(.06);auto range=limits(v,positive);int bins=90;if(integerBins){double step=std::max(1.,std::ceil((range.second-range.first)/90.));range.first=std::floor(range.first/step)*step-.5;range.second=std::ceil((range.second+.5)/step)*step-.5;bins=std::lround((range.second-range.first)/step);}
 std::array<TH1D,2> h={TH1D((name+"_data").c_str(),"",bins,range.first,range.second),TH1D((name+"_mc").c_str(),"",bins,range.first,range.second)};double ymax=0;
 for(int s=0;s<2;++s){h[s].Sumw2();for(size_t i=0;i<v[s].size();++i)h[s].Fill(v[s][i],weights?weights->at(s).at(i):1);double total=h[s].Integral(0,h[s].GetNbinsX()+1);if(total>0)h[s].Scale(1/total);h[s].SetLineColor(s?kRed+1:kBlack);h[s].SetLineWidth(2);h[s].SetStats(false);ymax=std::max(ymax,h[s].GetMaximum());root.cd();h[s].Write();}
 h[0].SetTitle((title+";"+axis+";Probability / bin").c_str());h[0].SetMaximum(1.4*ymax);auto*d=h[0].DrawCopy("HIST");auto*m=h[1].DrawCopy("HIST SAME");TLegend legend(.65,.74,.89,.87);legend.SetBorderSize(0);legend.AddEntry(d,"data, LG tail excluded","l");legend.AddEntry(m,"MC e_ssa threshold","l");legend.DrawClone();
}
void profile(TCanvas& c,int pad,const Samples& events,const std::array<std::vector<double>,2>& energy,const std::string& variable,const std::string& axis,TFile& root,const std::string& name){
 std::array<std::vector<double>,2> x={values(events[0],variable),values(events[1],variable)};auto range=limits(x);
 auto* plotPad=c.cd(pad);plotPad->SetLeftMargin(.14);plotPad->SetBottomMargin(.12);plotPad->SetRightMargin(.06);std::array<TProfile,2>h={TProfile((name+"_data").c_str(),"",24,range.first,range.second),TProfile((name+"_mc").c_str(),"",24,range.first,range.second)};
 double ymin=1,ymax=1;
 for(int s=0;s<2;++s){Distribution d(energy[s]);for(size_t i=0;i<x[s].size();++i)if(std::isfinite(x[s][i]))h[s].Fill(x[s][i],energy[s][i]/d.st.mean());h[s].SetLineColor(s?kRed+1:kBlack);h[s].SetMarkerColor(s?kRed+1:kBlack);h[s].SetMarkerStyle(s?24:20);h[s].SetStats(false);root.cd();h[s].Write();for(int b=1;b<=h[s].GetNbinsX();++b)if(h[s].GetBinEntries(b)>=5){ymin=std::min(ymin,h[s].GetBinContent(b)-h[s].GetBinError(b));ymax=std::max(ymax,h[s].GetBinContent(b)+h[s].GetBinError(b));}}
 h[0].SetTitle(("Event energy dependence;"+axis+";Mean E / sample mean E").c_str());h[0].SetMinimum(ymin-.10*(ymax-ymin+.02));h[0].SetMaximum(ymax+.35*(ymax-ymin+.02));auto*d=h[0].DrawCopy();auto*m=h[1].DrawCopy("SAME");TLegend legend(.73,.74,.89,.87);legend.SetBorderSize(0);legend.AddEntry(d,"data","lp");legend.AddEntry(m,"MC","lp");legend.DrawClone();
}
}
int main(int argc,char**argv)try{
 need(argc==3,"Usage: EventDistributionReview integrated_directory output_directory");
 fs::path input=fs::canonical(argv[1]),out=argv[2];fs::create_directories(out);
 gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);
 auto summary=table(out/"distribution_summary.tsv"),decomp=table(out/"variance_components.tsv"),correlations=table(out/"correlations.tsv"),impact=table(out/"crossvalidated_impact.tsv"),weights=table(out/"reweighting_impact.tsv"),models=table(out/"models.tsv"),metadata=table(out/"input_events.tsv");
 summary<<"energy_GeV\tselection\tsample\tobservable\tevents\teffective_events\tmean\tRMS\tRMS_over_mean\tq16\tq50\tq84\tq95\tq99\tsigma68_over_median\n";
 decomp<<"energy_GeV\tselection\tsample\tevents\tvariance_log_energy\tvariance_log_hits\tvariance_log_energy_per_hit\ttwice_log_covariance\texact_log_closure_error\tRMS_over_mean\tcandidate_common_gain_sigma\tcandidate_independent_existing_hit_noise_sigma_MeV\tdiagonal_layer_normalized_variance\toffdiagonal_layer_normalized_covariance\n";
 correlations<<"energy_GeV\tselection\tsample\tobservable\tPearson_energy\tPearson_log_energy\tPearson_hits\tPearson_energy_per_hit\n";
 impact<<"energy_GeV\tselection\tfold_scheme\tdata_grouping\tmc_grouping\tmodel\tdata_events\tmc_events\tdata_RMS\tmc_RMS\tvariance_gap_reduction\tdata_sigma68\tmc_sigma68\tdata_mean_shift\tmc_mean_shift\tdata_max_correction_factor\tmc_max_correction_factor\n";
 weights<<"energy_GeV\tselection\tfold_scheme\tmatching_variables\tdata_covered_events\tmc_covered_events\tdata_coverage\tmc_coverage\tdata_effective_events\tmc_effective_events\tdata_covered_baseline_RMS\tmc_covered_baseline_RMS\tdata_weighted_RMS\tmc_weighted_RMS\tcovered_variance_gap_reduction\tfull_cohort_variance_gap_reduction_diagnostic\tdata_weighted_mean_hits\tmc_weighted_mean_hits\tdata_weighted_sigma68\tmc_weighted_sigma68\tdata_max_weight\tmc_max_weight\tdiagnostic_status\n";
 models<<"energy_GeV\tselection\tsample\tfold_scheme\tmodel\tfold\tcoefficient_index\tcoefficient\n";
 metadata<<"path\tbytes\tfilesystem_mtime_ticks\tdata_rows\tmc_rows\n";
 std::map<std::tuple<int,int,std::string>,std::array<double,2>> layerVariance;
 {
  std::ifstream cov(input/"covariance.tsv");std::string line;std::getline(cov,line);while(std::getline(cov,line)){auto f=fields(line);if(f.size()<6)continue;layerVariance[{std::stoi(f[0]),f[1]=="mc"?1:0,f[2]}]={std::stod(f[3]),std::stod(f[4])};}
 }
 TFile root((out/"event_components.root").c_str(),"RECREATE");need(!root.IsZombie(),"Cannot create ROOT output");
 for(int beam:{5,20,100}){
  fs::path eventPath=input/"comparison/20mm"/(std::to_string(beam)+"GeV")/"events.tsv";
  const auto before=std::make_pair(fs::file_size(eventPath),fs::last_write_time(eventPath));auto all=readEvents(eventPath);
  metadata<<eventPath.string()<<'\t'<<before.first<<'\t'<<before.second.time_since_epoch().count()<<'\t'<<all[0].size()<<'\t'<<all[1].size()<<'\n';
  for(const std::string selection:{"cog20","cog20_energy_window"}){
   Samples ev;for(int s=0;s<2;++s)for(const auto&e:all[s])if(e.cog&&!e.tail&&(selection=="cog20"||e.window)){need(e.e>0&&e.n>0&&e.noT>0&&std::isfinite(e.e)&&std::isfinite(e.n),"Invalid selected event");ev[s].push_back(e);}
   need(ev[0].size()>200&&ev[1].size()>200,"Insufficient selected events");
   fs::path dest=out/"comparison"/selection/(std::to_string(beam)+"GeV");fs::create_directories(dest);
   const std::string prefix=std::to_string(beam)+"_"+selection+"_";
   std::array<std::vector<double>,2> original={values(ev[0],"energy"),values(ev[1],"energy")};
   std::array<Distribution,2> baseline={Distribution(original[0]),Distribution(original[1])};
   double gap=std::pow(baseline[0].st.resolution(),2)-std::pow(baseline[1].st.resolution(),2);
   for(int s=0;s<2;++s){
    const std::vector<std::string> observables={"energy","hits","central_hits","peripheral_hits","energy_per_hit_MeV","central_fraction","back5_fraction","edge_fraction","depth_layer","temperature_correction_log"};
    for(const auto& field:observables){Distribution d(values(ev[s],field));summary<<beam<<'\t'<<selection<<'\t'<<sampleName(s)<<'\t'<<field<<'\t'<<d.st.count<<'\t'<<d.st.neff()<<'\t'<<d.st.mean()<<'\t'<<d.st.rms()<<'\t'<<d.st.resolution()<<'\t'<<d.q16<<'\t'<<d.q50<<'\t'<<d.q84<<'\t'<<d.q95<<'\t'<<d.q99<<'\t'<<d.core<<'\n';}
    auto e=original[s],n=values(ev[s],"hits"),r=values(ev[s],"energy_per_hit_MeV");std::vector<double>le,ln,lr;double cross=0;Stats se,sn,sr;
    for(size_t i=0;i<e.size();++i){le.push_back(std::log(e[i]));ln.push_back(std::log(n[i]));lr.push_back(std::log(e[i]/n[i]));se.add(le.back());sn.add(ln.back());sr.add(lr.back());cross+=ln.back()*lr.back();}
    double covariance=cross/e.size()-sn.mean()*sr.mean(),closure=se.variance()-sn.variance()-sr.variance()-2*covariance;
    need(std::abs(closure)<1e-9,"Log variance identity failed");
    double sg=gap>0?std::sqrt(gap/(1+std::pow(baseline[1].st.resolution(),2))):NaN;
    double sh=gap>0?1000*baseline[1].st.mean()*std::sqrt(gap/Distribution(values(ev[1],"hits")).st.mean()):NaN;
    auto lv=layerVariance.find({beam,s,selection});
    if(ev[s].front().hasLayers){
     std::array<Stats,30> layers;double crosses[30][30]{};
     for(const auto& event:ev[s]){double layerSum=std::accumulate(event.layers.begin(),event.layers.end(),0.);need(std::abs(layerSum-event.e)<1e-8,"Layer energy closure failed");for(int l=0;l<30;++l){layers[l].add(event.layers[l]);for(int k=0;k<30;++k)crosses[l][k]+=event.layers[l]*event.layers[k];}}
     double diagonal=0,total=0;for(int l=0;l<30;++l)for(int k=0;k<30;++k){double covariance=crosses[l][k]/ev[s].size()-layers[l].mean()*layers[k].mean();total+=covariance;if(l==k)diagonal+=covariance;}
     need(std::abs(total-baseline[s].st.variance())<1e-8,"Layer covariance closure failed");
     layerVariance[{beam,s,selection}]={diagonal/std::pow(baseline[s].st.mean(),2),(total-diagonal)/std::pow(baseline[s].st.mean(),2)};lv=layerVariance.find({beam,s,selection});
     TH2D covarianceMap((prefix+sampleName(s)+"_layer_covariance").c_str(),"Layer energy covariance;Layer;Layer",30,-.5,29.5,30,-.5,29.5);for(int l=0;l<30;++l)for(int k=0;k<30;++k)covarianceMap.SetBinContent(l+1,k+1,(crosses[l][k]/ev[s].size()-layers[l].mean()*layers[k].mean())/std::pow(baseline[s].st.mean(),2));root.cd();covarianceMap.Write();
    }
    decomp<<beam<<'\t'<<selection<<'\t'<<sampleName(s)<<'\t'<<e.size()<<'\t'<<se.variance()<<'\t'<<sn.variance()<<'\t'<<sr.variance()<<'\t'<<2*covariance<<'\t'<<closure<<'\t'<<baseline[s].st.resolution()<<'\t'<<sg<<'\t'<<sh<<'\t'<<(lv==layerVariance.end()?NaN:lv->second[0])<<'\t'<<(lv==layerVariance.end()?NaN:lv->second[1])<<'\n';
    for(const std::string field:{"hits","energy_per_hit_MeV","central_fraction","back5_fraction","edge_fraction","depth_layer","temperature_C","lg_fraction","temperature_correction_log","x_mm","y_mm"}){
     auto v=values(ev[s],field);correlations<<beam<<'\t'<<selection<<'\t'<<sampleName(s)<<'\t'<<field<<'\t'<<correlation(e,v)<<'\t'<<correlation(le,v)<<'\t'<<correlation(n,v)<<'\t'<<correlation(r,v)<<'\n';
    }
   }
   std::map<std::string,Reweighted> weightedExamples;
   std::map<std::string,std::array<std::vector<double>,2>> correctedValues;correctedValues["baseline"]=original;
   for(bool grouped:{false,true}){
    std::array<FoldSet,2> partitions={folds(ev[0],grouped),folds(ev[1],grouped)};std::string scheme=grouped?"run_or_file_held_out":"event_held_out";
    for(const std::string model:{"nhit","position","position_nhit","position_leak","position_nhit_leak","temperature_proxy","temperature"}){
     if(model=="temperature"&&(!std::isfinite(ev[0].front().temperature)||!std::isfinite(ev[1].front().temperature)))continue;
     std::array<CV,2> cv={corrected(ev[0],partitions[0],model),corrected(ev[1],partitions[1],model)};
     std::array<Distribution,2> d={Distribution(cv[0].energy),Distribution(cv[1].energy)};
     impact<<beam<<'\t'<<selection<<'\t'<<scheme<<'\t'<<partitions[0].grouping<<'\t'<<partitions[1].grouping<<'\t'<<model<<'\t'<<ev[0].size()<<'\t'<<ev[1].size()<<'\t'<<d[0].st.resolution()<<'\t'<<d[1].st.resolution()<<'\t'<<(gap?1-(std::pow(d[0].st.resolution(),2)-std::pow(d[1].st.resolution(),2))/gap:NaN)<<'\t'<<d[0].core<<'\t'<<d[1].core<<'\t'<<d[0].st.mean()/baseline[0].st.mean()-1<<'\t'<<d[1].st.mean()/baseline[1].st.mean()-1<<'\t'<<cv[0].maxFactor<<'\t'<<cv[1].maxFactor<<'\n';
     for(int s=0;s<2;++s)for(size_t fold=0;fold<cv[s].coefficients.size();++fold)for(size_t j=0;j<cv[s].coefficients[fold].size();++j)models<<beam<<'\t'<<selection<<'\t'<<sampleName(s)<<'\t'<<scheme<<'\t'<<model<<'\t'<<fold<<'\t'<<j<<'\t'<<cv[s].coefficients[fold][j]<<'\n';
     if(!grouped)correctedValues[model]={std::move(cv[0].energy),std::move(cv[1].energy)};
    }
    for(const std::string mode:{"nhit","position","position_nhit","mc_to_data_nhit","mc_to_data_position","mc_to_data_position_nhit"}){
     auto rw=commonWeights(ev,partitions,beam,mode);if(!grouped&&mode.rfind("mc_to_data_",0)==0)weightedExamples[mode]=rw;std::array<std::vector<double>,2> maskWeights=rw.weights;
     for(auto& sample:maskWeights)for(double&w:sample)w=w>0?1:0;
     std::array<Distribution,2> covered={Distribution(original[0],maskWeights[0]),Distribution(original[1],maskWeights[1])};
     std::array<Distribution,2> d={Distribution(original[0],rw.weights[0]),Distribution(original[1],rw.weights[1])};
     std::array<Distribution,2> nh={Distribution(values(ev[0],"hits"),rw.weights[0]),Distribution(values(ev[1],"hits"),rw.weights[1])};
     double coveredGap=std::pow(covered[0].st.resolution(),2)-std::pow(covered[1].st.resolution(),2),newGap=std::pow(d[0].st.resolution(),2)-std::pow(d[1].st.resolution(),2);
     weights<<beam<<'\t'<<selection<<'\t'<<scheme<<'\t'<<mode<<'\t'<<rw.covered[0]<<'\t'<<rw.covered[1]<<'\t'<<double(rw.covered[0])/ev[0].size()<<'\t'<<double(rw.covered[1])/ev[1].size()<<'\t'<<rw.neff[0]<<'\t'<<rw.neff[1]<<'\t'<<covered[0].st.resolution()<<'\t'<<covered[1].st.resolution()<<'\t'<<d[0].st.resolution()<<'\t'<<d[1].st.resolution()<<'\t'<<(coveredGap?1-newGap/coveredGap:NaN)<<'\t'<<(gap?1-newGap/gap:NaN)<<'\t'<<nh[0].st.mean()<<'\t'<<nh[1].st.mean()<<'\t'<<d[0].core<<'\t'<<d[1].core<<'\t'<<*std::max_element(rw.weights[0].begin(),rw.weights[0].end())<<'\t'<<*std::max_element(rw.weights[1].begin(),rw.weights[1].end())<<'\t'<<(rw.neff[0]<100||rw.neff[1]<100||rw.neff[0]<.1*rw.covered[0]||rw.neff[1]<.1*rw.covered[1]?"unstable_importance_weights":double(rw.covered[0])/ev[0].size()<.5||double(rw.covered[1])/ev[1].size()<.5?"limited_common_support":"adequate_effective_sample")<<'\n';
    }
   }
   TCanvas canvas(("event_components_"+prefix).c_str(),"",2250,1500);canvas.SetBatch(true);canvas.SetCanvasSize(2250,1500);canvas.Divide(3,2);
   const std::array<std::string,6> distributionFields={"hits","central_hits","peripheral_hits","energy_per_hit_MeV","central_fraction","energy"};
   const std::array<std::string,6> distributionLabels={"N_{hit}","Central N_{hit}","Peripheral N_{hit}","E / N_{hit} [MeV]","Central energy / total energy","Energy / sample mean"};
   for(int panel=0;panel<6;++panel){auto v=std::array<std::vector<double>,2>{values(ev[0],distributionFields[panel]),values(ev[1],distributionFields[panel])};if(panel==5)for(int s=0;s<2;++s)for(double&x:v[s])x/=baseline[s].st.mean();overlay(canvas,panel+1,v,std::to_string(beam)+" GeV | "+distributionFields[panel],distributionLabels[panel],root,prefix+distributionFields[panel],nullptr,true,panel<3);}
   fastPng(canvas,(dest/"event_distributions.png").c_str());canvas.Clear();canvas.Divide(3,2);
   const std::array<std::string,6> axes={"hits","x_mm","y_mm","back5_fraction","depth_layer",std::isfinite(ev[0].front().temperature)?"temperature_C":"temperature_correction_log"};
   const std::array<std::string,6> labels={"N_{hit}","X shower CoG [mm]","Y shower CoG [mm]","Last 5-layer energy fraction","Mean shower depth [layer]",std::isfinite(ev[0].front().temperature)?"Energy-weighted event temperature [C]":"log(E / E without temperature correction)"};
   for(int panel=0;panel<6;++panel)profile(canvas,panel+1,ev,original,axes[panel],labels[panel],root,prefix+"profile_"+axes[panel]);
   fastPng(canvas,(dest/"event_response_profiles.png").c_str());canvas.Clear();canvas.Divide(3,2);
   const std::array<std::string,6> plotModels={"baseline","nhit","position","position_nhit","position_leak","position_nhit_leak"};
   for(int panel=0;panel<6;++panel){auto v=correctedValues.at(plotModels[panel]);for(int s=0;s<2;++s){Distribution d(v[s]);for(double&x:v[s])x/=d.st.mean();}overlay(canvas,panel+1,v,std::to_string(beam)+" GeV | event-held-out "+plotModels[panel],"Energy / sample mean",root,prefix+"cv_"+plotModels[panel]);canvas.cd(panel+1);Distribution d(correctedValues.at(plotModels[panel])[0]),m(correctedValues.at(plotModels[panel])[1]);text(.89,"RMS/mean data="+number(100*d.st.resolution(),2)+"%, MC="+number(100*m.st.resolution(),2)+"%");}
   fastPng(canvas,(dest/"diagnostic_resolution_components.png").c_str());canvas.Clear();
   canvas.Divide(3,2);
   const std::array<std::string,3> weightModes={"mc_to_data_nhit","mc_to_data_position","mc_to_data_position_nhit"};
   for(int panel=0;panel<3;++panel){
    const auto& rw=weightedExamples.at(weightModes[panel]);
    std::array<std::vector<double>,2> n={values(ev[0],"hits"),values(ev[1],"hits")};
    overlay(canvas,panel+1,n,std::to_string(beam)+" GeV | "+weightModes[panel],"N_{hit}",root,prefix+"weighted_N_"+weightModes[panel],&rw.weights,true,true);
    canvas.cd(panel+1);text(.89,"Effective events data="+number(rw.neff[0],0)+", MC="+number(rw.neff[1],0));
    auto e=original;for(int s=0;s<2;++s){Distribution d(e[s],rw.weights[s]);for(double&x:e[s])x/=d.st.mean();}
    overlay(canvas,panel+4,e,"Same weighted cohort","Energy / weighted sample mean",root,prefix+"weighted_E_"+weightModes[panel],&rw.weights,true);
    canvas.cd(panel+4);text(.89,"Covered data="+number(100.*rw.covered[0]/ev[0].size(),1)+"%, MC="+number(100.*rw.covered[1]/ev[1].size(),1)+"%");
   }
   fastPng(canvas,(dest/"event_reweighting.png").c_str());canvas.Clear();

   std::cout<<beam<<" GeV "<<selection<<" events "<<ev[0].size()<<"/"<<ev[1].size()<<", resolution "<<baseline[0].st.resolution()<<"/"<<baseline[1].st.resolution()<<std::endl;
  }
  need(before==std::make_pair(fs::file_size(eventPath),fs::last_write_time(eventPath)),"Input events changed while reading");
 }
 root.Close();return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
