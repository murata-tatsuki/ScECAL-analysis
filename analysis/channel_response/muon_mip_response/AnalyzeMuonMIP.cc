// Muon response diagnostics. Input files and calibration factors are read-only.
#include "MuonMIPCommon.hh"
#include "CurrentCalibration.hh"
#include "ExactFitCache.hh"
#include "MuonMIPPlotStyle.hh"
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <cstdlib>
#include <TBranch.h>
#include <unistd.h>
#include "../adc_energy/FastPng.hh"
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TParameter.h>
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <TMath.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace fs=std::filesystem;
namespace {
void need(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}
using Parameters=std::map<int,mmr::Parameter>;
using RootFile=std::unique_ptr<TFile>;
RootFile open(const fs::path& path,const char* mode="READ"){
 RootFile f(TFile::Open(path.c_str(),mode));
 need(f&&!f->IsZombie(),"Cannot open "+path.string());return f;
}
template<class T>void bind(TTree& tree,const char* name,T* value){
 need(tree.GetBranch(name),std::string("Missing ")+tree.GetName()+"/"+name);
 tree.SetBranchStatus(name,1);need(tree.SetBranchAddress(name,value)>=0,std::string("Branch type differs: ")+name);
}
void cacheEnabledBranches(TTree& tree){
 tree.SetCacheSize(64LL*1024*1024);
 tree.DropBranchFromCache("*",true);
 auto* branches=tree.GetListOfBranches();
 for(int i=0;i<branches->GetEntries();++i){
  auto* branch=dynamic_cast<TBranch*>(branches->At(i));
  if(branch&&tree.GetBranchStatus(branch->GetName()))tree.AddBranchToCache(branch->GetName(),true);
 }
 tree.StopCacheLearningPhase();
}
void read(TTree& tree,Long64_t entry){need(tree.GetEntry(entry)>0,"Cannot read "+std::string(tree.GetName())+" entry "+std::to_string(entry));}
std::ofstream table(const fs::path& path){fs::create_directories(path.parent_path());std::ofstream f(path);need(bool(f),"Cannot write "+path.string());f<<std::setprecision(12);return f;}
bool connected(int id){return id>=0&&id/100000<30&&id%100000/10000<6&&id%100<36&&(id%100000/10000!=5||id%100<30);}
int sipm(int id){return (id/100000<4||id/100000>27)?15:10;}
double correction(double temperature,int id){return 1-(temperature-20)*(sipm(id)==10?1.6/135:3.5/230);}
std::string sampleName(int sample){return sample?"mc":"data";}
struct Crossing {
 int sample=0,run=0,event=0,event_time=0,cellid=0,observed=0,ambiguous=0,calib_ok=0,response_ok=0;
 Long64_t event_index=0;
 double raw_hg=NAN,raw_lg=NAN,saved_hg=NAN,saved_energy=NAN,temperature=NAN,sx=NAN,sy=NAN,length=NAN,longitudinal=NAN,transverse=NAN,track_rms=NAN,truth_edep=NAN;
 void attach(TTree& t,bool response=false){
  t.SetBranchStatus("*",0);
  for(auto x:std::vector<std::pair<const char*,int*>>{{"sample",&sample},{"run",&run},{"cellid",&cellid},{"observed",&observed},{"ambiguous",&ambiguous},{"calib_ok",&calib_ok},{"response_ok",&response_ok}})bind(t,x.first,x.second);
  for(auto x:std::vector<std::pair<const char*,double*>>{{"raw_hg",&raw_hg},{"temperature",&temperature},{"sx",&sx},{"sy",&sy},{"length",&length},{"longitudinal",&longitudinal},{"transverse",&transverse}})bind(t,x.first,x.second);
  if(response)for(auto x:std::vector<std::pair<const char*,double*>>{{"raw_lg",&raw_lg},{"saved_hg",&saved_hg},{"saved_energy",&saved_energy},{"truth_edep",&truth_edep}})bind(t,x.first,x.second);
 }
};
using GeometryKey=std::tuple<int,int,int,int,int,int>;
struct Counts {std::array<long long,2> n{};};
GeometryKey geometryKey(const Crossing& h,double longBin,double transverseBin,double slopeBin){
 return {h.run,h.cellid,int(std::floor(h.longitudinal/longBin)),int(std::floor(h.transverse/transverseBin)),int(std::floor(h.sx/slopeBin)),int(std::floor(h.sy/slopeBin))};
}
enum Stage {Pedestal,Path,Temperature,SavedMIP,Truth,Linear,Stages};
const std::array<std::string,Stages> stageName={"pedestal_adc","path_adc","temperature_adc","saved_hg_mip","detected_truth_mev","linear_prediction_adc"};
const std::array<std::string,Stages> units={"ADC","ADC","ADC","MIP","MeV","ADC"};
const std::array<std::string,Stages> stageTitle={"Pedestal subtracted HG","HG, 2 mm path normalization","HG, path and 20 C normalization","Saved HG energy, path normalized","Muon truth, detected crossings","Linear ADC prediction, path normalized"};
double stageScale(int stage,double adcRef){return stage==SavedMIP?1:stage==Truth?.305:adcRef;}
struct Spectrum {
 std::array<std::unique_ptr<TH1D>,Stages> h;
 std::array<double,Stages> maxLower,minUpper;
 std::array<std::unique_ptr<TProfile>,4> profile;
 double sumw=0,sumw2=0;long long hits=0;
 Spectrum(){maxLower.fill(-INFINITY);minUpper.fill(INFINITY);}
};
struct Efficiency {long long crossings=0,observed=0,ambiguous=0,calib=0,response=0;};
struct Closure {long long checkedHG=0,checkedSelected=0,badHG=0,badSelected=0;double maxHG=0,maxSelected=0;};
struct DependencySpectrum {
 std::unique_ptr<TH1D> h;
 double lowerFactor=-INFINITY,upperFactor=INFINITY,sumx=0,sumw=0;
};
using DependencyKey=std::tuple<int,int,int,int>; // channel,matched,axis,bin
using DependencyGraphKey=std::tuple<int,int,int>; // channel,matched,axis
const std::array<std::string,4> dependencyAxis={"longitudinal","transverse","temperature","run"};
const std::array<std::string,4> dependencyAxisLabel={"Strip longitudinal position [mm]","Strip transverse position [mm]","Temperature [C]","Run"};
struct Channel {
 std::array<std::array<Spectrum,2>,2> spectra; // sample, geometry unmatched/matched
 std::array<Efficiency,2> efficiency;
 std::array<Closure,2> closure;
 std::array<std::unique_ptr<TProfile>,2> efficiencyLong;
 std::unique_ptr<TH1D> generatedTruth,generatedTruthDetected;
 double adcRef=0;
};
double refADC(int id,const std::array<Parameters,2>& pars){
 double a=pars[0].at(id).mip,b=pars[1].at(id).mip;
 return (std::isfinite(a)&&a>0&&std::isfinite(b)&&b>0)?std::sqrt(a*b):(std::isfinite(a)&&a>0?a:b);
}
void initialize(Spectrum& sp,int id,int sample,int matched,double adcRef){
 if(sp.h[0])return;
 for(int stage=0;stage<Stages;++stage){
  std::string name=sampleName(sample)+"_"+(matched?"matched_":"unmatched_")+std::to_string(id)+"_"+stageName[stage];
  double scale=stageScale(stage,adcRef);
  sp.h[stage]=std::make_unique<TH1D>(name.c_str(),(stageTitle[stage]+";"+units[stage]+";Weighted hits / bin").c_str(),200,0,8*scale);
  sp.h[stage]->SetDirectory(nullptr);sp.h[stage]->Sumw2();
 }
 for(int p=0;p<4;++p){
  std::string name=sampleName(sample)+"_"+std::to_string(id)+"_"+std::to_string(matched)+"_profile"+std::to_string(p);
  const double low=p==0?-22.5:p==1?-2.5:p==2?10:0;
  const double high=p==0?22.5:p==1?2.5:p==2?40:400;
  sp.profile[p]=std::make_unique<TProfile>(name.c_str(),"HG response / current MIP",p==3?400:30,low,high);
  sp.profile[p]->SetDirectory(nullptr);sp.profile[p]->Sumw2();
 }
}
struct Fit {
 std::string status="empty";
 double entries=0,sumw=0,sumw2=0,neff=0,fitNeff=0,underflow=0,overflow=0,low=NAN,high=NAN;
 double mpv=NAN,mpvError=NAN,landau=NAN,gaussian=NAN,area=NAN,chi2ndf=NAN,normalized=NAN;
 int rootStatus=-999,covariance=-1,ndf=0;bool valid=false;
};
double langau(double* x,double* p){
 const double invsqrt=.3989422804014327,mpShift=-.22278298;
 const double center=p[1]-mpShift*p[0],low=x[0]-5*p[3],high=x[0]+5*p[3],step=(high-low)/60;
 double sum=0;for(int i=0;i<60;++i){double xx=low+(i+.5)*step;sum+=TMath::Landau(xx,center,p[0],false)/p[0]*TMath::Gaus(x[0],xx,p[3],false);}
 return p[2]*step*sum*invsqrt/p[3];
}
unsigned long long rootFitCalls=0;
Fit fitHistogramUncached(const TH1D* h,double low,double high,double scale,double normalization,double minHits,double maxChi2){
 Fit result;if(!h)return result;
 result.entries=h->GetEntries();result.low=low;result.high=high;
 for(int b=0;b<=h->GetNbinsX()+1;++b){result.sumw+=h->GetBinContent(b);result.sumw2+=std::pow(h->GetBinError(b),2);}
 result.neff=result.sumw2>0?result.sumw*result.sumw/result.sumw2:0;
 result.underflow=result.sumw>0?h->GetBinContent(0)/result.sumw:0;
 result.overflow=result.sumw>0?h->GetBinContent(h->GetNbinsX()+1)/result.sumw:0;
 if(result.neff<minHits){result.status="insufficient_total_effective_hits";return result;}
 if(!std::isfinite(low)||!std::isfinite(high)){result.status="missing_threshold";return result;}
 const double bin=h->GetBinWidth(1);
 low=std::max(bin*.5,low);high=std::min(high,h->GetXaxis()->GetXmax()-bin*.5);
 result.low=low;result.high=high;
 if(high-low<12*bin){result.status="plateau_range_too_short";return result;}
 int first=h->FindFixBin(low),last=h->FindFixBin(high),peak=first;double w=0,w2=0;
 for(int b=first;b<=last;++b){w+=h->GetBinContent(b);w2+=std::pow(h->GetBinError(b),2);if(h->GetBinContent(b)>h->GetBinContent(peak))peak=b;}
 result.fitNeff=w2>0?w*w/w2:0;if(result.fitNeff<minHits){result.status="insufficient_plateau_effective_hits";return result;}
 double initial=std::clamp(h->GetBinCenter(peak),low+2*bin,high-2*bin),area=std::max(w*bin,1e-3);
 std::unique_ptr<TH1D> copy(static_cast<TH1D*>(h->Clone("fit_work")));copy->SetDirectory(nullptr);
 TF1 f("muon_langau",langau,low,high,4);f.SetParNames("Landau width","Landau MPV","Area","Gaussian sigma");
 f.SetParameters(.12*scale,initial,area,.08*scale);
 f.SetParLimits(0,.005*scale,scale);f.SetParLimits(1,low,high);f.SetParLimits(2,.001*area,100*area);f.SetParLimits(3,.003*scale,1.5*scale);
 ++rootFitCalls;
 auto fit=copy->Fit(&f,"QRS0N","",low,high);
 result.rootStatus=int(fit);if(fit.Get())result.covariance=fit->CovMatrixStatus();
 result.mpv=f.GetParameter(1);result.mpvError=f.GetParError(1);result.landau=f.GetParameter(0);result.gaussian=f.GetParameter(3);result.area=f.GetParameter(2);result.ndf=f.GetNDF();result.chi2ndf=result.ndf>0?f.GetChisquare()/result.ndf:NAN;result.normalized=result.mpv/normalization;
 if(result.rootStatus){result.status="fit_failed";return result;}
 if(result.covariance<2){result.status="poor_covariance";return result;}
 if(result.ndf<=0||!std::isfinite(result.chi2ndf)||result.chi2ndf>maxChi2){result.status="poor_chi2";return result;}
 if(!std::isfinite(result.mpv)||result.mpv<=0||!std::isfinite(result.mpvError)||result.mpvError<=0){result.status="invalid_mpv";return result;}
 if(result.mpv<=low+2*bin||result.mpv>=high-2*bin){result.status="mpv_near_plateau_boundary";return result;}
 for(int p=0;p<4;++p){double lo=0,hi=0;f.GetParLimits(p,lo,hi);double value=f.GetParameter(p);if(hi>lo&&(value-lo<1e-3*(hi-lo)||hi-value<1e-3*(hi-lo))){result.status="parameter_boundary";return result;}}
 result.status="ok";result.valid=true;return result;
}
std::unique_ptr<mmr::ExactFitCache<Fit>> fitCache;
Fit fitHistogram(const TH1D* h,double low,double high,double scale,double normalization,double minHits,double maxChi2){
 need(bool(fitCache),"Fit cache was not initialized");
 return fitCache->get(h,low,high,scale,normalization,minHits,maxChi2,[&]{
  return fitHistogramUncached(h,low,high,scale,normalization,minHits,maxChi2);
 });
}
using FitKey=std::tuple<int,int,int,int,int>; // channel,sample,matched,stage,plateau slot
std::unique_ptr<TF1> fitCurve(const Fit& result,const std::string& name,double normalization=1){
 if(!result.valid)return {};
 auto f=std::make_unique<TF1>(name.c_str(),langau,result.low,result.high,4);
 f->SetParNames("Landau width","Landau MPV","Area","Gaussian sigma");
 f->SetParameters(result.landau,result.mpv,result.area/normalization,result.gaussian);
 f->SetNpx(500);
 // Save a sampled function so ROOT can redraw it without this executable.
 f->Save(result.low,result.high,0,0,0,0);
 return f;
}
std::pair<double,double> bounds(const Spectrum& spectrum,const mmr::Parameter& p,int stage,double sigmas,double hgLimit){
 double low=spectrum.maxLower[stage],high=spectrum.minUpper[stage];
 if(stage==Truth||stage==Linear)return {0,high};
 if(!std::isfinite(p.threshold)||!std::isfinite(p.sigma)||p.sigma<0)return {NAN,NAN};
 // maxLower stores the largest multiplicative path/temperature conversion factor.
 low=(p.threshold+sigmas*p.sigma)*low;
 high=(hgLimit-p.ph)*high;
 return {std::max(0.,low),high};
}
void writeFit(std::ostream& out,int id,int sample,const std::string& condition,int stage,double sigma,const Fit& f){
 out<<id<<'\t'<<id/100000<<'\t'<<sipm(id)<<'\t'<<sampleName(sample)<<'\t'<<condition<<'\t'<<stageName[stage]<<'\t'<<units[stage]<<'\t'<<sigma<<'\t'<<f.entries<<'\t'<<f.sumw<<'\t'<<f.neff<<'\t'<<f.fitNeff<<'\t'<<f.underflow<<'\t'<<f.overflow<<'\t'<<f.low<<'\t'<<f.high<<'\t'<<f.mpv<<'\t'<<f.mpvError<<'\t'<<f.normalized<<'\t'<<f.landau<<'\t'<<f.gaussian<<'\t'<<f.area<<'\t'<<f.chi2ndf<<'\t'<<f.rootStatus<<'\t'<<f.covariance<<'\t'<<f.status<<'\n';
}
double median(std::vector<double> values){if(values.empty())return NAN;std::sort(values.begin(),values.end());size_t n=values.size();return n%2?values[n/2]:.5*(values[n/2-1]+values[n/2]);}
std::set<int> requestedChannels(const std::string& value){
 std::set<int> result;if(value=="all")return result;std::istringstream input(value);
 for(std::string field;std::getline(input,field,',');){int id=std::stoi(field);need(connected(id),"Invalid requested physical CellID "+field);result.insert(id);}need(!result.empty(),"Empty channels selection");return result;
}
void annotation(double y,const std::string& value,int color=kBlack){TLatex text;text.SetNDC();text.SetTextSize(.030);text.SetTextColor(color);text.DrawLatex(.13,y,value.c_str());}
void overlay(TCanvas& canvas,int pad,const TH1D* d,const TH1D* m,const std::string& title,const std::string& unit,std::vector<std::unique_ptr<TH1D>>& drawn,const std::string& firstLabel="data",const std::string& secondLabel="MC"){
 mmr_plot::drawResponsePanel(canvas,pad,d,m,title,unit,drawn,nullptr,nullptr,nullptr,nullptr,firstLabel,secondLabel);
}
mmr_plot::FitInfo plotFitInfo(const Fit& fit){
 return {fit.status,fit.valid,fit.mpv,fit.mpvError,fit.landau,fit.gaussian,fit.chi2ndf};
}
void profiles(TCanvas& canvas,const Channel& ch,int matched,const fs::path& path){
 canvas.Clear();canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(2,2);const std::array<std::string,4> axes={"Strip longitudinal position [mm]","Strip transverse position [mm]","Temperature [C]","Run"};
 for(int j=0;j<4;++j){canvas.cd(j+1);TProfile* first=nullptr;for(int s=0;s<2;++s){auto* profile=ch.spectra[s][matched].profile[j].get();if(!profile||!profile->GetEntries())continue;profile->SetTitle(("Channel HG response;"+axes[j]+";Path/temperature corrected ADC / current MIP").c_str());profile->SetStats(false);profile->SetLineColor(s?kRed+1:kBlack);profile->SetMarkerColor(s?kRed+1:kBlack);profile->SetMarkerStyle(s?24:20);profile->Draw(first?"SAME":"");if(!first)first=profile;}TLegend legend(.73,.73,.91,.87);for(int s=0;s<2;++s)if(ch.spectra[s][matched].profile[j]&&ch.spectra[s][matched].profile[j]->GetEntries())legend.AddEntry(ch.spectra[s][matched].profile[j].get(),s?"MC":"data","lp");legend.DrawClone();}
 fastPng(canvas,path.c_str());canvas.Clear();
}
void mpvProfiles(TCanvas& canvas,int cid,int matched,
 const std::map<DependencyGraphKey,std::array<std::unique_ptr<TGraphAsymmErrors>,2>>& graphs,
 const fs::path& path){
 std::array<std::array<const TGraphAsymmErrors*,2>,4> panels{};
 for(int axis=0;axis<4;++axis){auto found=graphs.find({cid,matched,axis});if(found!=graphs.end())for(int sample=0;sample<2;++sample)panels[axis][sample]=found->second[sample].get();}
 mmr_plot::drawMpvProfiles(canvas,cid,bool(matched),panels,path.string());
}
} // namespace

int main(int argc,char** argv)try{
 need(argc==4,"Usage: AnalyzeMuonMIP flat.root config.tsv output_directory");
 gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetCanvasPreferGL(false);
 mmr::Config config(argv[2]);fs::path out=argv[3];fs::create_directories(out);
 const std::string energyDir=std::to_string(config.integer("energy",100))+"GeV";
 const bool allow=config.integer("allow-recovered",0);
 std::array<Parameters,2> pars;
 for(int s=0;s<2;++s){std::string prefix=sampleName(s)+"-";pars[s]=mmr::loadParameters(config.text(prefix+"pedestal"),config.text(prefix+"mip"),config.text(prefix+"hl"),config.text(prefix+"threshold"),allow);}
 auto digi=mmr::loadMip(config.text("digi-mip",config.text("mc-mip")),allow);
 const double longBin=config.number("long-bin",5),transverseBin=config.number("transverse-bin",1.5),slopeBin=config.number("slope-bin",.01);
 const long long minBin=config.integer("min-bin",3);const double minFit=config.number("min-fit-hits",100),plateau=config.number("plateau-sigmas",3),maxChi=config.number("fit-chi2max",3),hgLimit=config.number("hg-limit",2600);
 need(longBin>0&&transverseBin>0&&slopeBin>0&&minBin>0&&minFit>=20&&plateau>=0&&maxChi>0&&hgLimit>0,"Invalid bin/fit configuration");
 const int workerCount=static_cast<int>(config.integer("workers",4));
 need(workerCount>=1&&workerCount<=64,"workers must be between 1 and 64");
 const long long cacheMiB=config.integer("fit-cache-mb",32);
 need(cacheMiB>=0&&cacheMiB<=4096,"Invalid fit cache memory budget");
 fitCache=std::make_unique<mmr::ExactFitCache<Fit>>(config.integer("fit-cache",1)?static_cast<size_t>(cacheMiB)*1024*1024:0);
 const auto requested=requestedChannels(config.text("channels","all"));
 const bool dependencyFits=config.integer("dependency-fits",1);
 const std::array<double,4> dependencyWidths={longBin,transverseBin,config.number("temperature-bin",1),1};
 need(dependencyWidths[2]>0,"temperature-bin must be positive");
 std::map<DependencyKey,std::array<DependencySpectrum,2>> dependency;
 auto input=open(argv[1]);need(allow||!input->TestBit(TFile::kRecovered),"Recovered flat ROOT");
 auto* crossings=dynamic_cast<TTree*>(input->Get("Crossings"));need(crossings,"Missing Crossings");Crossing row;row.attach(*crossings);cacheEnabledBranches(*crossings);
 std::map<GeometryKey,Counts> counts;std::map<int,Channel> channels;
 auto eligible=[&](){return row.response_ok&&row.observed&&!row.ambiguous&&row.calib_ok&&connected(row.cellid)&&row.sample>=0&&row.sample<2&&std::isfinite(row.raw_hg)&&row.raw_hg<hgLimit&&std::isfinite(row.length)&&row.length>0&&std::isfinite(row.sx)&&std::isfinite(row.sy)&&std::isfinite(row.longitudinal)&&std::isfinite(row.transverse)&&std::isfinite(row.temperature)&&correction(row.temperature,row.cellid)>0;};
 for(Long64_t i=0;i<crossings->GetEntries();++i){read(*crossings,i);if(eligible())++counts[geometryKey(row,longBin,transverseBin,slopeBin)].n[row.sample];}
 auto bins=table(out/"geometry_bins.tsv");bins<<"run\tCellID\tlong_bin\ttransverse_bin\tsx_bin\tsy_bin\tdata_hits\tmc_hits\tcommon_target\tdata_weight\tmc_weight\n";
 long long commonBins=0;double commonTarget=0;
 for(auto& [key,c]:counts){long long target=c.n[0]>=minBin&&c.n[1]>=minBin?std::min(c.n[0],c.n[1]):0;if(target){++commonBins;commonTarget+=target;}bins<<std::get<0>(key)<<'\t'<<std::get<1>(key)<<'\t'<<std::get<2>(key)<<'\t'<<std::get<3>(key)<<'\t'<<std::get<4>(key)<<'\t'<<std::get<5>(key)<<'\t'<<c.n[0]<<'\t'<<c.n[1]<<'\t'<<target<<'\t'<<(c.n[0]?double(target)/c.n[0]:0)<<'\t'<<(c.n[1]?double(target)/c.n[1]:0)<<'\n';}
 auto parameters=table(out/"inputs/calibration_parameters.tsv");parameters<<"CellID\tsample\tSiPM_um\tadopted_MIP_ADC\traw_MIP_ADC\tchi2_ndf\tMIP_fallback\tpedestal_HG\tpedestal_LG\tHL_gain\tHL_offset_HG_ADC\tthreshold_subtracted_ADC\tthreshold_sigma_ADC\tthreshold_fallback\tdigi_MIP_ADC\tpedestal_present\n";
 for(const auto& [id,p]:pars[0])for(int s=0;s<2;++s){auto& q=pars[s].at(id);parameters<<id<<'\t'<<sampleName(s)<<'\t'<<sipm(id)<<'\t'<<q.mip<<'\t'<<q.rawMip<<'\t'<<q.chi2ndf<<'\t'<<q.mipFallback<<'\t'<<q.ph<<'\t'<<q.pl<<'\t'<<q.gain<<'\t'<<q.offset<<'\t'<<q.threshold<<'\t'<<q.sigma<<'\t'<<q.thresholdFallback<<'\t'<<digi.at(id).mip<<'\t'<<q.pedPresent<<'\n';}
 std::array<long long,2> eligibleHits{},matchedHits{};std::array<double,2> matchedWeights{};
 row.attach(*crossings,true);cacheEnabledBranches(*crossings);
 for(Long64_t i=0;i<crossings->GetEntries();++i){
  read(*crossings,i);need(row.sample==0||row.sample==1,"Invalid sample");if(!connected(row.cellid))continue;
  auto& ch=channels[row.cellid];if(!ch.adcRef)ch.adcRef=refADC(row.cellid,pars);
  auto& eff=ch.efficiency[row.sample];++eff.crossings;eff.observed+=row.observed&&!row.ambiguous;eff.ambiguous+=row.ambiguous;eff.calib+=row.calib_ok;eff.response+=eligible();
  if(!ch.efficiencyLong[row.sample]){std::string n=sampleName(row.sample)+"_"+std::to_string(row.cellid)+"_efficiency_long";ch.efficiencyLong[row.sample]=std::make_unique<TProfile>(n.c_str(),"Conditional observed crossing fraction;Strip longitudinal position [mm];Observed / predicted crossing",30,-22.5,22.5);ch.efficiencyLong[row.sample]->SetDirectory(nullptr);}
  if(std::isfinite(row.longitudinal))ch.efficiencyLong[row.sample]->Fill(row.longitudinal,row.observed&&!row.ambiguous?1:0);
  auto& p=pars[row.sample].at(row.cellid);double C=correction(row.temperature,row.cellid),q=row.raw_hg-p.ph;
  auto& closure=ch.closure[row.sample];
  if(row.observed&&!row.ambiguous&&row.calib_ok&&std::isfinite(q)&&std::isfinite(C)&&C>0){
   double expected=.305/p.mip*C*q;if(std::isfinite(row.saved_hg)){double error=std::abs(expected-row.saved_hg);++closure.checkedHG;closure.maxHG=std::max(closure.maxHG,error);closure.badHG+=error>1e-6;}
   // Saved calibration uses the production HG/LG transition independently
   // of the diagnostic HG selection limit.
   double reconstructed=row.raw_hg<2600?expected:.305/p.mip*(C*(row.raw_lg-p.pl)*p.gain+p.offset);
   if(std::isfinite(reconstructed)&&std::isfinite(row.saved_energy)){double error=std::abs(reconstructed-row.saved_energy);++closure.checkedSelected;closure.maxSelected=std::max(closure.maxSelected,error);closure.badSelected+=error>1e-6;}
  }
  if(!eligible())continue;
  ++eligibleHits[row.sample];double path=2/row.length;
  const auto& count=counts.at(geometryKey(row,longBin,transverseBin,slopeBin));long long target=count.n[0]>=minBin&&count.n[1]>=minBin?std::min(count.n[0],count.n[1]):0;
  double weight=target?double(target)/count.n[row.sample]:0;if(weight){++matchedHits[row.sample];matchedWeights[row.sample]+=weight;}
  std::array<double,Stages> value={q,q*path,q*path*C,row.saved_hg/.305*path,row.truth_edep*path,row.truth_edep/.305*digi.at(row.cellid).mip*path};
  std::array<double,Stages> factor={1,path,path*C,path*C/p.mip,1,1};
  for(int match=0;match<2;++match){double w=match?weight:1;if(w<=0)continue;auto& sp=ch.spectra[row.sample][match];initialize(sp,row.cellid,row.sample,match,ch.adcRef);++sp.hits;sp.sumw+=w;sp.sumw2+=w*w;
   for(int st=0;st<Stages;++st){if((st==Truth||st==Linear)&&!row.sample)continue;if(!std::isfinite(value[st]))continue;sp.h[st]->Fill(value[st],w);
    if(st!=Truth&&st!=Linear){sp.maxLower[st]=std::max(sp.maxLower[st],factor[st]);sp.minUpper[st]=std::min(sp.minUpper[st],factor[st]);}
    else {sp.maxLower[st]=0;sp.minUpper[st]=sp.h[st]->GetXaxis()->GetXmax();}
   }
   std::array<double,4> pos={row.longitudinal,row.transverse,row.temperature,double(row.run)};for(int k=0;k<4;++k)sp.profile[k]->Fill(pos[k],value[Temperature]/p.mip,w);
   if(dependencyFits){
    // Inclusive-fit eligibility is applied later; prefill preserves each bin's
    // exact original event order and arithmetic, without a third ROOT scan.
    const double depFactor=path*C,depValue=q*depFactor;
    for(int axis=0;axis<4;++axis){
     int bin=int(std::floor(pos[axis]/dependencyWidths[axis]));
     auto& dep=dependency[{row.cellid,match,axis,bin}][row.sample];
     if(!dep.h){
      std::string name=sampleName(row.sample)+"_"+(match?"matched_":"unmatched_")+"cell"+std::to_string(row.cellid)+"_"+dependencyAxis[axis]+"_bin"+std::to_string(bin);
      dep.h=std::make_unique<TH1D>(name.c_str(),("20 C HG ADC within "+dependencyAxis[axis]+" bin;Path and temperature corrected ADC;Weighted hits / bin").c_str(),200,0,8*ch.adcRef);
      dep.h->SetDirectory(nullptr);dep.h->Sumw2();
     }
     dep.h->Fill(depValue,w);dep.sumw+=w;dep.sumx+=w*pos[axis];
     dep.lowerFactor=std::max(dep.lowerFactor,depFactor);dep.upperFactor=std::min(dep.upperFactor,depFactor);
    }
   }
  }
 }
 crossings->ResetBranchAddresses();
 auto* truths=dynamic_cast<TTree*>(input->Get("TruthHits"));need(truths,"Missing TruthHits: independent generated-muon reference required");
 int truthSample=0,id=0,detected=0;double energy=NAN,length=NAN;truths->SetBranchStatus("*",0);bind(*truths,"sample",&truthSample);bind(*truths,"cellid",&id);bind(*truths,"truth_edep",&energy);bind(*truths,"truth_length",&length);bind(*truths,"truth_detected",&detected);
 cacheEnabledBranches(*truths);
 long long generatedValid=0,generatedInvalidChord=0;
 for(Long64_t i=0;i<truths->GetEntries();++i){read(*truths,i);if(!connected(id)||truthSample!=1)continue;if(!std::isfinite(length)||length<=0||!std::isfinite(energy)||energy<=0){++generatedInvalidChord;continue;}
  auto& ch=channels[id];if(!ch.adcRef)ch.adcRef=refADC(id,pars);if(!ch.generatedTruth){ch.generatedTruth=std::make_unique<TH1D>(("generated_truth_"+std::to_string(id)).c_str(),"Generated muon truth, geometrically crossed strip;2 mm normalized edep [MeV];Hits / bin",200,0,8*.305);ch.generatedTruth->SetDirectory(nullptr);ch.generatedTruth->Sumw2();ch.generatedTruthDetected=std::make_unique<TH1D>(("generated_detected_"+std::to_string(id)).c_str(),"Generated truth with raw detection;2 mm normalized edep [MeV];Hits / bin",200,0,8*.305);ch.generatedTruthDetected->SetDirectory(nullptr);ch.generatedTruthDetected->Sumw2();}
  ch.generatedTruth->Fill(energy*2/length);if(detected)ch.generatedTruthDetected->Fill(energy*2/length);++generatedValid;
 }truths->ResetBranchAddresses();
 auto coverage=table(out/"geometry_coverage.tsv");coverage<<"sample\teligible_HG_hits\tmatched_hits\tmatched_hit_fraction\tmatched_sum_weights\tcommon_bins\tcommon_target\tgenerated_truth_valid\tgenerated_truth_invalid_primary_chord\n";
 for(int s=0;s<2;++s)coverage<<sampleName(s)<<'\t'<<eligibleHits[s]<<'\t'<<matchedHits[s]<<'\t'<<(eligibleHits[s]?double(matchedHits[s])/eligibleHits[s]:NAN)<<'\t'<<matchedWeights[s]<<'\t'<<commonBins<<'\t'<<commonTarget<<'\t'<<(s?generatedValid:0)<<'\t'<<(s?generatedInvalidChord:0)<<'\n';
 need(std::abs(matchedWeights[0]-matchedWeights[1])<1e-6*std::max(1.,commonTarget),"Geometry matching weights do not close");
 auto eventTable=table(out/"event_summary.tsv");eventTable<<"sample\tcollected_events\tcalib_paired_finite_events\tunavailable_event_energy\tmean_saved_energy_MeV\tmean_positive_saved_energy_MeV\tmean_saved_hit_count\n";
 if(auto* events=dynamic_cast<TTree*>(input->Get("Events"))){
  int sample=0,hits=0,paired=0;double energy=0,positive=0;
  events->SetBranchStatus("*",0);bind(*events,"sample",&sample);bind(*events,"calib_paired",&paired);bind(*events,"event_energy",&energy);bind(*events,"event_positive_energy",&positive);bind(*events,"event_hits",&hits);
  cacheEnabledBranches(*events);
  std::array<long long,2> total{},n{};std::array<double,2> es{},ps{},hs{};
  for(Long64_t i=0;i<events->GetEntries();++i){read(*events,i);need(sample>=0&&sample<2,"Invalid Events sample");++total[sample];if(!paired||!std::isfinite(energy)||!std::isfinite(positive)||hits<0)continue;++n[sample];es[sample]+=energy;ps[sample]+=positive;hs[sample]+=hits;}
  for(int s=0;s<2;++s)eventTable<<sampleName(s)<<'\t'<<total[s]<<'\t'<<n[s]<<'\t'<<total[s]-n[s]<<'\t'<<(n[s]?es[s]/n[s]:NAN)<<'\t'<<(n[s]?ps[s]/n[s]:NAN)<<'\t'<<(n[s]?hs[s]/n[s]:NAN)<<'\n';
  events->ResetBranchAddresses();
 }
 const Long64_t crossingEntries=crossings->GetEntries();
 bins.close();parameters.close();coverage.close();eventTable.close();
 input->Close();input.reset();
 // Matching-bin counts are no longer needed by any fit worker.
 counts.clear();
 struct rusage preparedUsage{};
 getrusage(RUSAGE_SELF,&preparedUsage);
 std::cout<<"Prepared "<<channels.size()<<" channels and "<<dependency.size()<<" dependency bins; peak preparation RSS="<<preparedUsage.ru_maxrss<<" KiB; starting "<<workerCount<<" isolated fit workers"<<std::endl;
 std::map<int,int> channelOwner;int nextOwner=0;
 for(const auto& item:channels){channelOwner[item.first]=nextOwner++%workerCount;}
 auto runShard=[&](const fs::path& shardOut,int worker){
 const fs::path out=shardOut;
 fs::create_directories(out);
 auto owns=[&](int cid){return channelOwner.at(cid)==worker;};
 fitCache->reset();rootFitCalls=0;
 auto efficiency=table(out/"efficiency.tsv");efficiency<<"CellID\tsample\tpredicted_leaveout_crossings\tunambiguous_observed\tambiguous\tcalib_available\tresponse_selected\tconditional_observed_fraction\n";
 auto closure=table(out/"calibration_closure.tsv");closure<<"CellID\tsample\tchecked_saved_HG\tsaved_HG_mismatches\tmax_HG_error_MeV\tchecked_selected_energy\tselected_energy_mismatches\tmax_selected_error_MeV\n";
 auto fitTable=table(out/"fits.tsv");fitTable<<"CellID\tlayer\tSiPM_um\tsample\tgeometry_condition\tstage\tunit\tplateau_sigmas\tentries\tsum_weights\teffective_hits\tfit_effective_hits\tunderflow_fraction\toverflow_fraction\tfit_low\tfit_high\tLandau_MPV\tMPV_error\tnormalized_MPV\tLandau_width\tGaussian_sigma\tArea\tchi2_ndf\tfit_status\tcovariance_status\tstatus\n";
 std::map<FitKey,Fit> fits;std::map<int,Fit> generatedFits;
 auto stability=table(out/"plateau_stability.tsv");stability<<"CellID\tsample\tgeometry_condition\tMPV_at_configured_plateau\tMPV_at_2sigma\tMPV_at_4sigma\tmax_relative_shift\tall_valid\tstable_within_3percent\n";
 auto spectra=open(out/"spectra.root","RECREATE");
 int summaryCell=0,summarySample=0,summaryMatched=0,summarySiPM=0,summaryFallback=0,summaryValid=0;
 double summaryNormalized=NAN;
 TTree workerSummary("WorkerSummary","Full precision channel fit summary for shard merge");
 workerSummary.SetDirectory(nullptr);
 workerSummary.Branch("CellID",&summaryCell);workerSummary.Branch("sample",&summarySample);
 workerSummary.Branch("matched",&summaryMatched);workerSummary.Branch("sipm_um",&summarySiPM);
 workerSummary.Branch("mip_fallback",&summaryFallback);workerSummary.Branch("valid",&summaryValid);
 workerSummary.Branch("normalized",&summaryNormalized);
 std::array<std::array<std::unique_ptr<TH2D>,2>,2> maps;
 for(int s=0;s<2;++s)for(int matched=0;matched<2;++matched){std::string name=sampleName(s)+"_"+(matched?"matched":"unmatched")+"_mpv_map";maps[s][matched]=std::make_unique<TH2D>(name.c_str(),"HG ADC MPV / current MIP;Local channel (36 chip + channel);Layer;Normalized MPV",216,-.5,215.5,30,-.5,29.5);maps[s][matched]->SetDirectory(nullptr);}
 std::array<long long,2> goodChannels{};std::map<std::tuple<int,int,bool,int>,std::vector<double>> groups;
 auto commonBounds=[&](const Channel& ch,int cid,int matched,int stage,double sigma){
  double low=-INFINITY,high=INFINITY;bool present=false;
  for(int s=0;s<2;++s){if((stage==Truth||stage==Linear)&&!s)continue;
   const auto& sp=ch.spectra[s][matched];
   if(!sp.h[stage]||sp.h[stage]->GetEntries()==0)continue;
   auto range=bounds(sp,pars[s].at(cid),stage,sigma,hgLimit);
   if(!std::isfinite(range.first)||!std::isfinite(range.second))return std::make_pair(double(NAN),double(NAN));
   present=true;low=std::max(low,range.first);high=std::min(high,range.second);
  }
  return present?std::make_pair(low,high):std::make_pair(double(NAN),double(NAN));
 };
 for(auto& [cid,ch]:channels){
  if(!owns(cid))continue;
  for(int s=0;s<2;++s){auto& e=ch.efficiency[s];efficiency<<cid<<'\t'<<sampleName(s)<<'\t'<<e.crossings<<'\t'<<e.observed<<'\t'<<e.ambiguous<<'\t'<<e.calib<<'\t'<<e.response<<'\t'<<(e.crossings?double(e.observed)/e.crossings:NAN)<<'\n';
   auto& c=ch.closure[s];closure<<cid<<'\t'<<sampleName(s)<<'\t'<<c.checkedHG<<'\t'<<c.badHG<<'\t'<<c.maxHG<<'\t'<<c.checkedSelected<<'\t'<<c.badSelected<<'\t'<<c.maxSelected<<'\n';
   for(int matched=0;matched<2;++matched){
    fs::path relative=fs::path(sampleName(s))/(matched?"matched":"unmatched")/energyDir;
    TDirectory* directory=spectra.get();for(const auto& part:relative){auto* next=directory->GetDirectory(part.c_str());if(!next)next=directory->mkdir(part.c_str());directory=next;}
    directory->cd();auto& sp=ch.spectra[s][matched];auto& p=pars[s].at(cid);
    for(int st=0;st<Stages;++st){if(st>=Truth&&!s)continue;if(sp.h[st])sp.h[st]->Write();
     auto range=commonBounds(ch,cid,matched,st,plateau);double scale=stageScale(st,ch.adcRef),normalization=st==SavedMIP?1:st==Truth?.305:st==Linear?digi.at(cid).mip:p.mip;
     Fit f;if(requested.empty()||requested.count(cid))f=fitHistogram(sp.h[st].get(),range.first,range.second,scale,normalization,minFit,maxChi);else f.status="not_requested";
     fits[{cid,s,matched,st,0}]=f;writeFit(fitTable,cid,s,matched?"matched":"unmatched",st,plateau,f);
     if(sp.h[st]){auto curve=fitCurve(f,std::string(sp.h[st]->GetName())+"_langau");if(curve)curve->Write();}
    }
    for(auto& profile:sp.profile)if(profile)profile->Write();
    if(ch.efficiencyLong[s]&&!matched)ch.efficiencyLong[s]->Write();
    auto f=fits.at({cid,s,matched,Temperature,0});std::array<Fit,2> variants;
    for(int k=0;k<2;++k){double sigma=k?4:2;auto range=commonBounds(ch,cid,matched,Temperature,sigma);
     if(!requested.empty()&&!requested.count(cid))variants[k].status="not_requested";
     else variants[k]=std::abs(sigma-plateau)<1e-12?f:fitHistogram(sp.h[Temperature].get(),range.first,range.second,ch.adcRef,p.mip,minFit,maxChi);
     fits[{cid,s,matched,Temperature,k+1}]=variants[k];if(std::abs(sigma-plateau)>1e-12)writeFit(fitTable,cid,s,matched?"matched":"unmatched",Temperature,sigma,variants[k]);}
    bool allValid=f.valid&&variants[0].valid&&variants[1].valid;double shift=allValid?std::max(std::abs(variants[0].mpv/f.mpv-1),std::abs(variants[1].mpv/f.mpv-1)):NAN;
    stability<<cid<<'\t'<<sampleName(s)<<'\t'<<(matched?"matched":"unmatched")<<'\t'<<f.mpv<<'\t'<<variants[0].mpv<<'\t'<<variants[1].mpv<<'\t'<<shift<<'\t'<<allValid<<'\t'<<(allValid&&shift<.03)<<'\n';
    summaryCell=cid;summarySample=s;summaryMatched=matched;summarySiPM=sipm(cid);
    summaryFallback=p.mipFallback;summaryValid=f.valid;summaryNormalized=f.normalized;workerSummary.Fill();
    if(f.valid){maps[s][matched]->SetBinContent(cid%100000/10000*36+cid%100+1,cid/100000+1,f.normalized);maps[s][matched]->SetBinError(cid%100000/10000*36+cid%100+1,cid/100000+1,f.mpvError/p.mip);groups[{s,sipm(cid),p.mipFallback,matched}].push_back(f.normalized);if(matched)++goodChannels[s];}
   }
  }
  if(ch.generatedTruth){spectra->cd();auto* d=spectra->GetDirectory(("mc/generated/"+energyDir).c_str());if(!d){auto* mc=spectra->GetDirectory("mc");auto* gen=mc->GetDirectory("generated");if(!gen)gen=mc->mkdir("generated");d=gen->mkdir(energyDir.c_str());}d->cd();ch.generatedTruth->Write();ch.generatedTruthDetected->Write();Fit f;
   if(requested.empty()||requested.count(cid))f=fitHistogram(ch.generatedTruth.get(),0,8*.305,.305,.305,minFit,maxChi);else f.status="not_requested";generatedFits[cid]=f;
   auto curve=fitCurve(f,std::string(ch.generatedTruth->GetName())+"_langau");if(curve)curve->Write();
   // Generated truth is a separate cohort, not geometry-matched reconstructed hits.
   fitTable<<cid<<'\t'<<cid/100000<<'\t'<<sipm(cid)<<"\tmc\tgenerated_independent\ttruth_mev\tMeV\tnan\t"<<f.entries<<'\t'<<f.sumw<<'\t'<<f.neff<<'\t'<<f.fitNeff<<'\t'<<f.underflow<<'\t'<<f.overflow<<'\t'<<f.low<<'\t'<<f.high<<'\t'<<f.mpv<<'\t'<<f.mpvError<<'\t'<<f.normalized<<'\t'<<f.landau<<'\t'<<f.gaussian<<'\t'<<f.area<<'\t'<<f.chi2ndf<<'\t'<<f.rootStatus<<'\t'<<f.covariance<<'\t'<<f.status<<'\n';
  }
  if(cid%100==0)std::cout<<"Fitted through CellID "<<cid<<'\n';
 }
 auto comparison=table(out/"comparison.tsv");comparison<<"CellID\tSiPM_um\tdata_MIP_fallback\tmc_MIP_fallback\tgeometry_condition\tdata_temperature_MPV_over_MIP\tmc_temperature_MPV_over_MIP\tdata_ADC_MPV_over_mc_ADC_MPV\tdata_calibrated_MPV_over_mc_calibrated_MPV\tdata_saved_MIP_MPV\tmc_saved_MIP_MPV\tmc_ADC_MPV_over_digi_MIP\tgenerated_truth_MPV_over_0305\tdetected_crossing_truth_MPV_over_0305\tdata_fit_status\tmc_fit_status\tboth_valid\n";
 std::map<std::tuple<int,bool,bool,int>,std::array<std::vector<double>,3>> commonGroups;
 std::array<long long,2> commonValid{};
 for(auto& [cid,ch]:channels)if(owns(cid))for(int matched=0;matched<2;++matched){const auto& d=fits.at({cid,0,matched,Temperature,0});const auto& m=fits.at({cid,1,matched,Temperature,0});const auto& ds=fits.at({cid,0,matched,SavedMIP,0});const auto& ms=fits.at({cid,1,matched,SavedMIP,0});const auto& truth=fits.at({cid,1,matched,Truth,0});comparison<<cid<<'\t'<<sipm(cid)<<'\t'<<pars[0].at(cid).mipFallback<<'\t'<<pars[1].at(cid).mipFallback<<'\t'<<(matched?"matched":"unmatched")<<'\t'<<d.normalized<<'\t'<<m.normalized<<'\t'<<(d.valid&&m.valid?d.mpv/m.mpv:NAN)<<'\t'<<(d.valid&&m.valid?d.normalized/m.normalized:NAN)<<'\t'<<(ds.valid?ds.normalized:NAN)<<'\t'<<(ms.valid?ms.normalized:NAN)<<'\t'<<(m.valid?m.mpv/digi.at(cid).mip:NAN)<<'\t'<<(generatedFits.count(cid)&&generatedFits.at(cid).valid?generatedFits.at(cid).normalized:NAN)<<'\t'<<(truth.valid?truth.normalized:NAN)<<'\t'<<d.status<<'\t'<<m.status<<'\t'<<(d.valid&&m.valid)<<'\n';
  if(d.valid&&m.valid){++commonValid[matched];auto& v=commonGroups[{sipm(cid),pars[0].at(cid).mipFallback,pars[1].at(cid).mipFallback,matched}];v[0].push_back(d.normalized);v[1].push_back(m.normalized);v[2].push_back(d.normalized/m.normalized);}
 }
 auto commonGroup=table(out/"common_group_summary.tsv");commonGroup<<"SiPM_um\tdata_MIP_adoption\tmc_MIP_adoption\tgeometry_condition\tcommon_valid_channel_count\tmedian_data_temperature_MPV_over_MIP\tmedian_mc_temperature_MPV_over_MIP\tmedian_channel_data_over_mc_ratio\n";
 for(auto& [key,v]:commonGroups)commonGroup<<std::get<0>(key)<<'\t'<<(std::get<1>(key)?"fallback":"individual")<<'\t'<<(std::get<2>(key)?"fallback":"individual")<<'\t'<<(std::get<3>(key)?"matched":"unmatched")<<'\t'<<v[0].size()<<'\t'<<median(v[0])<<'\t'<<median(v[1])<<'\t'<<median(v[2])<<'\n';
 auto fitCounts=table(out/"fit_status_counts.tsv");fitCounts<<"sample\tgeometry_condition\tstage\tstatus\tchannels\n";
 std::map<std::tuple<int,int,int,std::string>,long long> statusCounts;
 for(auto& [key,f]:fits)if(std::get<4>(key)==0)++statusCounts[{std::get<1>(key),std::get<2>(key),std::get<3>(key),f.status}];
 for(auto& [key,n]:statusCounts)fitCounts<<sampleName(std::get<0>(key))<<'\t'<<(std::get<1>(key)?"matched":"unmatched")<<'\t'<<stageName[std::get<2>(key)]<<'\t'<<std::get<3>(key)<<'\t'<<n<<'\n';
 auto group=table(out/"group_summary.tsv");group<<"sample\tSiPM_um\tMIP_adoption\tgeometry_condition\tvalid_channel_count\tmedian_temperature_MPV_over_MIP\n";
 for(auto& [key,v]:groups)group<<sampleName(std::get<0>(key))<<'\t'<<std::get<1>(key)<<'\t'<<(std::get<2>(key)?"fallback":"individual")<<'\t'<<(std::get<3>(key)?"matched":"unmatched")<<'\t'<<v.size()<<'\t'<<median(v)<<'\n';
 std::map<DependencyGraphKey,std::array<std::unique_ptr<TGraphAsymmErrors>,2>> dependencyGraphs;
 long long dependencyBinCount=0,dependencyValidFits=0;
 if(dependencyFits){
  const auto& widths=dependencyWidths;
  auto allowDependency=[&](int cid,int matched){return (!requested.empty()&&requested.count(cid))||fits.at({cid,0,matched,Temperature,0}).valid||fits.at({cid,1,matched,Temperature,0}).valid;};
  auto detail=table(out/"dependency_fits.tsv");
  detail<<"CellID\tsample\tgeometry_condition\taxis\tbin_index\tbin_low\tbin_high\tmean_axis_value\tADC_Landau_MPV\tMPV_error_ADC\tMPV_over_current_MIP\teffective_hits\tfit_effective_hits\tfit_low_ADC\tfit_high_ADC\tchi2_ndf\troot_fit_status\tcovariance_status\tstatus\n";
  for(auto& [key,samples]:dependency){
   const int cid=std::get<0>(key),matched=std::get<1>(key),axis=std::get<2>(key),bin=std::get<3>(key);
   if(!owns(cid)||!allowDependency(cid,matched))continue;
   ++dependencyBinCount;
   double low=-INFINITY,high=INFINITY;bool unknown=false;
   for(int s=0;s<2;++s)if(samples[s].h){const auto& p=pars[s].at(cid);if(!std::isfinite(p.threshold)||!std::isfinite(p.sigma)){unknown=true;continue;}low=std::max(low,(p.threshold+plateau*p.sigma)*samples[s].lowerFactor);high=std::min(high,(hgLimit-p.ph)*samples[s].upperFactor);}
   if(unknown)low=high=NAN;
   for(int s=0;s<2;++s){
    auto& sp=samples[s];if(!sp.h)continue;auto f=fitHistogram(sp.h.get(),std::max(0.,low),high,channels.at(cid).adcRef,pars[s].at(cid).mip,minFit,maxChi);
    detail<<cid<<'\t'<<sampleName(s)<<'\t'<<(matched?"matched":"unmatched")<<'\t'<<dependencyAxis[axis]<<'\t'<<bin<<'\t'<<bin*widths[axis]<<'\t'<<(bin+1)*widths[axis]<<'\t'<<(sp.sumw?sp.sumx/sp.sumw:NAN)<<'\t'<<f.mpv<<'\t'<<f.mpvError<<'\t'<<f.normalized<<'\t'<<f.neff<<'\t'<<f.fitNeff<<'\t'<<f.low<<'\t'<<f.high<<'\t'<<f.chi2ndf<<'\t'<<f.rootStatus<<'\t'<<f.covariance<<'\t'<<f.status<<'\n';
    TDirectory* directory=spectra.get();for(const auto& part:fs::path(sampleName(s))/(matched?"matched":"unmatched")/energyDir/"dependencies"){auto* next=directory->GetDirectory(part.c_str());if(!next)next=directory->mkdir(part.c_str());directory=next;}directory->cd();sp.h->Write();auto curve=fitCurve(f,std::string(sp.h->GetName())+"_langau");if(curve)curve->Write();
    if(f.valid){++dependencyValidFits;auto& graph=dependencyGraphs[{cid,matched,axis}][s];if(!graph){graph=std::make_unique<TGraphAsymmErrors>();graph->SetName((sampleName(s)+"_cell"+std::to_string(cid)+"_"+(matched?"matched_":"unmatched_")+dependencyAxis[axis]+"_mpv").c_str());graph->SetTitle(("Binwise Landau MPV;"+dependencyAxisLabel[axis]+";MPV / current MIP").c_str());}int point=graph->GetN();const double x=sp.sumx/sp.sumw,yError=f.mpvError/pars[s].at(cid).mip;
     graph->SetPoint(point,x,f.normalized);
     const double xLow=axis==3?0:std::max(0.,x-bin*widths[axis]),xHigh=axis==3?0:std::max(0.,(bin+1)*widths[axis]-x);
     graph->SetPointError(point,xLow,xHigh,yError,yError);
}
    sp.h.reset();
   }
  }
  // Graphs keep mean positions; asymmetric X bars span the original fit bins.
  // Run bars remain zero; exact bin extents and counts are retained in TSV.
  for(auto& [key,samples]:dependencyGraphs)for(auto& graph:samples)if(graph){TDirectory* directory=spectra.get();for(const auto& part:fs::path("comparison")/(std::get<1>(key)?"matched":"unmatched")/energyDir/"dependencies"){auto* next=directory->GetDirectory(part.c_str());if(!next)next=directory->mkdir(part.c_str());directory=next;}directory->cd();graph->Write();}
 }
 spectra->cd();for(auto& sample:maps)for(auto& m:sample)m->Write();
 if(workerCount>1)workerSummary.Write();
 TNamed("response_definition","HG pedestal subtraction;2 mm path normalization;C(T)=1-(T-20)*alpha;saved HG energy/0.305;matching uses run/channel/position/slopes only").Write();
 TNamed("efficiency_scope","Conditional on reconstructed leave-one-layer-out tracks;not unconditional hardware acquisition efficiency").Write();
 TNamed("truth_scope","Generated truth uses primary strip chord;invalid/noncrossed primary chords excluded;reconstructed detected truth is a separate conditional cohort").Write();
 TNamed("histogram_normalization","PNG histogram fractions use the sum including underflow and overflow; displayed bins do not include flow").Write();
 TNamed("dependency_scope","Binwise MPV fits use temperature-corrected HG;only inclusive-fit-valid channels or explicitly requested channels are revisited;mean profiles are distinct observables").Write();
 TParameter<double>("geometry_common_target",commonTarget).Write();spectra->Close();
 bool drawChannels=config.integer("channel-plots",1);
 TCanvas canvas("muon_response","",mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.SetBatch(true);canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);long long rendered=0;
 if(drawChannels)for(auto& [cid,ch]:channels){
  if(!owns(cid))continue;
  if(!requested.empty()&&!requested.count(cid))continue;
  bool valid=false;for(int s=0;s<2;++s)for(int m=0;m<2;++m)valid=valid||fits.at({cid,s,m,Temperature,0}).valid;
  if(requested.empty()&&!valid)continue;
  for(int matched=0;matched<2;++matched){if(!ch.spectra[0][matched].h[0]&&!ch.spectra[1][matched].h[0])continue;
   fs::path dest=out/"comparison"/(matched?"matched":"unmatched")/energyDir/("layer"+std::to_string(cid/100000))/("chip"+std::to_string(cid%100000/10000));fs::create_directories(dest);canvas.Clear();canvas.Divide(3,2);std::vector<std::unique_ptr<TH1D>> drawn;
   for(int stage=0;stage<4;++stage){
    std::array<mmr_plot::FitInfo,2> info;std::array<std::unique_ptr<TF1>,2> curves;
    for(int sample=0;sample<2;++sample){const auto& fit=fits.at({cid,sample,matched,stage,0});info[sample]=plotFitInfo(fit);curves[sample]=fitCurve(fit,"draw_fit_"+std::to_string(stage)+"_"+std::to_string(sample));}
    mmr_plot::drawResponsePanel(canvas,stage+1,ch.spectra[0][matched].h[stage].get(),ch.spectra[1][matched].h[stage].get(),"Cell "+std::to_string(cid)+" | "+stageTitle[stage],units[stage],drawn,&info[0],&info[1],curves[0].get(),curves[1].get());
   }
   overlay(canvas,5,ch.generatedTruth.get(),ch.spectra[1][matched].h[Truth].get(),"MC generated vs reconstructed detected truth","MeV",drawn,"Generated truth","Detected crossings");canvas.cd(5);annotation(.92,"Independent truth uses primary strip chord");
   overlay(canvas,6,ch.spectra[1][matched].h[Linear].get(),ch.spectra[1][matched].h[Temperature].get(),"MC linear prediction vs final HG","ADC",drawn,"Linear prediction","Final HG (20 C)");canvas.cd(6);annotation(.92,"Same reconstructed crossing cohort");
   fastPng(canvas,(dest/("cell"+std::to_string(cid)+"_response.png")).c_str());canvas.Clear();drawn.clear();
   profiles(canvas,ch,matched,dest/("cell"+std::to_string(cid)+"_mean_profiles.png"));
   if(dependencyFits)mpvProfiles(canvas,cid,matched,dependencyGraphs,dest/("cell"+std::to_string(cid)+"_profiles.png"));
   ++rendered;
  }
 }
 canvas.Clear();canvas.SetCanvasSize(mmr_plot::canvasWidth,mmr_plot::canvasHeight);canvas.Divide(2,2);for(int s=0;s<2;++s)for(int m=0;m<2;++m){canvas.cd(1+2*s+m);auto* h=maps[s][m].get();h->SetMinimum(.5);h->SetMaximum(1.5);h->SetStats(false);h->SetTitle((sampleName(s)+", "+(m?"matched":"unmatched")+";Local channel (36 chip + channel);Layer;ADC MPV / MIP").c_str());h->Draw("COLZ");}
 fs::create_directories(out/"comparison");fastPng(canvas,(out/"comparison/muon_mip_channel_map.png").c_str());canvas.Clear();
 long long badSavedHG=0,badSavedSelected=0;for(const auto& [cid,ch]:channels)if(owns(cid))for(const auto& c:ch.closure){badSavedHG+=c.badHG;badSavedSelected+=c.badSelected;}
 auto status=table(out/"analysis_status.tsv");status<<"geometry_common_bins\tdata_matched_valid_channels\tmc_matched_valid_channels\tcommon_matched_valid_channels\tcommon_unmatched_valid_channels\tchannel_conditions_rendered\tcrossings_read\tgenerated_truth_valid\tdependency_fits_enabled\tdependency_bins\tdependency_valid_fits\tsaved_HG_mismatches\tselected_energy_mismatches\tdiagnostic_status\n"<<commonBins<<'\t'<<goodChannels[0]<<'\t'<<goodChannels[1]<<'\t'<<commonValid[1]<<'\t'<<commonValid[0]<<'\t'<<rendered<<'\t'<<crossingEntries<<'\t'<<generatedValid<<'\t'<<dependencyFits<<'\t'<<dependencyBinCount<<'\t'<<dependencyValidFits<<'\t'<<badSavedHG<<'\t'<<badSavedSelected<<'\t'<<(badSavedHG||badSavedSelected?"saved_calibration_mismatch":commonValid[1]?"paired_fits_available":"no_common_valid_matched_fit")<<'\n';
 auto performance=table(out/"performance.tsv");
 performance<<"worker_index\tfit_requests\tcache_hits\tcache_misses\troot_fit_calls\tcache_evictions\tfit_cache_bytes\tdependency_prepared_bins\n";
 const auto& cacheStats=fitCache->stats();
 performance<<worker<<'\t'<<cacheStats.requests<<'\t'<<cacheStats.hits<<'\t'<<cacheStats.misses<<'\t'<<rootFitCalls<<'\t'<<cacheStats.evictions<<'\t'<<fitCache->bytes()<<'\t'<<dependency.size()<<'\n';
 std::cout<<"Muon response fit worker "<<worker<<" complete: "<<out<<", matched fit channels data/MC "<<goodChannels[0]<<'/'<<goodChannels[1]<<std::endl;
 };
 if(workerCount==1){runShard(out,0);return 0;}
 const fs::path binaryDir=fs::read_symlink("/proc/self/exe").parent_path();
 const fs::path project=binaryDir.parent_path().parent_path().parent_path();
 const std::string configuredScratch=config.text("scratch-root");
 const fs::path scratchRoot=configuredScratch.empty()?project/".scratch/muon_mip_response_runs":fs::absolute(configuredScratch);
 fs::create_directories(scratchRoot);
 std::string pattern=(scratchRoot/("fit_workers_"+std::to_string(getpid())+"_XXXXXX")).string();
 std::vector<char> mutablePattern(pattern.begin(),pattern.end());mutablePattern.push_back('\0');
 char* created=mkdtemp(mutablePattern.data());
 need(created!=nullptr,"Cannot create isolated worker temporary directory in "+scratchRoot.string());
 const fs::path parts=created;
 std::cout<<"Fit worker temporary outputs: "<<parts<<std::endl;
 std::vector<pid_t> children;
 for(int worker=0;worker<workerCount;++worker){
  std::cout.flush();std::cerr.flush();
  const pid_t pid=fork();
  if(pid<0){
   for(pid_t child:children){int status=0;waitpid(child,&status,0);}
   throw std::runtime_error("Cannot create isolated fit worker");
  }
  if(pid==0){
   try{runShard(parts/("worker_"+std::to_string(worker)),worker);std::cout.flush();std::cerr.flush();_exit(0);}
   catch(const std::exception& e){std::cerr<<"Fit worker "<<worker<<": "<<e.what()<<std::endl;_exit(1);}
  }
  children.push_back(pid);
 }
 bool failed=false;
 for(pid_t child:children){int status=0;need(waitpid(child,&status,0)==child,"Cannot wait for fit worker");failed=failed||!WIFEXITED(status)||WEXITSTATUS(status)!=0;}
 need(!failed,"At least one fit worker failed; partial worker outputs are preserved in "+parts.string());
 const fs::path merger=fs::read_symlink("/proc/self/exe").parent_path()/"MergeMuonShards";
 need(fs::is_regular_file(merger),"Missing fit shard merger: "+merger.string());
 const pid_t mergerPid=fork();need(mergerPid>=0,"Cannot create fit shard merger");
 if(mergerPid==0){execl(merger.c_str(),merger.c_str(),parts.c_str(),argv[2],argv[3],static_cast<char*>(nullptr));std::cerr<<"Cannot execute fit shard merger"<<std::endl;_exit(127);}
 int mergeStatus=0;need(waitpid(mergerPid,&mergeStatus,0)==mergerPid&&WIFEXITED(mergeStatus)&&WEXITSTATUS(mergeStatus)==0,"Fit shard merge failed; partial worker outputs are preserved in "+parts.string());
 // Remove only this process's unique temporary subtree after the final ROOT,
 // scientific tables, full-precision summaries, and PNGs have been merged.
 fs::remove_all(parts);
 std::error_code scratchCleanup;
 fs::remove(scratchRoot,scratchCleanup); // succeeds only if the parent is empty
 std::cout<<"Muon response complete: "<<out<<", isolated workers="<<workerCount<<std::endl;
 return 0;
}catch(const std::exception& e){std::cerr<<"AnalyzeMuonMIP: "<<e.what()<<std::endl;return 1;}
