#pragma once
#include "../channel_response/muon_mip_response/CurrentCalibration.hh"
#include "../channel_response/muon_mip_response/TrackGeometry.hh"
#include <TH1D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TLeaf.h>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <numeric>
#include <climits>
namespace mip {
namespace fs=std::filesystem;
using mmr::require;
constexpr double nan=std::numeric_limits<double>::quiet_NaN();
inline std::set<int> selection(const std::string&s) {
 std::set<int> out;if(s=="all") {for(int id:mmr::channels())out.insert(id);return out;}
 for(auto&v:mmr::split(s,',')){size_t n;int id=std::stoi(v,&n);require(n==v.size()&&mmr::connected(id)&&mmr::physical(id)==id,"Invalid physical CellID: "+v);require(out.insert(id).second,"Duplicate selected CellID");}
 require(!out.empty(),"Empty channel selection");return out;
}
inline std::string channelName(int id){return "cell"+std::to_string(id);}
inline TDirectory* directory(TDirectory&d,const std::string&s){auto*p=d.GetDirectory(s.c_str());if(!p)p=d.mkdir(s.c_str());require(p,"Cannot create ROOT directory "+s);return p;}
inline void fresh(const std::string&p){require(!fs::exists(p),"Refusing to overwrite "+p);}
struct Meta {
 int CellID=0,selected=0,old_fallback=1,ped_present=0,threshold_fallback=1;
 double pedestal=nan,threshold=nan,threshold_sigma=nan,old_mip=nan,old_raw_mip=nan;
 Long64_t predicted=0,observed=0,unique=0,duplicate=0,usable=0,switched=0,forced4000=0,saturated=0,lg_saturated=0,invalid_adc=0,missing_temperature=0,invalid_temperature=0,gain0=0;
 void branches(TTree&t,bool read=false){
#define B(v) if(read)mmr::bind(&t,#v,&v);else t.Branch(#v,&v)
 B(CellID);B(selected);B(old_fallback);B(ped_present);B(threshold_fallback);B(pedestal);B(threshold);B(threshold_sigma);B(old_mip);B(old_raw_mip);B(predicted);B(observed);B(unique);B(duplicate);B(usable);B(switched);B(forced4000);B(saturated);B(lg_saturated);B(invalid_adc);B(missing_temperature);B(invalid_temperature);B(gain0);
#undef B
 }
};
inline std::map<int,Meta> readMeta(TFile&f){std::map<int,Meta> a;auto*t=mmr::tree(f,"Channels");Meta m;m.branches(*t,true);for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);require(a.emplace(m.CellID,m).second,"Duplicate channel metadata");}t->ResetBranchAddresses();return a;}
struct Condition {
 int CellID=0,run=0,mode=0,bin=0;double factor=1,lo=1,hi=1;Long64_t predicted=0;
 std::string name()const{return "r"+std::to_string(run)+"_m"+std::to_string(mode)+"_b"+std::to_string(bin);}
 void branches(TTree&t,bool read=false){
#define B(v) if(read)mmr::bind(&t,#v,&v);else t.Branch(#v,&v)
 B(CellID);B(run);B(mode);B(bin);B(factor);B(lo);B(hi);B(predicted);
#undef B
 }
};
inline std::map<int,std::vector<Condition>> readConditions(TFile&f){std::map<int,std::vector<Condition>>a;auto*t=mmr::tree(f,"Conditions");Condition c;c.branches(*t,true);for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);a[c.CellID].push_back(c);}t->ResetBranchAddresses();return a;}
inline std::unique_ptr<TH1D> hist(const std::string&name,int bins,double low,double high){auto h=std::make_unique<TH1D>(name.c_str(),name.c_str(),bins,low,high);h->SetDirectory(nullptr);return h;}
inline std::unique_ptr<TH1D> getHist(TFile&f,const std::string&path){auto*h=dynamic_cast<TH1D*>(f.Get(path.c_str()));require(h,"Missing histogram "+path);h->SetDirectory(nullptr);return std::unique_ptr<TH1D>(h);}
inline std::unique_ptr<TH1D> rebinForPlot(const TH1D&h,int group){
 int n=h.GetNbinsX()/group;auto out=hist(std::string(h.GetName())+"_display",n,h.GetXaxis()->GetXmin(),h.GetXaxis()->GetXmin()+n*group*h.GetBinWidth(1));
 for(int i=0;i<=h.GetNbinsX()+1;++i){int b=i==0?0:std::min(n+1,(i-1)/group+1);out->AddBinContent(b,h.GetBinContent(i));}out->SetEntries(h.GetEntries());return out;
}
// Decision codes: 0 undecided, 1 direct measurement accepted, 2 rejected.
struct Result {
 int CellID=0,decision=0,status=-1,covariance=0,NDF=0,chnEntries=0,conditions=0,boundary=0,variant=0,model_outlier_hits=0;
 double LandauMPV=nan,LandauWidth=nan,GauSigma=nan,LandauMPV_error=nan,LandauWidth_error=nan,GauSigma_error=nan,ChiSquare=nan;
 double convolution_peak=nan,fit_low=nan,fit_high=nan,plateau=nan,relative_error=nan,stability=nan,correlation=nan,objective=nan,edm=nan,bin_width=nan;
 int boundary_mask=0,model=0;
 double raw_fit_low=0;
 double cov[9]={};
 std::string reason="not_selected";
 void branches(TTree&t,bool read=false){
#define B(v) if(read)mmr::bind(&t,#v,&v);else t.Branch(#v,&v)
 B(CellID);B(decision);B(status);B(covariance);B(NDF);B(chnEntries);B(conditions);B(boundary);B(variant);B(model_outlier_hits);B(LandauMPV);B(LandauWidth);B(GauSigma);B(LandauMPV_error);B(LandauWidth_error);B(GauSigma_error);B(ChiSquare);B(convolution_peak);B(fit_low);B(fit_high);B(plateau);B(relative_error);B(stability);B(correlation);B(objective);B(edm);B(bin_width);
#undef B
 // Optional fields preserve read compatibility with earlier saved fits.
#define O(v) if(!read)t.Branch(#v,&v);else if(t.GetBranch(#v))mmr::bind(&t,#v,&v)
 O(boundary_mask);O(model);O(raw_fit_low);
#undef O
 if(!read){t.Branch("cov",cov,"cov[9]/D");t.Branch("reason",&reason);}else mmr::bind(&t,"cov",cov);
 }
};
inline void resultHeader(std::ostream&o){o<<"CellID\tvariant\tdecision\treason\tentries\tconditions\tstatus\tcovariance\tMPV\tMPV_error\tLandauWidth\tLandauWidth_error\tGauSigma\tGauSigma_error\tconvolution_peak\tchi2\tNDF\tchi2_ndf\tfit_low\tfit_high\tbin_width\tplateau\trelative_error\tstability\tmax_MPV_correlation\tboundary\tdeviance\tedm\tmodel_outlier_hits\tboundary_mask\tmodel\traw_fit_low\n";}
inline void resultRow(std::ostream&o,const Result&r){o<<std::setprecision(12)<<r.CellID<<'\t'<<r.variant<<'\t'<<r.decision<<'\t'<<r.reason<<'\t'<<r.chnEntries<<'\t'<<r.conditions<<'\t'<<r.status<<'\t'<<r.covariance<<'\t'<<r.LandauMPV<<'\t'<<r.LandauMPV_error<<'\t'<<r.LandauWidth<<'\t'<<r.LandauWidth_error<<'\t'<<r.GauSigma<<'\t'<<r.GauSigma_error<<'\t'<<r.convolution_peak<<'\t'<<r.ChiSquare<<'\t'<<r.NDF<<'\t'<<(r.NDF>0?r.ChiSquare/r.NDF:nan)<<'\t'<<r.fit_low<<'\t'<<r.fit_high<<'\t'<<r.bin_width<<'\t'<<r.plateau<<'\t'<<r.relative_error<<'\t'<<r.stability<<'\t'<<r.correlation<<'\t'<<r.boundary<<'\t'<<r.objective<<'\t'<<r.edm<<'\t'<<r.model_outlier_hits<<'\t'<<r.boundary_mask<<'\t'<<r.model<<'\t'<<r.raw_fit_low<<'\n';}
inline std::map<int,Result> readResults(TFile&f){auto*t=mmr::tree(f,"Diagnostics");Result r;r.branches(*t,true);std::string*reason=nullptr;mmr::bind(t,"reason",&reason);std::map<int,Result>a;for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);r.reason=*reason;require(a.emplace(r.CellID,r).second,"Duplicate result CellID");}t->ResetBranchAddresses();delete reason;return a;}
}
