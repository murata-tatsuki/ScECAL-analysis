#pragma once
#include "../muon_mip_response/CurrentCalibration.hh"
#include "../muon_mip_response/TrackGeometry.hh"
#include <TH1D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <cstdint>
namespace spr {
using mmr::require; namespace fs=std::filesystem;
constexpr double qnan=std::numeric_limits<double>::quiet_NaN();
constexpr int schema=1;
inline std::ofstream table(const fs::path&p){fs::create_directories(p.parent_path());std::ofstream o(p);require(bool(o),"Cannot write "+p.string());o<<std::setprecision(12);return o;}
inline std::set<int> selection(const std::string&s){std::set<int>v;if(s=="all"){for(int id:mmr::channels())v.insert(id);}else for(auto&a:mmr::split(s,',')){int id=std::stoi(a);require(mmr::connected(id)&&mmr::physical(id)==id,"Invalid channel "+a);v.insert(id);}require(!v.empty(),"Empty channel selection");return v;}
inline uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
inline uint64_t identityHash(int run,int time,int event){return mix(uint32_t(run))^mix(uint32_t(time)+0x100000000ULL)^mix(uint32_t(event)+0x200000000ULL);}
inline std::string sample(int s){return s?"mc":"data";}
inline std::string condition(int c){return c==0?"unadjusted":c==1?"position_matched":"data_mc_matched";}
inline const char* stage(int s){return s==0?"pedestal_adc":s==1?"path_adc":"temperature_adc";}
inline std::string base(int s,int c,const mmr::Config&cfg){return sample(s)+"/"+condition(c)+"/"+cfg.text("energy","100")+"GeV";}
inline std::string channelDir(int id){return "layer"+std::to_string(id/100000)+"/chip"+std::to_string(id%100000/10000);}
inline TDirectory* directory(TDirectory&d,const std::string&path){auto*p=&d;for(auto&s:mmr::split(path,'/')){auto*q=p->GetDirectory(s.c_str());if(!q)q=p->mkdir(s.c_str());require(q,"Cannot create ROOT directory");p=q;}return p;}
inline std::unique_ptr<TH1D> hist(const std::string&name){auto h=std::make_unique<TH1D>(name.c_str(),";ADC;Absolute hits / 8 ADC",400,-200,3000);h->SetDirectory(nullptr);h->Sumw2();return h;}
struct Crossing {
 int sample=0,run=0,source_run=0,event=0,event_time=0,cycle=0,cellid=0,full_cellid=-1,observed=0,multiplicity=0,usable=0,temperature_state=0,calib_paired=0,gain=0;
 Long64_t entry=0;double raw_hg=qnan,raw_lg=qnan,temperature=qnan,saved_temperature=qnan,saved_energy=qnan,sx=qnan,sy=qnan,length=qnan,longitudinal=qnan,transverse=qnan,track_rms=qnan,raw_adc=qnan,path_adc=qnan,current_adc=qnan;
 void branches(TTree&t,bool reading=false){
#define B(x) if(reading)mmr::bind(&t,#x,&x);else t.Branch(#x,&x)
 B(sample);B(run);B(source_run);B(event);B(event_time);B(cycle);B(cellid);B(full_cellid);B(observed);B(multiplicity);B(usable);B(temperature_state);B(calib_paired);B(gain);B(entry);B(raw_hg);B(raw_lg);B(temperature);B(saved_temperature);B(saved_energy);B(sx);B(sy);B(length);B(longitudinal);B(transverse);B(track_rms);B(raw_adc);B(path_adc);B(current_adc);
#undef B
 }
};
struct Key {
 int sample=0,condition=0,split=0,cellid=0,axis=0,bin=0,stage=0;
 auto tie()const{return std::tie(sample,condition,split,cellid,axis,bin,stage);}
 bool operator<(const Key&v)const{return tie()<v.tie();}
 std::string name()const{return "s"+std::to_string(sample)+"c"+std::to_string(condition)+"d"+std::to_string(split)+"i"+std::to_string(cellid)+"a"+std::to_string(axis)+"b"+std::to_string(bin)+"t"+std::to_string(stage);}
};
inline int bins(int axis,const mmr::Config&c){return c.integer(axis?"v-bins":"u-bins",3);}
inline double half(int axis){return axis?2.5:22.5;}
inline int posBin(double x,int axis,const mmr::Config&c){int b=int(std::floor((x+half(axis))/(2*half(axis))*bins(axis,c)));return b>=0&&b<bins(axis,c)?b:-1;}
inline double center(int b,int a,const mmr::Config&c){return -half(a)+(b+.5)*2*half(a)/bins(a,c);}
struct Bin {
 Key k;Long64_t predicted=0,observed=0,unique=0,duplicate=0,usable=0,forced=0,saturated=0,switched=0,invalid_temperature=0,negative=0,zero=0;
 double sum=0,sum2=0,fmin=INFINITY,fmax=0,threshold=qnan,sigma=qnan,pedestal=qnan;
 int threshold_fallback=1;
 void branches(TTree&t,bool reading=false){
#define B(x) if(reading)mmr::bind(&t,#x,&x);else t.Branch(#x,&x)
 B(k.sample);B(k.condition);B(k.split);B(k.cellid);B(k.axis);B(k.bin);B(k.stage);B(predicted);B(observed);B(unique);B(duplicate);B(usable);B(forced);B(saturated);B(switched);B(invalid_temperature);B(negative);B(zero);B(sum);B(sum2);B(fmin);B(fmax);B(threshold);if(reading)mmr::bind(&t,t.GetBranch("threshold_sigma")?"threshold_sigma":"sigma",&sigma);else t.Branch("threshold_sigma",&sigma);B(pedestal);B(threshold_fallback);
#undef B
 }
};
struct Fit {
 Bin b;int status=-999,covariance=0,ndf=0,valid=0;double low=qnan,high=qnan,nfit=0,mpv=qnan,mpv_error=qnan,width=qnan,width_error=qnan,sigma=qnan,sigma_error=qnan,area=qnan,chi2=qnan,relative=qnan,relative_error=qnan,stability=qnan;
 std::string reason="unmeasured";
 void branches(TTree&t,bool reading=false){b.branches(t,reading);
#define B(x) if(reading)mmr::bind(&t,#x,&x);else t.Branch(#x,&x)
 B(status);B(covariance);B(ndf);B(valid);B(low);B(high);B(nfit);B(mpv);B(mpv_error);B(width);B(width_error);if(reading)mmr::bind(&t,"gaussian_sigma",&sigma);else t.Branch("gaussian_sigma",&sigma);B(sigma_error);B(area);B(chi2);B(relative);B(relative_error);B(stability);
#undef B
 if(!reading)t.Branch("reason",&reason);
 }
};
inline std::map<Key,Fit> readFits(TFile&f){auto*t=mmr::tree(f,"Fits");Fit x;x.branches(*t,true);std::string*r=nullptr;mmr::bind(t,"reason",&r);std::map<Key,Fit>v;for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);x.reason=*r;require(v.emplace(x.b.k,x).second,"Duplicate fit key");}t->ResetBranchAddresses();delete r;return v;}
inline double mean(const Bin&b){return b.usable?b.sum/b.usable:qnan;}
inline double meanError(const Bin&b){return b.usable>1?std::sqrt(std::max(0.,(b.sum2-b.sum*b.sum/b.usable)/(b.usable*(b.usable-1)))):qnan;}
}
