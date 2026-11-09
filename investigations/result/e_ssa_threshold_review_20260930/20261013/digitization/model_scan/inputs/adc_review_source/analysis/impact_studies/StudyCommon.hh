#pragma once
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;
inline void require(bool ok, const std::string& message) { if (!ok) throw std::runtime_error(message); }
inline std::unique_ptr<TFile> openFile(const std::string& path) {
 auto f=std::unique_ptr<TFile>(TFile::Open(path.c_str())); require(f && !f->IsZombie(), "Cannot open "+path); return f;
}
template<class T> T* object(TFile& f,const std::string& key) { auto*p=dynamic_cast<T*>(f.Get(key.c_str())); require(p,"Missing "+key+" in "+f.GetName());return p; }
template<class T> void bind(TTree* t,const char* key,T* address) { t->SetBranchStatus(key,1);require(t->SetBranchAddress(key,address)>=0,std::string("Invalid branch ")+key); }
inline int physical(int id) { return id/100000*100000+(id%100000/10000)*10000+id%100; }
inline bool connected(int id) { return id>=0 && id/100000<30 && id%100000/10000<6 && id%100<36 && !(id%100000/10000==5 && id%100>29); }
inline int group(int id) { int layer=id/100000;return layer>3 && layer<28 ? 1 : 0; }
struct Stats {
 long long n=0; double mean=0,m2=0;
 void add(double x) { require(std::isfinite(x),"Nonfinite input");++n;double d=x-mean;mean+=d/n;m2+=d*(x-mean); }
 double rms()const{return n?std::sqrt(std::max(0.,m2/n)):0;}
 double sem()const{return n>1?std::sqrt(std::max(0.,m2/(n*(n-1.)))):0;}
};
struct Parameter { double gain=1,intercept=0,mip=1,ph=0,pl=0; bool mipFallback=false,hlFallback=false,ped=false; double offset()const{return .305*intercept/mip;} };
// Reproduce InterFactorExtract_simulation/MIPFactorExtract_simulation, including
// their quality cuts, two SiPM groups, and sentinel=1 fallback convention.
inline std::map<int,Parameter> parameters(const std::string& hl,const std::string& mip,const std::string& ped) {
 std::map<int,Parameter> a;
 for(int l=0;l<30;++l)for(int c=0;c<6;++c)for(int ch=0;ch<(c==5?30:36);++ch)a[l*100000+c*10000+ch]=Parameter{};
 // Calibration.h uses array={1}; Init resets gain/MIP but NOT intercept.
 // Thus only element [0][0][0] starts at 1, and all other intercepts at 0.
 a.at(0).intercept=1;
 double sg[2]={},sb[2]={},sm[2]={};int ng[2]={},nm[2]={};
 {auto f=openFile(hl);auto*t=object<TTree>(*f,"InterCalib");int id;double s,b;bind(t,"CellID",&id);bind(t,"Slope",&s);bind(t,"Intercept",&b);
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id<0||id/100000>=30)continue;require(std::isfinite(s)&&std::isfinite(b),"Invalid HL parameter");if(s<.02||s>.05)continue;int g=group(id);sg[g]+=1/s;sb[g]+=-b/s;++ng[g];if(a.count(id)){a[id].gain=1/s;a[id].intercept=-b/s;}}
 }
 {auto f=openFile(mip);auto*t=object<TTree>(*f,"MIP_Fit");int id,ndf;double mpv,chi;bind(t,"CellID",&id);bind(t,"LandauMPV",&mpv);bind(t,"ChiSquare",&chi);bind(t,"NDF",&ndf);
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id<0||id/100000>=30)continue;int g=group(id);require(std::isfinite(chi)&&std::isfinite(mpv),"Invalid MIP fit");if(chi/ndf>(g?2.:1.7))continue;require(mpv!=0,"Zero MIP parameter");sm[g]+=mpv;++nm[g];if(a.count(id))a[id].mip=mpv;}
 }
 for(int g=0;g<2;++g)require(ng[g]&&nm[g],"No accepted constants in a SiPM group");
 for(auto&[id,p]:a){int g=group(id);p.hlFallback=p.gain==1||p.intercept==1;p.mipFallback=p.mip==1;if(p.gain==1)p.gain=sg[g]/ng[g];if(p.intercept==1)p.intercept=sb[g]/ng[g];if(p.mip==1)p.mip=sm[g]/nm[g];}
 {auto f=openFile(ped);auto*t=object<TTree>(*f,"ChnLevel");std::vector<int>*ids=nullptr;std::vector<double>*ph=nullptr,*pl=nullptr;bind(t,"CellID",&ids);bind(t,"PedHighMean",&ph);bind(t,"PedLowMean",&pl);
  // Calibration.h reads pedestal branches through vector<float>, so retain
  // the same rounding even when the parameter ROOT stores vector<double>.
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);require(ids->size()==ph->size()&&ids->size()==pl->size(),"Pedestal size mismatch");for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!a.count(id))continue;auto&p=a[id];p.ph=static_cast<float>(ph->at(j));p.pl=static_cast<float>(pl->at(j));p.ped=true;}}
  t->ResetBranchAddresses();delete ids;delete ph;delete pl;
 }
 return a;
}
inline std::map<std::string,std::string> options(int argc,char**argv,std::map<std::string,std::string> o) {
 for(int i=1;i<argc;i+=2){require(i+1<argc,"Option needs a value");std::string key=argv[i];require(key.rfind("--",0)==0&&o.count(key.substr(2)),"Unknown option "+key);o[key.substr(2)]=argv[i+1];}return o;
}
