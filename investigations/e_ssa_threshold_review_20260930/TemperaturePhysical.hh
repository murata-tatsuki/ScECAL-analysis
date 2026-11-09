#pragma once
// Diagnostic fixed-optical-yield model. Not a replacement for measured curves.
#include <array>
#include <cstdlib>
#include <limits>
#include <memory>
#ifndef PHYS_PDE_LINEAR
#define PHYS_PDE_LINEAR 0
#endif
#ifndef PHYS_TEMP_REAL
#define PHYS_TEMP_REAL 0
#endif
namespace TemperaturePhysical {
using Matrix=std::vector<std::vector<double>>;
inline std::vector<Matrix> bank;
inline Matrix current;
inline double alpha(int l){return l>3&&l<28?1.6/135:3.5/230;}
inline double u20(int l){return l>3&&l<28?4.5:4.;}
inline double p0(int l){return l>3&&l<28?.1:.15;}
inline double ratioU(int l,double t){return 1.-.06*(t-20.)/u20(l);}
inline double gainRatio(int l,double t){return 1.-alpha(l)*(t-20.);}
inline double pct(int l,double t){double r=ratioU(l,t);return p0(l)*r*r;}
inline double pdeRatio(int l,double t){return PHYS_PDE_LINEAR?ratioU(l,t):1.;}
inline double local(int l,double x,double y,const Matrix&m){
 const double sx[16]={39.5986,-106.197,89.9922,-22.1234,-103.886,-.508,-60.2488,59.1312,28.2702,-50.3428,100.99,-106.68,-106.299,105.918,-47.0916,53.2638};
 const double sy[16]={-34.3916,-34.0868,-11.1506,-13.1064,8.0772,7.366,34.5186,33.7312,54.0766,75.819,101.473,101.498,-101.752,-102.032,-70.6882,-70.3326};
 if(l<0||l>=30||m.size()!=30||m[l].size()!=16)throw std::runtime_error("Bad thermal matrix");
 double d[16],dm=std::numeric_limits<double>::infinity();for(int s=0;s<16;++s){double a=l%2?-sx[s]:-sy[s],b=l%2?-sy[s]:-sx[s];d[s]=std::hypot(x-a,y-b);if(d[s]==0)return m[l][s];dm=std::min(dm,d[s]);}
 long double sum=0,temp=0;for(int s=0;s<16;++s){double w=dm/d[s];sum+=w;temp+=w*m[l][s];}return temp/sum;
}
inline void initialize(Extract&ex,const char*led){
 std::unique_ptr<TFile>f(TFile::Open(led,"READ"));auto*t=dynamic_cast<TTree*>(f->Get("T_Event"));if(!t)throw std::runtime_error("Missing LED tree");int id=0;double temp=0;t->SetBranchAddress("cellIDs",&id);t->SetBranchAddress("temp",&temp);double measured[30][6][36]{};double sum[2]{},n[2]{};
 for(Long64_t j=0;j<t->GetEntries();++j){t->GetEntry(j);int l=id/1000000,p=id%1000000/10000,c=id%100;if(l>=30||p>=6||c>=36||!std::isfinite(temp))continue;measured[l][p][c]=temp;int g=l>3&&l<28;sum[g]+=temp;++n[g];}
 for(int l=0;l<30;++l)for(int p=0;p<6;++p)for(int c=0;c<36;++c){int g=l>3&&l<28;double mt=measured[l][p][c];if(mt==0)mt=sum[g]/n[g];double r=gainRatio(l,mt);if(!(r>0))throw std::runtime_error("Invalid LED thermal gain conversion");ex.SiPMGain[l][p][c]/=r;}t->ResetBranchAddresses();
 const char*path=std::getenv("TEMPERATURE_BANK");if(!path)throw std::runtime_error("TEMPERATURE_BANK required");std::unique_ptr<TFile>b(TFile::Open(path,"READ"));auto*bt=dynamic_cast<TTree*>(b->Get("TemperatureBank"));Matrix*m=nullptr;bt->SetBranchAddress("Temperature",&m);for(Long64_t j=0;j<bt->GetEntries();++j){bt->GetEntry(j);bank.push_back(*m);}bt->ResetBranchAddresses();delete m;
}
inline void set(size_t entry){size_t offset=std::stoul(std::getenv("TEMPERATURE_OFFSET")?std::getenv("TEMPERATURE_OFFSET"):"0");if(entry+offset>=bank.size())throw std::runtime_error("Temperature bank exhausted");current=PHYS_TEMP_REAL?bank[entry+offset]:Matrix(30,std::vector<double>(16,20.));}
}
