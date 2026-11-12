#pragma once
#include "../resolution/HLTailSelection.hh"
#include "../impact_studies/StudyCommon.hh"
#include <TROOT.h>
#include "../calibration_parameters/EBUdecode.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TBox.h>
#include <TLatex.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <random>
#include <regex>
#include <set>
constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
inline std::set<int> readMask(const std::string& path) {
 // Same extraction as hg_lg_tail; normalize memory-cell digits before lookup.
 std::ifstream in(path);require(bool(in),"Cannot read mask "+path);
 std::string s((std::istreambuf_iterator<char>(in)),{});
 auto a=s.find("const unordered_set<int> bad_channels = {");require(a!=std::string::npos,"Missing mask");
 a=s.find('{',a);auto b=s.find("};",a);require(b!=std::string::npos,"Invalid mask");s=s.substr(a,b-a);
 std::set<int> out;std::regex re("[0-9]+");for(std::sregex_iterator i(s.begin(),s.end(),re),end;i!=end;++i)out.insert(physical(std::stoi(i->str())));return out;
}
struct StripGeometry { double x,y,z,hx,hy,bx,by;bool edgeX,edgeY,masked; };
inline std::map<int,StripGeometry> geometry(double width,const std::set<int>&mask) {
 require(width>0&&width<110,"edge-mm must be in (0,110)");std::map<int,StripGeometry> result;
 for(int l=0;l<30;++l)for(int c=0;c<6;++c)for(int h=0;h<(c==5?30:36);++h){
  int id=l*100000+c*10000+h;auto*p=EBUdecode(l,c,h);StripGeometry s{p[0],p[1],p[2],l%2?2.5:22.5,l%2?22.5:2.5,l%2?111.15:113.3,l%2?113.3:111.15,false,false,bool(mask.count(id))};
  s.edgeX=std::abs(s.x)+s.hx>s.bx-width;s.edgeY=std::abs(s.y)+s.hy>s.by-width;result[id]=s;
 }return result;
}
struct Event {
 int file_index=0,run=0,time=0,event=0,hits=0,negative_hits=0,masked_hits=0,disconnected_hits=0,max_layer=-1;
 Long64_t entry=0;bool tail_rich=false,positive=false,cog_valid=false,pass20=false;
 double energy=0,negative_energy=0,back3=0,back5=0,edge=0,edge_x=0,edge_y=0,depth=NaN,depth_layer=NaN;
 double indicator_fraction[5]{},layer_energy[30]{},layer_fraction[30]{},cog_x[30]{},cog_y[30]{},width_x[30]{},width_y[30]{};int layer_hits[30]{};
 void branches(TTree&t,bool write) {
#define FIELD(x,type) if(write)t.Branch(#x,&x,#x "/" type);else bind(&t,#x,&x)
 FIELD(file_index,"I");FIELD(run,"I");FIELD(time,"I");FIELD(event,"I");FIELD(hits,"I");FIELD(negative_hits,"I");FIELD(masked_hits,"I");FIELD(disconnected_hits,"I");FIELD(max_layer,"I");FIELD(entry,"L");FIELD(tail_rich,"O");FIELD(positive,"O");FIELD(cog_valid,"O");FIELD(pass20,"O");
 FIELD(energy,"D");FIELD(negative_energy,"D");FIELD(back3,"D");FIELD(back5,"D");FIELD(edge,"D");FIELD(edge_x,"D");FIELD(edge_y,"D");FIELD(depth,"D");FIELD(depth_layer,"D");
#undef FIELD
#define ARRAY(x,type) if(write)t.Branch(#x,x,#x "[30]/" type);else bind(&t,#x,x)
 if(write)t.Branch("indicator_fraction",indicator_fraction,"indicator_fraction[5]/D");else bind(&t,"indicator_fraction",indicator_fraction);
 ARRAY(layer_energy,"D");ARRAY(layer_fraction,"D");ARRAY(layer_hits,"I");ARRAY(cog_x,"D");ARRAY(cog_y,"D");ARRAY(width_x,"D");ARRAY(width_y,"D");
#undef ARRAY
 }
};
inline const std::array<std::string,5> indicators={"back3","back5","edge","edge_x","edge_y"};
inline double absolute(const Event&e,int k){return k==0?e.back3:k==1?e.back5:k==2?e.edge:k==3?e.edge_x:e.edge_y;}
inline double fraction(const Event&e,int k){return e.energy>0?absolute(e,k)/e.energy:NaN;}
inline std::string energyLabel(double en){std::ostringstream text;text<<en;return text.str();}
inline double beamEnergy(const std::string& text){
 const std::set<std::string> supported={"0.5","1","2","3","4","5","10","20","30","40","50","60","70","80","100","120","150","200","250"};
 require(supported.count(text),"Unsupported beam energy: "+text);return std::stod(text);
}
inline fs::path samplePath(const fs::path&base,int s,int cog,double en){return base/(s?"mc":"data")/"exclude"/(std::to_string(cog)+"mm")/(energyLabel(en)+"GeV");}
inline double quantile(const std::vector<double>&v,double p){if(v.empty())return NaN;double k=p*(v.size()-1);size_t i=k;return v[i]+(k-i)*(v[std::min(i+1,v.size()-1)]-v[i]);}
struct Metrics {size_t n=0;double mean=NaN,sem=NaN,rms=NaN,median=NaN,w68=NaN,gaus=NaN,gaus_err=NaN,fit_lo=NaN,fit_hi=NaN;int fit_status=-1;};
inline Metrics metrics(std::vector<double>v){
 Metrics m;m.n=v.size();if(v.empty())return m;Stats s;for(double x:v)s.add(x);std::sort(v.begin(),v.end());m.mean=s.mean;m.sem=s.sem();m.rms=s.mean>0?s.rms()/s.mean:NaN;m.median=quantile(v,.5);m.w68=m.median>0?(quantile(v,.84)-quantile(v,.16))/(2*m.median):NaN;
 if(v.size()<100||s.rms()<=0)return m;
 // Fixed prescription, never an event cut: binned Gaussian over q16..q84.
 m.fit_lo=quantile(v,.16);m.fit_hi=quantile(v,.84);if(m.fit_hi<=m.fit_lo)return m;
 TH1D h("fit_work","",100,0,std::max(v.back()*1.001,1e-6));h.SetDirectory(nullptr);for(double x:v)h.Fill(x);
 TF1 f("core_gaussian","gaus",m.fit_lo,m.fit_hi);f.SetParameters(h.GetMaximum(),m.median,(m.fit_hi-m.fit_lo)/2);
 auto r=h.Fit(&f,"QSNR");m.fit_status=int(r);
 if(m.fit_status==0&&r.Get()&&r->CovMatrixStatus()>=2&&f.GetParameter(1)>0&&f.GetParameter(2)>0){double mu=f.GetParameter(1),sig=f.GetParameter(2);m.gaus=sig/mu;double var=r->CovMatrix(2,2)/(mu*mu)+sig*sig*r->CovMatrix(1,1)/std::pow(mu,4)-2*sig*r->CovMatrix(1,2)/std::pow(mu,3);m.gaus_err=std::sqrt(std::max(0.,var));}else if(m.fit_status==0)m.fit_status=100;
 return m;
}
inline void metricHeader(std::ostream&o){o<<"n\tmean_GeV\tmean_sem\tgaussian_sigma_mean\tgaussian_error\tcentral68_halfwidth_median\trms_mean\tmedian\tfit_status\tfit_lo\tfit_hi\tstatistics";}
inline void metricRow(std::ostream&o,const Metrics&m){o<<m.n<<'\t'<<m.mean<<'\t'<<m.sem<<'\t'<<m.gaus<<'\t'<<m.gaus_err<<'\t'<<m.w68<<'\t'<<m.rms<<'\t'<<m.median<<'\t'<<m.fit_status<<'\t'<<m.fit_lo<<'\t'<<m.fit_hi<<'\t'<<(m.n<100?"LOW_STAT":"OK");}
