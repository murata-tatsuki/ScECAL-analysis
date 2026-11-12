#pragma once
#include "../muon_mip_response/CurrentCalibration.hh"
#include "../muon_mip_response/TrackGeometry.hh"
#include <TNamed.h>
#include <TParameter.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <tuple>
#include <climits>
namespace ts {
using namespace mmr;
constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
constexpr const char* schema="time_stability_1";
struct Moment {
 long long n=0; long double s=0,q=0;
 void add(double x){if(std::isfinite(x)){++n;s+=x;q+=x*x;}}
 double mean()const{return n?double(s/n):NaN;}
 double var()const{return n>1?double(std::max(0.L,(q-s*s/n)/(n-1))):NaN;}
 double sd()const{return std::sqrt(var());}
 double se()const{return n>1?sd()/std::sqrt(n):NaN;}
};
struct Base {
 int sample=0,particle=0,file=0,run=0,source_run=0,event=0,time=0,segment=0,selected=0;
 Long64_t entry=0;double seconds=NaN;
 void io(TTree&t,bool read){
#define B(x) if(read)bind(&t,#x,&x);else t.Branch(#x,&x)
 B(sample);B(particle);B(file);B(run);B(source_run);B(event);B(time);B(segment);B(selected);B(entry);B(seconds);
#undef B
 }
};
struct Event:Base {
 int paired=0,tail=0,cog_ok=0,recorded=0,strips=0,virtual_hits=-1,raw_hits=0,raw_strips=0,negative=0,zero=0,duplicate=0;
 double energy=NaN,signed_energy=NaN,temperature=NaN,cogx=NaN,cogy=NaN;
 double layer[30]={};
 void io(TTree&t,bool read){Base::io(t,read);
#define B(x) if(read)bind(&t,#x,&x);else t.Branch(#x,&x)
 B(paired);B(tail);B(cog_ok);B(recorded);B(strips);B(virtual_hits);B(raw_hits);B(raw_strips);B(negative);B(zero);B(duplicate);B(energy);B(signed_energy);B(temperature);B(cogx);B(cogy);
#undef B
 if(read)bind(&t,"layer",layer);else t.Branch("layer",layer,"layer[30]/D");
 }
};
struct Hit:Base {
 int cell=0,memory=-1,observed=0,unique=0,valid=0,temp_measured=0;
 double before=NaN,after=NaN,temperature=NaN,sensor=NaN,factor=NaN,low=NaN,high=NaN;
 double x=NaN,y=NaN,sx=NaN,sy=NaN,length=NaN,longitudinal=NaN,transverse=NaN;
 void io(TTree&t,bool read){Base::io(t,read);
#define B(x) if(read)bind(&t,#x,&x);else t.Branch(#x,&x)
 B(cell);B(memory);B(observed);B(unique);B(valid);B(temp_measured);B(before);B(after);B(temperature);B(sensor);B(factor);B(low);B(high);B(x);B(y);B(sx);B(sy);B(length);B(longitudinal);B(transverse);
#undef B
 }
};
struct Input {int sample,particle,run;std::string raw,calib,truth,evidence;};
inline std::vector<Input> inputs(const std::string&p){std::ifstream f(p);require(bool(f),"Cannot open manifest");std::vector<Input>v;std::set<std::string>seen;
 for(std::string line;std::getline(f,line);){if(line.empty()||line[0]=='#')continue;auto a=split(line);if(a[0]=="sample")continue;require(a.size()==7,"Manifest: sample particle run decode calib truth evidence (7 TAB columns)");Input i{std::stoi(a[0]),a[1]=="muon"?0:a[1]=="electron"?1:a[1]=="pedestal"?2:-1,std::stoi(a[2]),a[3],a[4],a[5],a[6]};require((i.sample==0||i.sample==1)&&i.particle>=0,"Bad sample/particle");require(seen.insert(fs::canonical(i.raw).string()).second,"Repeated raw input");for(auto&p:{i.raw,i.calib,i.truth,i.evidence})if(p!="-")require(fs::is_regular_file(p),"Missing input "+p);require(i.particle==2||i.calib!="-","Calibration required");require(!i.sample||i.truth!="-","MC truth required (muon must be separate production)");v.push_back(i);}require(!v.empty(),"Empty manifest");return v;}
inline std::set<int> mask(const Config&c){std::set<int>s;if(c.text("mask","-")=="-")return s;std::ifstream f(c.text("mask"));require(bool(f),"Missing mask");for(std::string line;std::getline(f,line);)if(!line.empty()&&line[0]!='#'){int id=std::stoi(line);require(connected(id)&&physical(id)==id,"Mask must contain physical connected CellIDs");s.insert(id);}return s;}
inline std::vector<std::pair<int,int>> groups(int id){return {{0,id},{1,id/100000*6+id%100000/10000},{2,id/100000}};}
inline std::string groupName(int g){return g==0?"channel":g==1?"chip":g==2?"layer":"event";}
inline std::string particleName(int p){return p==0?"muon":p==1?"electron":"pedestal";}
inline std::string sampleName(int s){return s?"mc":"data";}
inline std::string condition(int s,int p,const Config&c){return p==1?(c.integer(s?"mc-tail":"data-tail")?"electron_exclude":"electron_keep"):p==0?"muon":"pedestal";}
inline fs::path subdir(int s,int p,const Config&c){fs::path pth=sampleName(s);pth/=condition(s,p,c);if(p==1)pth/=c.text("cog","200")+"mm";pth/=(p==0?c.text("muon-energy","100"):c.text("energy","100"))+"GeV";return pth;}
inline std::string axisName(int axis){return axis==0?"event_order":axis==1?"elapsed_seconds":axis==2?"run":axis==3?"temperature_C":axis==4?"timestamp_ticks":"pooled";}
struct Axis {int axis,run,file,segment,bin;double x;};
inline std::vector<Axis> axes(const Base&b,double temperature,const Config&c){
 std::vector<Axis> a={{5,-1,-1,0,0,0},{0,b.run,b.file,b.segment,int(b.entry/c.integer("event-bin",250)),double(b.entry)},{2,b.run,-1,0,b.run,double(b.run)}};
 if(std::isfinite(b.seconds))a.push_back({1,b.run,b.file,b.segment,int(std::floor(b.seconds/c.number("time-bin",60))),b.seconds});
 if(b.sample==0){a.push_back({4,b.run,b.file,b.segment,int(std::floor(double(b.time)/c.number("tick-bin",60))),double(b.time)});if(std::isfinite(temperature))a.push_back({3,b.run,-1,0,int(std::floor(temperature/c.number("temperature-bin",.25))),temperature});}
 return a;
}
}
