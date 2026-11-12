#pragma once
#include "../muon_mip_response/TrackGeometry.hh"
#include <TNamed.h>
#include <TH1D.h>
#include <Math/QuantFuncMathCore.h>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <tuple>
namespace de {
using namespace mmr;
using mmr_track::Settings;
constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
constexpr const char* schema="detection_efficiency_v1";
inline std::string slurp(const std::string&p){std::ifstream f(p);require(bool(f),"Cannot read "+p);return {std::istreambuf_iterator<char>(f),{}};}
inline std::ofstream table(const std::string&p,const std::string&h){require(!fs::exists(p),"Refusing overwrite: "+p);std::ofstream f(p);require(bool(f),"Cannot write "+p);f<<std::setprecision(12)<<h<<'\n';return f;}
inline void fresh(const std::string&p){require(!fs::exists(p),"Refusing overwrite: "+p);}
inline std::vector<int> integers(const std::string&s){std::vector<int> r;for(auto&v:split(s,',')){size_t n;int x=std::stoi(v,&n);require(n==v.size(),"Invalid integer list");r.push_back(x);}return r;}
struct Variant{std::string name;Settings cuts;};
inline std::vector<Variant> variants(const Config&c){Settings s(c);std::vector<Variant>v{{"nominal",s}};if(c.integer("stability",1)){auto a=s;a.longMargin=0;a.transMargin=0;v.push_back({"margin_loose",a});a=s;a.longMargin=std::max(2.,s.longMargin);a.transMargin=std::max(.75,s.transMargin);v.push_back({"margin_tight",a});a=s;a.minLayers=std::max(9,s.minLayers);a.rms=std::min(2.,s.rms);a.angle=std::min(.035,s.angle);v.push_back({"quality_tight",a});}return v;}
inline std::vector<std::vector<int>> groups(const Config&c){std::vector<std::vector<int>> v{{}};auto s=c.text("groups","0,1,2,3;10,11,12,13;26,27,28,29");if(s=="none")return v;for(auto&t:split(s,';')){auto a=integers(t);std::set<int> u(a.begin(),a.end());require(a.size()>=2&&a.size()==u.size(),"Groups need distinct >=2 layers");for(int l:a)require(l>=0&&l<30,"Invalid group layer");v.push_back(a);}return v;}
using Masks=std::map<std::pair<int,int>,std::string>;
inline Masks masks(const Config&c){Masks m;auto p=c.text("mask");if(p.empty())return m;std::ifstream f(p);require(bool(f),"Cannot open mask");for(std::string line;std::getline(f,line);){if(line.empty()||line[0]=='#')continue;auto v=split(line);require(v.size()==3,"Mask: sample TAB physical_CellID TAB reason");int s=std::stoi(v[0]),id=std::stoi(v[1]);require((s==0||s==1)&&connected(id)&&physical(id)==id,"Invalid mask channel");require(m.emplace(std::make_pair(s,id),v[2]).second,"Duplicate mask channel");}return m;}
struct Geometry{
 std::array<std::array<std::array<int,5>,42>,30> grid{};
 std::array<double,30> z{};std::map<int,mmr_track::StripRectangle> strips;
 Geometry(){for(auto&l:grid)for(auto&r:l)r.fill(-1);for(int id:channels()){auto p=geometry(id);int l=id/100000;int r=std::lround(20.5-(l%2?p[0]:p[1])/5.3),c=std::lround(2-(l%2?p[1]:p[0])/45.4);require(r>=0&&r<42&&c>=0&&c<5&&grid[l][r][c]<0,"Geometry grid mismatch");grid[l][r][c]=id;strips.emplace(id,mmr_track::StripRectangle{id,l,p});z[l]=p[2];}}
 int cell(int l,int r,int c)const{return r>=0&&r<42&&c>=0&&c<5?grid[l][r][c]:-1;}
};
#define CROSS_I(X) X(sample) X(run) X(source_run) X(file) X(time) X(cycle) X(trigger) X(variant) X(mode) X(layer) X(cellid) X(track_ok) X(chord) X(ambiguous) X(track_ambiguous) X(masked) X(eligible) X(multiplicity) X(memory_extra) X(repeated_records) X(neighbor_only) X(nx) X(ny)
#define CROSS_D(X) X(x) X(y) X(sx) X(sy) X(rms) X(longitudinal) X(transverse) X(length) X(sigma_x) X(sigma_y)
struct Crossing{
#define DECL_I(n) int n=0;
 CROSS_I(DECL_I)
#undef DECL_I
#define DECL_D(n) double n=NaN;
 CROSS_D(DECL_D)
#undef DECL_D
 Long64_t entry=0;
 void branches(TTree&t){
#define B(n) t.Branch(#n,&n);
 CROSS_I(B) CROSS_D(B) B(entry)
#undef B
 }
 void bindTo(TTree*t){
#define B(n) bind(t,#n,&n);
 CROSS_I(B) CROSS_D(B) B(entry)
#undef B
 }
};
struct Count{
 Long64_t attempts=0,track=0,geometry=0,n=0,k=0,unique=0,duplicate=0,memory_extra=0,repeated=0,neighbor=0,ambiguous=0,masked=0;
 void add(const Crossing&x){++attempts;track+=x.track_ok;geometry+=x.chord;ambiguous+=x.ambiguous||x.track_ambiguous;masked+=x.masked;if(!x.eligible)return;++n;k+=x.multiplicity>0;unique+=x.multiplicity==1;duplicate+=x.multiplicity>1;memory_extra+=x.memory_extra;repeated+=x.repeated_records;neighbor+=x.neighbor_only;}
};
struct Interval{double p=NaN,lo=NaN,hi=NaN;};
inline Interval wilson(double k,double n,double z){if(n<=0)return {};double p=k/n,d=1+z*z/n,c=(p+z*z/(2*n))/d,h=z*std::sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return {p,std::max(0.,c-h),std::min(1.,c+h)};}
// Newcombe difference of two independent Wilson intervals (method 10).
inline Interval difference(Interval a,Interval b){if(!std::isfinite(a.p)||!std::isfinite(b.p))return {};double d=a.p-b.p;return {d,std::max(-1.,d-std::hypot(a.p-a.lo,b.hi-b.p)),std::min(1.,d+std::hypot(a.hi-a.p,b.p-b.lo))};}
inline std::string state(Long64_t n,int minimum){return n==0?"unmeasured":n<minimum?"insufficient":"ok";}
inline uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}
inline uint64_t rank(const Crossing&x){return mix(uint64_t(x.run)^mix(uint32_t(x.time))^mix(uint32_t(x.cycle)+0x100000000ULL)^mix(uint32_t(x.trigger)+0x200000000ULL));}
// Integer binning uses only predicted geometry / non-target tracking, never ADC or hit outcome.
using Bin=std::array<int,8>;
inline Bin bin(const Crossing&x,const Config&c){return {x.variant,x.run,x.cellid,int(std::floor(x.longitudinal/c.number("long-bin",10))),int(std::floor(x.transverse/c.number("trans-bin",1))),int(std::floor(x.sx/c.number("slope-bin",.02))),int(std::floor(x.sy/c.number("slope-bin",.02))),int(std::floor(x.rms/c.number("rms-bin",1)))};}
inline std::string condition(int v,int matched,const Config&c){return (matched?"common_":"unadjusted_")+variants(c).at(v).name;}
} // namespace de
namespace de {
#define RESULT_I(X) X(sample) X(variant) X(matched) X(level) X(id) X(status)
#define RESULT_L(X) X(attempts) X(track) X(geometry) X(predicted) X(observed) X(unique) X(duplicate) X(memory_extra) X(repeated) X(neighbor) X(ambiguous) X(masked) X(events)
#define RESULT_D(X) X(eff) X(lo) X(hi) X(unique_eff) X(unique_lo) X(unique_hi) X(cluster_se)
struct Result{
#define D(n) int n=0;
 RESULT_I(D)
#undef D
#define D(n) Long64_t n=0;
 RESULT_L(D)
#undef D
#define D(n) double n=NaN;
 RESULT_D(D)
#undef D
 void branches(TTree&t){
#define B(n) t.Branch(#n,&n);
 RESULT_I(B) RESULT_L(B) RESULT_D(B)
#undef B
 }
 void bindTo(TTree*t){
#define B(n) bind(t,#n,&n);
 RESULT_I(B) RESULT_L(B) RESULT_D(B)
#undef B
 }
};
}
