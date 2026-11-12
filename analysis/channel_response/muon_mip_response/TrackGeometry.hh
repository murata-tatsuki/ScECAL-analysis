#pragma once
#include "MuonMIPCommon.hh"
#include <limits>
namespace mmr_track {
using mmr::require;
constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
struct Settings {
 Long64_t maxEvents;int minLayers;double angle,rms,longMargin,transMargin;bool allow;
 Settings(const mmr::Config&c):maxEvents(c.integer("max-events",0)),minLayers(c.integer("min-layers",7)),angle(c.number("max-angle",.05)),rms(c.number("max-track-rms",3)),longMargin(c.number("long-margin",1)),transMargin(c.number("trans-margin",.25)),allow(c.integer("allow-recovered",0)!=0){require(maxEvents>=0&&minLayers>=3&&minLayers<=14,"Invalid sample/layer limits");require(angle>0&&angle<.5&&rms>0&&longMargin>=0&&longMargin<22.5&&transMargin>=0&&transMargin<2.5,"Invalid track/margin limits");}
};
struct Point{int layer;double z,q;};
struct Fit{bool valid=false;double a=NaN,b=NaN,rms=NaN;int n=0;};
inline Fit ols(const std::vector<Point>&p){Fit f;if(p.size()<3)return f;double mz=0,mq=0;for(auto&v:p){mz+=v.z;mq+=v.q;}mz/=p.size();mq/=p.size();double zz=0,zq=0;for(auto&v:p){zz+=(v.z-mz)*(v.z-mz);zq+=(v.z-mz)*(v.q-mq);}if(zz<=0)return f;f.b=zq/zz;f.a=mq-f.b*mz;double ss=0;for(auto&v:p){double r=v.q-f.a-f.b*v.z;ss+=r*r;}f.rms=std::sqrt(ss/p.size());f.n=p.size();f.valid=true;return f;}
inline bool good(const Fit&f,const Settings&c){return f.valid&&f.n>=c.minLayers&&std::isfinite(f.rms)&&f.rms<=c.rms&&std::abs(std::atan(f.b))<=c.angle;}
inline Fit refine(const std::vector<Point>&points,double a,double b,const Settings&c){Fit f;for(int iteration=0;iteration<4;++iteration){std::map<int,std::pair<double,Point>>nearest;for(auto&p:points){double r=std::abs(p.q-a-b*p.z);if(r>std::max(5.3,2*c.rms))continue;auto it=nearest.find(p.layer);if(it==nearest.end()||r<it->second.first)nearest[p.layer]={r,p};}if(nearest.size()<size_t(c.minLayers))return {};std::vector<Point>selected;for(auto&v:nearest)selected.push_back(v.second.second);f=ols(selected);if(!f.valid||std::abs(std::atan(f.b))>c.angle)return {};a=f.a;b=f.b;}f.valid=good(f,c);return f;}
// Fine-coordinate duplicates collapse before fitting. One point per layer is
// chosen by geometric residual, never by ADC or the target-layer response.
inline Fit track(const std::array<std::vector<double>,30>&fine,const std::array<double,30>&z,int parity,int excluded,const Settings&c){std::vector<Point>points,initial;for(int l=parity;l<30;l+=2){if(l==excluded||fine[l].empty())continue;initial.push_back({l,z[l],fine[l][fine[l].size()/2]});for(double q:fine[l])points.push_back({l,z[l],q});}if(initial.size()<size_t(c.minLayers))return {};auto first=ols(initial);if(good(first,c)&&points.size()==initial.size()&&std::all_of(points.begin(),points.end(),[&](const Point&p){return std::abs(p.q-first.a-first.b*p.z)<=std::max(5.3,2*c.rms);}))return first;if(first.valid){auto fast=refine(points,first.a,first.b,c);if(good(fast,c))return fast;}Fit best;size_t n=points.size();for(size_t seed=0;seed<32;++seed){size_t i=(seed*37+seed/7)%n,j=(seed*83+n/2+1)%n;if(points[i].layer==points[j].layer)continue;double b=(points[j].q-points[i].q)/(points[j].z-points[i].z);if(std::abs(std::atan(b))>c.angle)continue;auto f=refine(points,points[i].q-b*points[i].z,b,c);if(good(f,c)&&(!best.valid||f.n>best.n||(f.n==best.n&&f.rms<best.rms)))best=f;}return best;}
struct StripRectangle{int id,layer;std::array<double,3>p;};
// Full z +/-1 mm segment stays inside the 45x5 mm active rectangle.
inline bool chord(const StripRectangle&s,double x,double y,double sx,double sy,const Settings&c,double&len,double&longAt,double&transAt){if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(sx)||!std::isfinite(sy))return false;longAt=s.layer%2?y-s.p[1]:x-s.p[0];transAt=s.layer%2?x-s.p[0]:y-s.p[1];double ls=s.layer%2?sy:sx,ts=s.layer%2?sx:sy;if(std::abs(longAt)+std::abs(ls)>22.5-c.longMargin||std::abs(transAt)+std::abs(ts)>2.5-c.transMargin)return false;len=2*std::sqrt(1+sx*sx+sy*sy);return true;}
}
