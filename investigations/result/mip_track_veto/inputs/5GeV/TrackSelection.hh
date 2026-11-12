#pragma once
#include "Common.hh"

// Operational through-going, sparse, straight-track tag, NOT particle identification.
// Uses saved layer moments: only the short strip coordinate is fitted.
struct TrackConfig {
 bool enabled=false,v2=false; double beam=0,lowSignalMax=.1;
 double maxEnergy=.03,maxHitsPerLayer=2,maxRms=5.3;
 int maxHits=60,minLayers=15,minEndLayers=2;
 static TrackConfig parse(const std::map<std::string,std::string>& o) {
  TrackConfig c;require(o.at("track-cut")=="on"||o.at("track-cut")=="off"||o.at("track-cut")=="v2","track-cut must be on, off or v2");c.enabled=o.at("track-cut")!="off";c.v2=o.at("track-cut")=="v2";c.beam=beamEnergy(o.at("energy"));
  auto number=[&](const char* key){size_t n=0;double x=std::stod(o.at(key),&n);require(n==o.at(key).size()&&std::isfinite(x),std::string("Invalid ")+key);return x;};
  auto integer=[&](const char* key){double x=number(key);require(x==std::floor(x)&&x>=0&&x<=100000,std::string("Invalid integer ")+key);return int(x);};
  c.lowSignalMax=number("low-signal-max-energy-gev");require(c.lowSignalMax>0,"Invalid low-signal threshold");c.maxEnergy=number("track-max-energy-gev");c.maxHitsPerLayer=number("track-max-hits-per-layer");c.maxRms=number("track-max-rms-mm");
  c.maxHits=integer("track-max-hits");c.minLayers=integer("track-min-layers");c.minEndLayers=integer("track-min-end-layers");
  require(c.maxEnergy>0&&c.maxHits>0&&c.maxHitsPerLayer>=1&&c.maxRms>0&&c.minLayers>=8&&c.minLayers<=30&&c.minEndLayers>=1&&c.minEndLayers<=5,"Invalid track threshold range");return c;
 }
 std::string caption()const {return v2?"; v2: low-signal + track veto":enabled?"; track veto":"";}
 std::string definition()const {
  if(v2){TrackConfig old=*this;old.v2=false;std::ostringstream s;s<<std::setprecision(17)<<"track_cut=v2; reject = old_track_tag OR ((beam_energy_GeV >= 10) AND (0 < total_calibrated_energy_GeV <= "<<lowSignalMax<<")); beam_energy_GeV="<<beam<<"; low-signal veto does NOT require hits, straightness or penetration; diagnostic ENERGY selection, NOT muon/pion PID or leakage identification; old_track_tag defined by: "<<old.definition();return s.str();}
  std::ostringstream s;s<<std::setprecision(17)<<"track_cut="<<(enabled?"on":"off")<<"; tag = (0 < energy_GeV <= "<<maxEnergy<<") && (hits <= "<<maxHits<<") && (active_layers >= "<<minLayers<<") && (active_front5 >= "<<minEndLayers<<") && (active_back5 >= "<<minEndLayers<<") && (hits/active_layers <= "<<maxHitsPerLayer<<") && (nx >= 4 && ny >= 4) && (rms_x_mm <= "<<maxRms<<" && rms_y_mm <= "<<maxRms<<"); active means layer_energy>0; hits includes accepted zero-energy and repeated memory hits; unweighted layer fits x(z) odd layers, y(z) even layers; rms includes within-layer short-coordinate variance; tag is NOT muon/pion identification; on rejects tag after low-LG-tail and positive-energy requirements, before CoG";return s.str();
 }
};
struct TrackTag {
 int active=0,front=0,back=0,nx=0,ny=0;
 double occupancy=NaN,rmsX=NaN,rmsY=NaN;
 bool low=false,sparse=false,spanning=false,straight=false,oldTrack=false,lowSignal=false,tag=false;
};
inline std::array<double,30> trackLayerZ(){std::array<double,30> z{};for(int l=0;l<30;++l)z[l]=EBUdecode(l,0,0)[2];return z;}
inline TrackTag trackTag(const Event&e,const TrackConfig&c,const std::array<double,30>&z) {
 TrackTag t;
 for(int l=0;l<30;++l)if(e.layer_energy[l]>0){++t.active;if(l<5)++t.front;if(l>=25)++t.back;}
 t.occupancy=t.active?double(e.hits)/t.active:NaN;
 auto fit=[&](bool x,int&n){
  double sz=0,sq=0;for(int l=x?1:0;l<30;l+=2)if(e.layer_energy[l]>0){++n;sz+=z[l];sq+=x?e.cog_x[l]:e.cog_y[l];}
  if(n<4)return NaN;
  const double mz=sz/n,mq=sq/n;double zz=0,zq=0;
  for(int l=x?1:0;l<30;l+=2)if(e.layer_energy[l]>0){double dz=z[l]-mz;zz+=dz*dz;zq+=dz*((x?e.cog_x[l]:e.cog_y[l])-mq);}
  if(zz<=0)return NaN;
  const double slope=zq/zz;double ss=0;
  for(int l=x?1:0;l<30;l+=2)if(e.layer_energy[l]>0){double q=x?e.cog_x[l]:e.cog_y[l],w=x?e.width_x[l]:e.width_y[l];if(!std::isfinite(q)||!std::isfinite(w))return NaN;double residual=q-mq-slope*(z[l]-mz);ss+=residual*residual+w*w;}
  return std::sqrt(ss/n); // geometric RMS, not a fit uncertainty or chi-square
 };
 t.rmsX=fit(true,t.nx);t.rmsY=fit(false,t.ny);
 t.low=e.energy>0&&e.energy<=c.maxEnergy;
 t.sparse=e.hits<=c.maxHits&&t.occupancy<=c.maxHitsPerLayer;
 t.spanning=t.active>=c.minLayers&&t.front>=c.minEndLayers&&t.back>=c.minEndLayers;
 t.straight=std::isfinite(t.rmsX)&&std::isfinite(t.rmsY)&&t.rmsX<=c.maxRms&&t.rmsY<=c.maxRms;
 t.oldTrack=t.low&&t.sparse&&t.spanning&&t.straight;t.lowSignal=c.beam>=10&&e.energy>0&&e.energy<=c.lowSignalMax;t.tag=t.oldTrack||(c.v2&&t.lowSignal);return t;
}
