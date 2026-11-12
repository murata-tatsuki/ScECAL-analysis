#pragma once
#include "Common.hh"

namespace selected_gain {
struct Spectrum {
 std::array<const Json*,4> methods{};
 const Json& base()const{for(auto* p:methods)if(p)return *p;throw std::runtime_error("Empty spectrum");}
 const Json* primary(const Settings& cfg)const{return methods[cfg.primary=="multi_gaussian"?0:1];}
};
inline bool candidate(const Json* r){return r&&str(*r,"status")=="candidate"&&std::isfinite(num(*r,"gain"))&&num(*r,"gain")>0;}
inline Json evidence(const Spectrum& s,const Settings& cfg) {
 const auto& b=s.base();auto cond=conditionFields(b);const auto* primary=s.primary(cfg);
 Json e={{"condition",cond["condition"]},{"primary_method",cfg.primary},{"qualified",0},{"used",0}};
 for(auto k:{"sample","source_raw","source_sample"})e[k]=str(b,k);
 for(auto k:{"cellid","sid","dac","temperature_C","temperature_RMS_C","entries"})e[k]=num(b,k);
 std::vector<std::string> why;
 if(cfg.dead.count(physicalCell(b)))why.push_back("known_dead_channel");
 if(num(b,"memory_cell")!=0||num(b,"hit_tag")!=1||num(b,"gain_tag")!=1)why.push_back("unsupported_tag_or_memory_cell");
 if(!std::isfinite(num(b,"temperature_C"))||num(b,"temperature_bin")<0||!(num(b,"valid_temperature_entries")>0))why.push_back("temperature_missing");
 if(!std::isfinite(num(b,"dac"))||num(b,"dac")<0||num(b,"dac")!=std::floor(num(b,"dac")))why.push_back("invalid_led_dac");
 if(!candidate(primary))why.push_back("primary_not_candidate");
 if(!candidate(s.methods[2]))why.push_back("fft_not_candidate");
 double largest=0;
 for(int i=0;i<4;++i){const auto* r=s.methods[i];e[Methods[i]+"_gain"]=r?num(*r,"gain"):NaN;e[Methods[i]+"_status"]=r?str(*r,"status"):"missing";
  if(candidate(primary)&&candidate(r)&&r!=primary){largest=std::max(largest,std::abs(num(*r,"gain")-num(*primary,"gain")));if(relative(num(*r,"gain"),num(*primary,"gain"))>cfg.agreement)why.push_back("method_conflict_"+Methods[i]);}
 }
 if(primary){e["primary_gain"]=num(*primary,"gain");e["fit_error"]=num(*primary,"fit_error");e["stability_ADC"]=num(*primary,"stability_ADC");
  if(!std::isfinite(num(*primary,"fit_error"))||num(*primary,"fit_error")<=0)why.push_back("invalid_fit_error");
  if(!std::isfinite(num(*primary,"stability_ADC"))||num(*primary,"stability_ADC")<0)why.push_back("invalid_stability");
  if(!std::isfinite(num(*primary,"deviance"))||num(*primary,"deviance")<0||!std::isfinite(num(*primary,"ndf"))||!(num(*primary,"ndf")>0))why.push_back("invalid_fit_quality");
 }
 e["method_spread_ADC"]=largest;e["qualified"]=why.empty()?1:0;e["status"]=why.empty()?"qualified":"rejected";e["reasons"]=reasons(why);return e;
}
struct GroupResult {Json decision;std::vector<Json> evidence,beam;};
inline GroupResult selectGroup(const std::vector<const Spectrum*>& spectra,const std::vector<const Spectrum*>& beam,const Settings& cfg) {
 require(!spectra.empty(),"Empty selection group");GroupResult out;auto& d=out.decision;const auto& b=spectra.front()->base();d=conditionFields(b);int cell=physicalCell(b);
 d.update({{"cellid",cell},{"accepted",0},{"status","unresolved"},{"method",cfg.primary},{"n_spectra",spectra.size()},{"gain_error_kind","conservative_quality_scale_not_statistical"},{"quality_statistic","representative_fit_poisson_deviance"},{"temperature_error_kind","observed_spread_not_standard_error"}});
 std::map<int,std::vector<double>> byDac;std::vector<double> gains,errors,temperatures;double stability=0,methodSpread=0,tempRms=0;std::vector<const Json*> qualified;
 for(auto* s:spectra){require(str(conditionFields(s->base()),"condition")==str(d,"condition"),"Mixed selection conditions");out.evidence.push_back(evidence(*s,cfg));const auto& e=out.evidence.back();if(num(e,"qualified")!=1)continue;
  const auto* r=s->primary(cfg);double g=num(*r,"gain");qualified.push_back(r);gains.push_back(g);errors.push_back(num(*r,"fit_error"));temperatures.push_back(num(*r,"temperature_C"));byDac[int(num(s->base(),"dac"))].push_back(g);
  stability=std::max(stability,num(*r,"stability_ADC"));methodSpread=std::max(methodSpread,num(e,"method_spread_ADC"));double rms=num(*r,"temperature_RMS_C");if(std::isfinite(rms))tempRms=std::max(tempRms,rms);
 }
 d["n_qualified"]=qualified.size();d["n_dacs"]=byDac.size();std::vector<std::string> why;
 if(cfg.dead.count(cell))why.push_back("known_dead_channel");
 if(byDac.size()<size_t(cfg.minDacs))why.push_back("too_few_qualified_dacs");
 if(!qualified.empty()) {
  std::vector<double> dacGains;for(auto& [dac,v]:byDac)dacGains.push_back(median(v));double g=median(dacGains),temp=median(temperatures),spread=0,tdev=0;
  for(double x:gains)spread=std::max(spread,std::abs(x-g));
  for(double x:temperatures)tdev=std::max(tdev,std::abs(x-temp));
  double tspan=*std::max_element(temperatures.begin(),temperatures.end())-*std::min_element(temperatures.begin(),temperatures.end());
  if(spread/g>cfg.spread)why.push_back("dac_or_repeat_inconsistent");
  if(tspan>cfg.temperatureSpan)why.push_back("temperature_span_too_wide");
  const Json* rep=nullptr;for(auto* r:qualified)if(!rep||std::abs(num(*r,"gain")-g)<std::abs(num(*rep,"gain")-g)||(std::abs(num(*r,"gain")-g)==std::abs(num(*rep,"gain")-g)&&num(*r,"fit_error")<num(*rep,"fit_error")))rep=r;
  d.update({{"gain",g},{"gain_error",std::max({median(errors),spread,stability,methodSpread})},{"temperature_C",temp},{"temperature_error_C",std::max(tdev,tempRms)},{"temperature_span_C",tspan},{"fit_error_scale",median(errors)},{"dac_spread_ADC",spread},{"stability_scale_ADC",stability},{"method_spread_ADC",methodSpread},
   {"representative_sample",str(*rep,"sample")},{"representative_sid",num(*rep,"sid")},{"representative_source_raw",str(*rep,"source_raw")},{"representative_gain",num(*rep,"gain")},{"representative_deviance",num(*rep,"deviance")},{"representative_ndf",num(*rep,"ndf")}});
  // Beam spectra are diagnostic corroboration only: the run settings are not
  // necessarily the LED settings, even at the same temperature.
  for(auto* s:beam){const auto& bb=s->base();const auto* r=s->primary(cfg);if(!candidate(r)||!candidate(s->methods[2])||relative(num(*r,"gain"),num(*s->methods[2],"gain"))>cfg.agreement)continue;
   if(str(bb,"beam")!=str(b,"beam")||num(bb,"memory_cell")!=num(b,"memory_cell")||num(bb,"hit_tag")!=num(b,"hit_tag")||num(bb,"gain_tag")!=num(b,"gain_tag"))continue;
   double dt=num(bb,"temperature_C")-temp;if(!std::isfinite(dt)||std::abs(dt)>cfg.beamTemperature)continue;
   auto sa=str(b,"settings"),sb=str(bb,"settings");if(sa!="unconfirmed"&&sb!="unconfirmed"&&sa!=sb)continue;
   double diff=relative(num(*r,"gain"),g);out.beam.push_back({{"condition",d["condition"]},{"cellid",cell},{"sample",str(bb,"sample")},{"kind",str(bb,"kind")},{"sid",num(bb,"sid")},{"settings",sb},{"gain",num(*r,"gain")},{"temperature_C",num(bb,"temperature_C")},{"delta_temperature_C",dt},{"relative_difference",diff},{"agrees",diff<=cfg.agreement?1:0},{"status",diff<=cfg.agreement?"supports":"disagrees"},{"condition_status",sa=="unconfirmed"||sb=="unconfirmed"?"settings_unconfirmed_different_run":"declared_settings_match_different_run"}});
  }
 }
 int agree=0;for(auto& r:out.beam)agree+=num(r,"agrees")==1;d["n_beam_checks"]=out.beam.size();d["n_beam_agree"]=agree;d["beam_check_status"]=out.beam.empty()?"not_available":agree==int(out.beam.size())?"supports":agree==0?"disagrees":"mixed";
 d["accepted"]=why.empty()?1:0;d["status"]=why.empty()?"selected":"unresolved";d["reasons"]=reasons(why);
 for(auto& e:out.evidence)e["used"]=why.empty()&&num(e,"qualified")==1?1:0;
 return out;
}
inline std::vector<GroupResult> selectChannel(const Json& rows,const Settings& cfg) {
 std::map<std::pair<std::string,long long>,Spectrum> spectra;
 for(auto& r:rows){require(num(r,"bootstrap_requested")==0,"Expected bootstrap-free comparison");auto method=std::find(Methods.begin(),Methods.end(),str(r,"method"));require(method!=Methods.end(),"Unknown method in comparison");double sid=num(r,"sid");require(std::isfinite(sid)&&sid>=0&&sid==std::floor(sid),"Invalid spectrum id");auto& s=spectra[{str(r,"sample"),static_cast<long long>(sid)}];size_t index=method-Methods.begin();require(!s.methods[index],"Duplicate spectrum/method result");
  if(std::any_of(s.methods.begin(),s.methods.end(),[](const Json* p){return p!=nullptr;})){const auto& b=s.base();for(auto key:{"cellid","temperature_bin","hit_tag","gain_tag","memory_cell"})require(num(b,key)==num(r,key),"Method metadata mismatch");for(auto key:{"kind","beam","period","settings","dac"})require(str(b,key)==str(r,key),"Method condition mismatch");}
  s.methods[index]=&r;
 }
 std::map<std::string,std::vector<const Spectrum*>> groups;std::vector<const Spectrum*> beam;
 for(auto& [key,s]:spectra){auto kind=str(s.base(),"kind");require(kind=="led"||kind=="electron"||kind=="muon","Unknown data kind");if(kind=="led")groups[str(conditionFields(s.base()),"condition")].push_back(&s);else beam.push_back(&s);}
 std::vector<GroupResult> result;for(auto& [key,group]:groups)result.push_back(selectGroup(group,beam,cfg));return result;
}
}
