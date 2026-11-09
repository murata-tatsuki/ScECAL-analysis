#pragma once
#include "../sipm_gain/compare/CompareIO.hh"
#include <TBranch.h>
#include <TChain.h>
#include <TH1D.h>
#include <TStyle.h>
#include <array>
#include <cmath>
#include <numeric>

namespace selected_gain {
using namespace gain;
using namespace gain::comparison;
inline const std::string DefaultReference = "/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/20230504_1425_Calib_LED_all.root";
inline const Schema DecisionSchema{
 {"condition", "beam", "period", "settings", "status", "reasons", "method", "representative_sample", "representative_source_raw", "gain_error_kind", "quality_statistic", "temperature_error_kind", "beam_check_status"},
 {"cellid", "temperature_bin", "hit_tag", "gain_tag", "memory_cell", "accepted", "gain", "gain_error", "temperature_C", "temperature_error_C", "temperature_span_C", "n_spectra", "n_qualified", "n_dacs", "fit_error_scale", "dac_spread_ADC", "stability_scale_ADC", "method_spread_ADC", "representative_sid", "representative_gain", "representative_deviance", "representative_ndf", "n_beam_checks", "n_beam_agree"}};
inline const Schema EvidenceSchema{
 {"condition", "sample", "source_raw", "source_sample", "status", "reasons", "primary_method", "multi_gaussian_status", "envelope_comb_status", "fft_status", "autocorrelation_status"},
 {"cellid", "sid", "dac", "temperature_C", "temperature_RMS_C", "entries", "qualified", "used", "primary_gain", "fit_error", "stability_ADC", "method_spread_ADC", "multi_gaussian_gain", "envelope_comb_gain", "fft_gain", "autocorrelation_gain"}};
inline const Schema BeamSchema{
 {"condition", "sample", "kind", "settings", "condition_status", "status"},
 {"cellid", "sid", "gain", "temperature_C", "delta_temperature_C", "relative_difference", "agrees"}};
inline const std::vector<std::string> OutputTrees = {"Decisions", "Evidence", "BeamChecks"};
inline const std::array<Schema,3> OutputSchemas = {DecisionSchema, EvidenceSchema, BeamSchema};

inline std::string digest(const std::string& s) {
 unsigned char bytes[EVP_MAX_MD_SIZE]; unsigned n=0;
 require(EVP_Digest(s.data(),s.size(),bytes,&n,EVP_sha256(),nullptr)==1,"Cannot hash selection metadata");
 std::ostringstream out; out<<std::hex<<std::setfill('0'); for(unsigned i=0;i<n;++i)out<<std::setw(2)<<unsigned(bytes[i]); return out.str();
}
inline double median(std::vector<double> v) {
 require(!v.empty(),"Median of empty selection"); for(double x:v)require(std::isfinite(x),"Nonfinite median input");
 std::sort(v.begin(),v.end()); size_t n=v.size(); return n%2?v[n/2]:(v[n/2-1]+v[n/2])/2;
}
inline double relative(double a,double b) { return a>0&&b>0?2*std::abs(a-b)/(a+b):INFINITY; }
inline std::string reasons(const std::vector<std::string>& v) { std::string s; for(auto& r:v){if(!s.empty())s+=';';s+=r;}return s; }
inline std::string integerText(double n) { require(std::isfinite(n)&&n==std::floor(n),"Invalid integer metadata");return std::to_string(static_cast<long long>(n)); }
inline Json conditionFields(const Json& row) {
 Json c; for(auto k:{"beam","period","settings"})c[k]=str(row,k);
 for(auto k:{"temperature_bin","hit_tag","gain_tag","memory_cell"})c[k]=num(row,k);
 require(str(c,"beam")=="ps"||str(c,"beam")=="sps","Unknown beam"); token(str(c,"period"));
 c["condition"]=str(c,"beam")+"_"+str(c,"period")+"_t"+integerText(num(c,"temperature_bin"))+"_h"+integerText(num(c,"hit_tag"))+"_g"+integerText(num(c,"gain_tag"))+"_m"+integerText(num(c,"memory_cell"))+"_s"+digest(str(c,"settings")).substr(0,16);
 return c;
}
inline Json temperatureInterval(double bin,double width) {
 require(std::isfinite(width)&&width>0&&std::isfinite(bin)&&bin==std::floor(bin)&&bin>=-1,"Invalid temperature interval");
 return {{"valid",bin>=0},{"width_C",width},{"lower_C",bin>=0?Json(bin*width):Json()},{"upper_C",bin>=0?Json((bin+1)*width):Json()}};
}
inline std::string temperatureDirectory(const Json& interval) {
 if(!interval.value("valid",false))return "temperature_unavailable";
 const double width=num(interval,"width_C");int precision=1;
 while(precision<8&&std::abs(width*std::pow(10.,precision)-std::round(width*std::pow(10.,precision)))>1e-8)++precision;
 std::ostringstream s;s<<std::fixed<<std::setprecision(precision)<<"temperature_"<<num(interval,"lower_C")<<"-"<<num(interval,"upper_C")<<"C";return s.str();
}
inline fs::path conditionPath(const Json& c,const Json& interval) {
 return fs::path("data")/token(str(c,"beam"))/token(str(c,"period"))/temperatureDirectory(interval)/("hit_"+integerText(num(c,"hit_tag"))+"_gain_"+integerText(num(c,"gain_tag"))+"_memory_"+integerText(num(c,"memory_cell")))/("settings_"+digest(str(c,"settings")).substr(0,16));
}
inline Json rejectAbsent(const Json& condition,int cell,bool dead) {
 Json r=condition;r["cellid"]=cell;r["accepted"]=0;r["status"]="unresolved";r["reasons"]=dead?"known_dead_channel":"no_led_data_in_condition";return r;
}
struct Settings {
 int minDacs=3; double agreement=.05, spread=.05, temperatureSpan=.5, beamTemperature=.5;
 std::string primary="multi_gaussian"; std::set<int> dead;
 explicit Settings(const Json& j) {
  minDacs=j.at("min_dacs");agreement=j.at("agreement_fraction");spread=j.at("spread_fraction");temperatureSpan=j.at("temperature_span_C");beamTemperature=j.at("beam_temperature_tolerance_C");primary=j.at("primary_method");
  for(auto& x:j.at("dead_channels"))dead.insert(x.get<int>());
 }
};
struct Paths {
 fs::path project,here,input,out; std::string label,outputLabel,comparison;
 Paths(const Options& o):here(fs::canonical("/proc/self/exe").parent_path()) {
  project=o.has("project")?fs::canonical(o.required("project")):here.parent_path().parent_path();
  label=token(o.required("label"));outputLabel=token(o.get("output-label",label));comparison=token(o.get("comparison-name","comparison"));
  input=project/"analysis/result/sipm_gain"/label/comparison;out=project/"analysis/result/sipm_gain_selection"/outputLabel;
 }
};
inline Json sourceChecks(const Paths& p) {
 Json s=Json::object();for(auto name:{"Common.hh","Selection.hh","Select.cc","Merge.cc","Plots.cc","Temperature.cc","Makefile"})s["analysis/sipm_gain_selection/"+std::string(name)]=sha256(p.here/name);
 for(auto name:{"common/GainIO.hh","common/GainMethods.hh","common/Pipeline.hh","common/LEDGroups.hh","compare/CompareIO.hh"})s["analysis/sipm_gain/"+std::string(name)]=sha256(p.here.parent_path()/"sipm_gain"/name);
 return s;
}
inline Json referenceSchema(const fs::path& path) {
 auto f=openRoot(path);auto* t=tree(*f,"T_Event");Json schema=Json::object();
 for(auto* o:*t->GetListOfBranches()){auto* b=dynamic_cast<TBranch*>(o);require(b&&b->GetListOfLeaves()->GetEntries()==1,"Unsupported reference branch");auto* l=dynamic_cast<TLeaf*>(b->GetListOfLeaves()->At(0));require(l&&l->GetLenStatic()==1,"Reference must have scalar branches");schema[b->GetName()]=l->GetTypeName();}
 Json expected={{"cellIDs","Int_t"},{"gain","Double_t"},{"gainError","Double_t"},{"chiSquare","Double_t"},{"NDF","Double_t"},{"temp","Double_t"},{"tempError","Double_t"}};
 require(schema==expected,"Reference T_Event schema differs from the supported seven calibration branches");return schema;
}
inline Json newConfiguration(const Options& o,const Paths& p) {
 Json c={{"format_version",1},{"label",p.label},{"output_label",p.outputLabel},{"comparison_name",p.comparison},{"comparison_path",fs::canonical(p.input).string()},{"comparison_configuration_sha256",sha256(p.input/"inputs/configuration.json")},{"comparison_numeric_marker_sha256",sha256(p.input/"NUMERIC_COMPLETE")},
  {"min_dacs",o.integer("min-dacs",3)},{"agreement_fraction",o.real("agreement",.05)},{"spread_fraction",o.real("max-spread",.05)},{"temperature_span_C",o.real("max-temperature-span",.5)},{"beam_temperature_tolerance_C",o.real("beam-temperature-tolerance",.5)},{"primary_method",o.get("primary-method","multi_gaussian")},{"dead_channels",Json::array()},{"bootstrap",0},{"sources",sourceChecks(p)}};
 require(num(c,"min_dacs")>=2&&num(c,"min_dacs")<=10000,"--min-dacs must be at least 2");
 require(num(c,"agreement_fraction")>0&&num(c,"agreement_fraction")<.5&&num(c,"spread_fraction")>0&&num(c,"spread_fraction")<.5,"Invalid agreement/spread fraction");
 require(num(c,"temperature_span_C")>0&&num(c,"beam_temperature_tolerance_C")>=0,"Invalid temperature tolerance");
 require(str(c,"primary_method")=="multi_gaussian"||str(c,"primary_method")=="envelope_comb","Primary method must be a Gaussian fit");
 if(o.has("dead-channels")) {
  auto path=fs::canonical(o.required("dead-channels"));std::ifstream in(path);require(bool(in),"Cannot read dead-channel list");std::string line;std::set<int> ids;
  while(std::getline(in,line)){line=trim(line.substr(0,line.find('#')));if(line.empty()||line=="cellid")continue;size_t used=0;int id=std::stoi(line,&used);require(used==line.size(),"Dead-channel list must contain one physical CellID per line");physicalCell({{"cellid",id}});require(ids.insert(id).second,"Duplicate dead-channel CellID");}
  for(int id:ids)c["dead_channels"].push_back(id);
  c["dead_channel_file"]=fingerprint(path);c["dead_channel_sha256"]=sha256(path);
 }
 auto ref=fs::canonical(o.get("reference-root",DefaultReference));c["reference"]=fingerprint(ref);c["reference_sha256"]=sha256(ref);c["calibration_schema"]=referenceSchema(ref);
 return c;
}
inline Json configuration(const Options& o,const Paths& p) {
 require(fs::is_regular_file(p.input/"NUMERIC_COMPLETE"),"Complete comparison merge first");
 auto summary=readJson(p.input/"summary.json");require(!summary.value("partial",true),"Partial comparison cannot be exported as calibration candidates");
 auto complete=readJson(p.input/"NUMERIC_COMPLETE");
 require(complete.at("files").at("summary.json")==sha256(p.input/"summary.json"),"Comparison summary changed");
 require(complete.at("configuration_sha256")==sha256(p.input/"inputs/configuration.json"),"Comparison configuration changed");
 fs::create_directories(p.out/"inputs");Lock lock(p.out/".configuration.lock",true);auto path=p.out/"inputs/configuration.json";
 if(!fs::exists(path)){
  require(o.get("stage","select")=="select","Run selection before merge/plots");auto c=newConfiguration(o,p);
  for(auto it=c["sources"].begin();it!=c["sources"].end();++it){auto dest=p.out/"inputs/source"/it.key();fs::create_directories(dest.parent_path());fs::copy_file(p.project/it.key(),dest,fs::copy_options::overwrite_existing);require(sha256(dest)==it.value(),"Source changed while snapshotting");}
  atomicJson(path,c);
 }
 auto c=readJson(path);require(c.at("format_version")==1&&str(c,"label")==p.label&&str(c,"output_label")==p.outputLabel&&str(c,"comparison_name")==p.comparison,"Selection input/label differs");
 require(str(c,"comparison_path")==fs::canonical(p.input).string(),"Comparison path differs");
 if(o.get("stage")=="plots"||o.get("stage")=="temperature-paths") {
  // Rendering consumes sealed numerical results. Keep their original source
  // snapshot intact; the renderer records its own sources with the plots.
  for(auto it=c.at("sources").begin();it!=c.at("sources").end();++it)
   require(sha256(p.out/"inputs/source"/it.key())==it.value(),"Frozen numerical source changed: "+it.key());
 } else require(c.at("sources")==sourceChecks(p),"Selection source changed; use a new --output-label");
 require(c.at("comparison_configuration_sha256")==sha256(p.input/"inputs/configuration.json")&&c.at("comparison_numeric_marker_sha256")==sha256(p.input/"NUMERIC_COMPLETE"),"Comparison changed; use a new --output-label");
 const std::map<std::string,std::string> numeric={{"min-dacs","min_dacs"},{"agreement","agreement_fraction"},{"max-spread","spread_fraction"},{"max-temperature-span","temperature_span_C"},{"beam-temperature-tolerance","beam_temperature_tolerance_C"}};
 for(auto& [option,key]:numeric)if(o.has(option))require(o.real(option,0)==num(c,key),"Selection option differs: "+option);
 if(o.has("primary-method"))require(o.get("primary-method")==str(c,"primary_method"),"Primary method differs");
 if(o.has("reference-root"))require(sha256(o.required("reference-root"))==c.at("reference_sha256"),"Reference differs");
 if(o.has("dead-channels"))require(c.contains("dead_channel_sha256")&&sha256(o.required("dead-channels"))==c.at("dead_channel_sha256"),"Dead-channel mask differs");
 return c;
}
inline Json metadataInfo(const Json& config,const std::vector<std::string>& command) {
 return {{"configuration",config},{"command",commandLine(command)},{"gainError_definition","max(median conditional fit_error, max within-condition gain deviation, max stability_ADC, max supported-method difference); not a statistical standard error"},{"chiSquare_definition","Poisson deviance of the representative accepted LED fit; not a chi-square of the median"},{"NDF_definition","ndf of the same representative LED fit"},{"tempError_definition","max observed temperature deviation and input temperature RMS; not a standard error"},{"settings_note","Unknown settings remain flagged; no temperature correction or cross-period pooling"},{"calibration_deployment",false}};
}
inline std::string configHash(const Paths& p){return sha256(p.out/"inputs/configuration.json");}
void selectLayer(const Paths&,const Json&,int,const std::vector<std::string>&);
void merge(const Paths&,const Json&,const std::vector<std::string>&);
double temperatureWidth(const Paths&,const Json&);
void migrateTemperaturePaths(const Paths&,const Json&,const std::vector<std::string>&);
void plots(const Paths&,const Json&,int,int,const std::string&,const std::vector<std::string>&);
}
