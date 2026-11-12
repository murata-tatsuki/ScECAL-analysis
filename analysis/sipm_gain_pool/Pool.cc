#include "../sipm_gain_selection/Common.hh"
#include <TH1D.h>

using namespace selected_gain;
namespace {
const std::vector<std::string> Fields={"sid","run","cellid","hit_tag","gain_tag","temperature_bin","memory_cell","time_bin","entries","valid_temperature_entries","overflow","temperature_C","temperature_RMS_C","pooled_runs"};
const Schema PoolSchema{{},Fields};
const Schema ContributionSchema{{"sample","beam","kind","energy_GeV","source_root"},{"sid","source_sid","cellid","temperature_bin","entries","valid_temperature_entries","temperature_C","temperature_RMS_C","overflow"}};
struct Locations {
 fs::path project,source,out,fit;std::string from,label;
 Locations(const Options& o){project=o.has("project")?fs::canonical(o.required("project")):fs::canonical("/proc/self/exe").parent_path().parent_path().parent_path();from=token(o.required("from-label"));label=token(o.required("label"));require(from!=label,"Use a new pooled label");source=project/"analysis/result/sipm_gain"/from;fit=project/"analysis/result/sipm_gain"/label;out=project/"analysis/result/sipm_gain_pool"/label;}
};
Json sources(const Locations& p){Json j=Json::object();for(const auto& name:{"analysis/sipm_gain_pool/Pool.cc","analysis/sipm_gain_pool/Makefile","analysis/sipm_gain_selection/Common.hh","analysis/sipm_gain/compare/CompareIO.hh","analysis/sipm_gain/common/GainIO.hh","analysis/sipm_gain/common/GainMethods.hh","analysis/sipm_gain/common/Pipeline.hh","analysis/sipm_gain/common/LEDGroups.hh"})j[name]=sha256(p.project/name);return j;}
std::set<int> integers(const std::string& text,int maximum){std::set<int> values;if(text=="all"){for(int i=0;i<maximum;++i)values.insert(i);return values;}for(const auto& s:split(text)){size_t n=0;int v=std::stoi(s,&n);require(n==s.size()&&v>=0&&v<maximum,"Invalid integer list");require(values.insert(v).second,"Repeated integer");}require(!values.empty(),"Empty integer list");return values;}
Json pooledRow(const Locations& p,const std::string& beam,const std::string& kind){auto sample=beam+"_"+kind+"_temperature_pool";return {{"sample",sample},{"kind",kind},{"beam",beam},{"energy_GeV","mixed"},{"period","temperature_pool"},{"led_group","-1"},{"dac","-1"},{"settings","user_confirmed_common_within_"+beam},{"raw",(p.fit/"data"/kind/sample/"collection/spectra.root").string()}};}
void prepare(const Locations& p,const Options& o,const std::vector<std::string>& command){
 require(o.get("settings-policy")=="common-within-beam","Declare --settings-policy common-within-beam only after confirming bias/ASIC settings within each beam");
 fs::create_directories(p.out/"inputs");Lock lock(p.out/".lock");
 if(fs::exists(p.out/"PREPARE_COMPLETE")){auto cfg=readJson(p.out/"inputs/configuration.json");verifyComplete(p.out,"PREPARE_COMPLETE",sha256(p.out/"inputs/configuration.json"));require(cfg.at("sources")==sources(p),"Pooling source changed; use a new label");
  for(const auto& k:{"beam","kind","layers","channels"}){if(o.has(k))require(str(cfg,k)==o.get(k),"Pooling option differs");}
  if(o.has("min-runs"))require(num(cfg,"min_runs")==o.integer("min-runs",2),"min-runs differs");std::cout<<"Already prepared"<<std::endl;return;}
 require(!fs::exists(p.fit)||(!fs::exists(p.fit/"data")&&fs::exists(p.out/"inputs/configuration.json")),"Output fit label already exists; choose a new label");
 auto beam=o.get("beam","all"),kind=o.get("kind","all"),layers=o.get("layers","all"),channels=o.get("channels","all");
 require(beam=="all"||beam=="ps"||beam=="sps","Invalid beam");require(kind=="all"||kind=="electron"||kind=="muon","Invalid particle");integers(layers,30);
 if(channels!="all")for(int cell:integers(channels,3000000))physicalCell({{"cellid",cell}});
 int minRuns=o.integer("min-runs",2);require(minRuns>=2,"At least two runs are required");
 Json input=Json::array(),rows=Json::array();std::set<std::pair<std::string,std::string>> groups;double width=NaN;
 for(const auto& unit:readJson(p.source/"inputs/analysis_units.json")){const auto b=str(unit,"beam"),k=str(unit,"kind");if(k=="led"||(beam!="all"&&b!=beam)||(kind!="all"&&k!=kind))continue;
  require((b=="ps"||b=="sps")&&(k=="electron"||k=="muon"),"Unexpected source group");
  auto directory=p.source/"data"/k/token(str(unit,"sample"))/"collection";Lock sourceLock(directory/".lock",false,true);
  auto meta=readJson(directory/"inputs/collection.json");require(meta.at("manifest_row")==unit,"Source manifest changed");require(num(meta.at("collection"),"memory_cell_selection")==0,"Source memory selection differs");
  double w=num(meta.at("collection"),"temperature_bin");require(std::isfinite(w)&&w>0&&(!std::isfinite(width)||width==w),"Different temperature grids cannot be merged");width=w;
  auto file=directory/"spectra.root";auto stamp=fingerprint(file);auto checksum=sha256(file);require(checksum==trim(readText(directory/"COLLECT_COMPLETE"))&&stamp==fingerprint(file),"Source collection changed or incomplete");
  input.push_back({{"unit",unit},{"root",file.string()},{"fingerprint",stamp},{"root_sha256",checksum},{"marker_sha256",sha256(directory/"COLLECT_COMPLETE")},{"metadata_sha256",sha256(directory/"inputs/collection.json")}});groups.emplace(b,k);
  std::cout<<"Verified "<<str(unit,"sample")<<std::endl;
 }
 require(!input.empty(),"No matching input runs");for(const auto& [b,k]:groups)rows.push_back(pooledRow(p,b,k));
 Json cfg={{"format_version",1},{"from_label",p.from},{"label",p.label},{"beam",beam},{"kind",kind},{"layers",layers},{"channels",channels},{"min_runs",minRuns},{"temperature_width_C",width},{"settings_policy","user_confirmed_common_within_beam"},{"energy_pooling","all energies within each particle and beam; energy retained in Contributions"},{"adc_alignment","uncorrected ADC bins; no run-dependent pedestal shift applied"},{"input",input},{"manifest",rows},{"sources",sources(p)},{"command",commandLine(command)}};
 atomicJson(p.out/"inputs/configuration.json",cfg);
 std::vector<std::string> files={"inputs/configuration.json"};for(auto it=cfg["sources"].begin();it!=cfg["sources"].end();++it){auto rel=fs::path("inputs/source")/it.key();fs::create_directories((p.out/rel).parent_path());if(!fs::exists(p.out/rel))fs::copy_file(p.project/it.key(),p.out/rel);require(sha256(p.out/rel)==it.value(),"Source changed while snapshotting");files.push_back(rel.string());}
 fs::create_directories(p.fit/"inputs");writeTable(p.fit/"inputs/manifest.tsv",rows,ManifestFields);writeJson(p.fit/"inputs/analysis_units.json",analysisUnits(rows));
 seal(p.out,"PREPARE_COMPLETE",sha256(p.out/"inputs/configuration.json"),files,{{"input_runs",input.size()},{"groups",groups.size()}});
 std::cout<<"Prepared "<<input.size()<<" input runs in "<<groups.size()<<" beam/particle groups"<<std::endl;
}
Json configuration(const Locations& p){auto c=readJson(p.out/"inputs/configuration.json");verifyComplete(p.out,"PREPARE_COMPLETE",sha256(p.out/"inputs/configuration.json"));require(str(c,"from_label")==p.from&&str(c,"label")==p.label&&c.at("sources")==sources(p),"Pooling source/label changed");return c;}
fs::path layerPath(const Locations& p,const std::string& beam,const std::string& kind,int layer){return p.out/"data"/beam/kind/layerName(layer);}
void checkGroup(const Json& cfg,const std::string& beam,const std::string& kind){bool found=false;for(const auto& row:cfg.at("manifest"))found|=str(row,"beam")==beam&&str(row,"kind")==kind;require(found,"Group is not in prepared configuration");}
void collect(const Locations& p,const Json& cfg,const Options& o,const std::vector<std::string>& command){
 const auto beam=o.required("beam"),kind=o.required("kind");checkGroup(cfg,beam,kind);const int layer=o.integer("layer",-1);require(integers(str(cfg,"layers"),30).count(layer),"Layer is not in prepared selection");
 auto out=layerPath(p,beam,kind,layer);fs::create_directories(out);Lock lock(out/".lock");const auto hash=sha256(p.out/"inputs/configuration.json");if(fs::exists(out/"POOL_COMPLETE")){verifyComplete(out,"POOL_COMPLETE",hash);std::cout<<"Already pooled "<<out<<std::endl;return;}
 std::set<int> cells;if(str(cfg,"channels")!="all")cells=integers(str(cfg,"channels"),3000000);
 struct Acc {Vec counts=Vec(4096,0.);double entries=0,nt=0,sumT=0,sumT2=0,overflow=0;std::set<std::string> runs;Json contributions=Json::array();};
 std::map<std::pair<int,int>,Acc> accum;long long inputSpectra=0;double inputCounts=0;
 for(const auto& item:cfg.at("input")){const auto& unit=item.at("unit");if(str(unit,"beam")!=beam||str(unit,"kind")!=kind)continue;fs::path file=str(item,"root");auto directory=file.parent_path();Lock inputLock(directory/".lock",false,true);
  require(fingerprint(file)==item.at("fingerprint")&&sha256(directory/"COLLECT_COMPLETE")==item.at("marker_sha256")&&sha256(directory/"inputs/collection.json")==item.at("metadata_sha256"),"Input changed since prepare");
  auto root=openRoot(file);auto* treeIn=tree(*root,"Spectra");
  for(Long64_t i=0;i<treeIn->GetEntries();++i){require(treeIn->GetEntry(i)>0,"Cannot read source spectrum metadata");int cell=leaf(*treeIn,"cellid");if(cell/100000!=layer||(!cells.empty()&&!cells.count(cell)))continue;
   if(leaf(*treeIn,"memory_cell")!=0||leaf(*treeIn,"hit_tag")!=1||leaf(*treeIn,"gain_tag")!=1)continue;
   int bin=leaf(*treeIn,"temperature_bin");double nt=leaf(*treeIn,"valid_temperature_entries"),temp=leaf(*treeIn,"temperature_C"),rms=leaf(*treeIn,"temperature_RMS_C"),width=num(cfg,"temperature_width_C");
   if(bin<0||nt<=0||!std::isfinite(temp)||!std::isfinite(rms))continue;
   require(temp>=bin*width-1e-7&&temp<(bin+1)*width+1e-7,"Temperature does not match source interval");
   int sid=leaf(*treeIn,"sid");std::unique_ptr<TH1D> h(dynamic_cast<TH1D*>(root->Get(("s"+std::to_string(sid)).c_str())));require(h&&h->GetNbinsX()==4096&&h->GetXaxis()->GetXmin()==-.5&&h->GetXaxis()->GetXmax()==4095.5,"Unexpected ADC grid");h->SetDirectory(nullptr);if(h->Integral()<=0)continue;
   auto& a=accum[{cell,bin}];for(int j=0;j<4096;++j){double v=h->GetBinContent(j+1);require(std::isfinite(v)&&v>=0,"Invalid ADC counts");a.counts[j]+=v;inputCounts+=v;}
   double entries=leaf(*treeIn,"entries"),overflow=leaf(*treeIn,"overflow");a.entries+=entries;a.nt+=nt;a.sumT+=nt*temp;a.sumT2+=nt*(rms*rms+temp*temp);a.overflow+=overflow;a.runs.insert(str(unit,"sample"));++inputSpectra;
   a.contributions.push_back({{"sample",str(unit,"sample")},{"beam",beam},{"kind",kind},{"energy_GeV",str(unit,"energy_GeV")},{"source_root",file.string()},{"source_sid",sid},{"cellid",cell},{"temperature_bin",bin},{"entries",entries},{"valid_temperature_entries",nt},{"temperature_C",temp},{"temperature_RMS_C",rms},{"overflow",overflow}});
  }
  require(fingerprint(file)==item.at("fingerprint"),"Input changed while pooling");
 }
 auto root=openRoot(out/"spectra.root","RECREATE");RootRows spectra(*root,"Spectra",PoolSchema),contributions(*root,"Contributions",ContributionSchema);Json index=Json::array();double keptCounts=0,droppedCounts=0;int sid=0,dropped=0;
 for(auto& [key,a]:accum){double count=std::accumulate(a.counts.begin(),a.counts.end(),0.);if(a.runs.size()<size_t(num(cfg,"min_runs"))){droppedCounts+=count;++dropped;continue;}auto [cell,bin]=key;double temp=a.sumT/a.nt;
  Json r={{"sid",sid},{"run",-1},{"cellid",cell},{"hit_tag",1},{"gain_tag",1},{"temperature_bin",bin},{"memory_cell",0},{"time_bin",0},{"entries",a.entries},{"valid_temperature_entries",a.nt},{"overflow",a.overflow},{"temperature_C",temp},{"temperature_RMS_C",std::sqrt(std::max(0.,a.sumT2/a.nt-temp*temp))},{"pooled_runs",a.runs.size()}};spectra.add(r);
  TH1D h(("s"+std::to_string(sid)).c_str(),";HG ADC;entries",4096,-.5,4095.5);h.SetDirectory(nullptr);for(int j=0;j<4096;++j){h.SetBinContent(j+1,a.counts[j]);h.SetBinError(j+1,std::sqrt(a.counts[j]));}h.SetEntries(a.entries);root->cd();require(h.Write()>0,"Cannot save pooled histogram");
  for(auto& x:a.contributions){x["sid"]=sid;contributions.add(x);}auto interval=temperatureInterval(bin,num(cfg,"temperature_width_C"));r["temperature_lower_C"]=interval["lower_C"];r["temperature_upper_C"]=interval["upper_C"];index.push_back(r);keptCounts+=count;++sid;
 }
 require(inputCounts==keptCounts+droppedCounts,"ADC count conservation failed");spectra.finish();contributions.finish();Json info={{"configuration_sha256",hash},{"beam",beam},{"kind",kind},{"layer",layer},{"input_spectra",inputSpectra},{"pooled_spectra",sid},{"insufficient_run_spectra",dropped},{"input_adc_counts",inputCounts},{"pooled_adc_counts",keptCounts},{"excluded_adc_counts",droppedCounts},{"command",commandLine(command)}};metadata(*root,info);closeRoot(*root);writeTable(out/"spectra.tsv",index);atomicJson(out/"summary.json",info);seal(out,"POOL_COMPLETE",hash,{"spectra.root","spectra.tsv","summary.json"},info);std::cout<<"Pooled "<<sid<<" spectra for "<<beam<<" "<<kind<<" "<<layerName(layer)<<std::endl;
}
void mergePool(const Locations& p,const Json& cfg,const Options& o,const std::vector<std::string>& command){
 const auto beam=o.required("beam"),kind=o.required("kind");checkGroup(cfg,beam,kind);const auto hash=sha256(p.out/"inputs/configuration.json");auto unit=analysisUnits(Json::array({pooledRow(p,beam,kind)})).at(0);auto directory=p.fit/"data"/kind/str(unit,"sample")/"collection";fs::create_directories(directory);Lock outputLock(directory/".lock");auto marker=p.out/"data"/beam/kind;Json checks=Json::object();
 auto layers=integers(str(cfg,"layers"),30);for(int l:layers){auto input=layerPath(p,beam,kind,l);verifyComplete(input,"POOL_COMPLETE",hash);checks[layerName(l)]=sha256(input/"POOL_COMPLETE");}
 if(fs::exists(marker/"MERGE_COMPLETE")){auto saved=verifyComplete(marker,"MERGE_COMPLETE",hash);require(saved.at("layers")==checks,"Pooling layers changed");require(trim(readText(directory/"COLLECT_COMPLETE"))==sha256(directory/"spectra.root"),"Final collection changed");std::cout<<"Already merged "<<beam<<" "<<kind<<std::endl;return;}
 require(!fs::exists(directory/"COLLECT_COMPLETE"),"Existing untracked final collection");auto output=openRoot(directory/"spectra.root","RECREATE");RootRows spectra(*output,"Spectra",PoolSchema),contributions(*output,"Contributions",ContributionSchema);int next=0;
 for(int l:layers){auto input=openRoot(layerPath(p,beam,kind,l)/"spectra.root");auto* t=tree(*input,"Spectra");Reader reader(*t,PoolSchema);std::map<int,int> sids;
  for(Long64_t i=0;i<t->GetEntries();++i){auto r=reader.at(i);int old=num(r,"sid"),sid=next++;require(sids.emplace(old,sid).second,"Duplicate source sid");r["sid"]=sid;spectra.add(r);std::unique_ptr<TH1D> h(dynamic_cast<TH1D*>(input->Get(("s"+std::to_string(old)).c_str())));require(bool(h),"Missing pooled histogram");h->SetDirectory(nullptr);h->SetName(("s"+std::to_string(sid)).c_str());output->cd();require(h->Write()>0,"Cannot merge pooled histogram");}
  auto* c=tree(*input,"Contributions");Reader cr(*c,ContributionSchema);for(Long64_t i=0;i<c->GetEntries();++i){auto r=cr.at(i);r["sid"]=sids.at(int(num(r,"sid")));contributions.add(r);}
 }
 spectra.finish();contributions.finish();metadata(*output,{{"pool_configuration",cfg},{"layers",checks},{"command",commandLine(command)}});closeRoot(*output);
 Json record={{"manifest_row",unit},{"collection",{{"memory_cell_selection",0},{"temperature_bin",num(cfg,"temperature_width_C")},{"time_bin",0},{"bootstrap_storage",false},{"channels",str(cfg,"channels")},{"events_per_file",0}}},{"pool_configuration_sha256",hash},{"pool_configuration_path",(p.out/"inputs/configuration.json").string()},{"pool_layers",checks},{"min_pooled_runs",cfg.at("min_runs")},{"settings_policy",cfg.at("settings_policy")},{"command",commandLine(command)}};
 atomicJson(directory/"inputs/collection.json",record);writeText(directory/"COLLECT_COMPLETE",sha256(directory/"spectra.root")+'\n');
 std::vector<std::string> outputs;for(auto name:{"spectra.root","inputs/collection.json","COLLECT_COMPLETE"})outputs.push_back(fs::relative(directory/name,marker).string());seal(marker,"MERGE_COMPLETE",hash,outputs,{{"layers",checks},{"spectra",next}});std::cout<<"Merged "<<next<<" spectra; ready for the existing Analyze fit stage"<<std::endl;
}
}
int main(int argc,char** argv)try{
 Options o;std::vector<std::string> command;for(int i=0;i<argc;++i)command.emplace_back(argv[i]);const std::set<std::string> allowed={"stage","from-label","label","project","settings-policy","beam","kind","layer","layers","channels","min-runs"};
 for(int i=1;i<argc;++i){std::string key=argv[i];if(key=="--help"){std::cout<<"Pool --from-label full_grouped --label temperature_pool --stage prepare|collect|merge\nprepare: --settings-policy common-within-beam [--beam ps|sps --kind electron|muon --layers all --channels all --min-runs 2]\ncollect: --beam ps|sps --kind electron|muon --layer 0..29\nmerge: --beam ps|sps --kind electron|muon\nEach stage is run separately. Fits use analysis/sipm_gain/Analyze after all merges.\n";return 0;}require(key.rfind("--",0)==0&&allowed.count(key.substr(2)),"Unknown option");key=key.substr(2);require(!o.has(key)&&i+1<argc&&std::string(argv[i+1]).rfind("--",0)!=0,"Duplicate option or missing value");o.values[key]=argv[++i];}
 Locations p(o);gROOT->SetBatch(kTRUE);auto stage=o.required("stage");if(stage=="prepare")require(!o.has("layer"),"prepare uses --layers, not --layer");else for(auto key:{"settings-policy","layers","channels","min-runs"})require(!o.has(key),"Set pooling options during prepare");if(stage=="merge")require(!o.has("layer"),"merge includes every prepared layer");if(stage=="prepare")prepare(p,o,command);else{require(stage=="collect"||stage=="merge","Unknown stage");Lock global(p.out/".lock",false,true);auto c=configuration(p);if(stage=="collect")collect(p,c,o,command);else mergePool(p,c,o,command);}return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
