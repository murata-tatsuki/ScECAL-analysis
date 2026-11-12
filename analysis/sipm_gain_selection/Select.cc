#include "Selection.hh"

namespace selected_gain {
void selectLayer(const Paths& p,const Json& config,int layer,const std::vector<std::string>& command) {
 const auto name=layerName(layer);auto input=p.input/"layers"/name;
 Lock inputLock(input/".lock",false,true);
 const auto inputSHA=sha256(input/"NUMERIC_COMPLETE");auto merged=readJson(p.input/"NUMERIC_COMPLETE");
 require(merged.at("layers").at(name)==inputSHA,"Comparison layer differs from completed merge: "+name);
 auto info=readJson(input/"NUMERIC_COMPLETE");
 require(info.at("configuration_sha256")==config.at("comparison_configuration_sha256"),"Comparison layer configuration differs");
 require(info.at("files").at("numeric.root")==sha256(input/"numeric.root"),"Comparison layer ROOT changed");
 auto out=p.out/"layers"/name;fs::create_directories(out);Lock lock(out/".lock");auto signature=configHash(p);
 if(fs::exists(out/"SELECTION_COMPLETE")){
  auto saved=verifyComplete(out,"SELECTION_COMPLETE",signature);require(saved.at("input_marker_sha256")==inputSHA,"Selection input changed");
  std::cout<<"Already complete: "<<name<<std::endl;return;
 }
 Settings settings(config);auto source=openRoot(input/"numeric.root");auto* t=tree(*source,"Gain");Reader reader(*t,gainSchema());
 auto target=openRoot(out/"selection.root","RECREATE");RootRows decisions(*target,"Decisions",DecisionSchema),evidenceRows(*target,"Evidence",EvidenceSchema),beamRows(*target,"BeamChecks",BeamSchema);
 TsvRows decisionTsv(out/"decisions.tsv",DecisionSchema.columns()),evidenceTsv(out/"evidence.tsv",EvidenceSchema.columns()),beamTsv(out/"beam_checks.tsv",BeamSchema.columns());
 long long accepted=0;std::set<int> channels,acceptedChannels;std::set<std::string> conditions;
 auto consume=[&](const Json& rows){
  for(auto& group:selectChannel(rows,settings)){
   const auto& d=group.decision;decisions.add(d);decisionTsv.add(d);conditions.insert(str(d,"condition"));channels.insert(physicalCell(d));
   if(num(d,"accepted")==1){++accepted;acceptedChannels.insert(physicalCell(d));}
   for(auto& e:group.evidence){evidenceRows.add(e);evidenceTsv.add(e);}
   for(auto& b:group.beam){beamRows.add(b);beamTsv.add(b);}
  }
 };
 int previous=-1;Json rows=Json::array();
 for(Long64_t i=0;i<t->GetEntries();++i){auto r=reader.at(i);int cell=physicalCell(r);require(cell/100000==layer&&cell>=previous,"Input Gain must be sorted by physical channel within layer");
  if(cell!=previous&&!rows.empty()){consume(rows);rows=Json::array();}
  previous=cell;rows.push_back(std::move(r));
 }
 if(!rows.empty())consume(rows);
 Json stats={{"layer",layer},{"input_rows",t->GetEntries()},{"channel_conditions",decisions.size()},{"selected_channel_conditions",accepted},{"unresolved_channel_conditions",decisions.size()-accepted},{"led_channels",channels.size()},{"selected_unique_channels",acceptedChannels.size()},{"conditions",conditions.size()},{"evidence_rows",evidenceRows.size()},{"beam_checks",beamRows.size()},{"input_marker_sha256",inputSHA},{"command",commandLine(command)}};
 decisions.finish();evidenceRows.finish();beamRows.finish();decisionTsv.finish();evidenceTsv.finish();beamTsv.finish();
 auto meta=metadataInfo(config,command);meta["summary"]=stats;metadata(*target,meta);closeRoot(*target);
 atomicJson(out/"summary.json",stats);
 seal(out,"SELECTION_COMPLETE",signature,{"selection.root","decisions.tsv","evidence.tsv","beam_checks.tsv","summary.json"},stats);
 std::cout<<"Completed selection: "<<name<<", selected "<<accepted<<" / "<<decisions.size()<<" channel-conditions ("<<acceptedChannels.size()<<" unique channels)"<<std::endl;
}
}

int main(int argc,char** argv)try{
 using namespace selected_gain;
 const std::set<std::string> allowed={"label","output-label","comparison-name","project","stage","layer","reference-root","dead-channels","primary-method","min-dacs","agreement","max-spread","max-temperature-span","beam-temperature-tolerance","plot-limit","plot-scope"};
 Options options;std::vector<std::string> command;for(int i=0;i<argc;++i)command.emplace_back(argv[i]);
 for(int i=1;i<argc;++i){std::string key=argv[i];if(key=="--help"){
   std::cout<<"Select --label LABEL --stage select|merge|plots|temperature-paths [--layer 0..29]\n"
     "  --output-label NAME         separate output conditions / validation\n"
     "  --comparison-name NAME     default: comparison\n"
     "  --primary-method METHOD    multi_gaussian (default) or envelope_comb\n"
     "  --min-dacs N               >= 2, default 3 (2 is a relaxed selection)\n"
     "  --agreement FRACTION       default 0.05, symmetric method difference\n"
     "  --max-spread FRACTION      default 0.05, max gain deviation / median\n"
     "  --max-temperature-span C   default 0.5\n"
     "  --beam-temperature-tolerance C  default 0.5 (diagnostics only)\n"
     "  --dead-channels FILE       one physical CellID per line\n"
     "  --reference-root FILE      legacy seven-branch schema reference\n"
     "  --plot-limit N             max channel plots per layer, 0 = all\n"
     "  --plot-scope selected|all   default: selected (all overviews retained)\n"
     "  --project DIRECTORY        project root\n"
     "Stages are separate; merge requires all 30 completed layers.\n";return 0;
  }
  require(key.rfind("--",0)==0,"Expected --option");key=key.substr(2);require(allowed.count(key),"Unknown option: "+key);require(!options.has(key),"Duplicate option: "+key);require(i+1<argc&&std::string(argv[i+1]).rfind("--",0)!=0,"Missing value: "+key);options.values[key]=argv[++i];
 }
 const auto stage=options.required("stage");require(stage=="select"||stage=="merge"||stage=="plots"||stage=="temperature-paths","Unknown stage");
 require((stage!="merge"&&stage!="temperature-paths")||!options.has("layer"),"Merge requires all layers; omit --layer");
 require(stage=="plots"||!options.has("plot-limit"),"--plot-limit is only for plots");
 require(stage=="plots"||!options.has("plot-scope"),"--plot-scope is only for plots");
 int first=0,last=30;if(options.has("layer")){auto l=options.integer("layer",-1);require(l>=0&&l<30,"Layer must be 0..29");first=l;last=first+1;}
 auto limit=options.integer("plot-limit",0);require(limit>=0&&limit<=1000000,"Invalid plot limit");
 auto scope=options.get("plot-scope","selected");require(scope=="selected"||scope=="all","Invalid plot scope");
 Paths paths(options);require(fs::is_directory(paths.input),"Comparison directory not found");
 Lock inputLock(paths.input/".lock",false,true);fs::create_directories(paths.out);
 // A merge excludes selection/plots; independent layer jobs may run together.
 Lock outputLock(paths.out/".lock",false,stage!="merge"&&stage!="temperature-paths");auto config=configuration(options,paths);
 gROOT->SetBatch(kTRUE);
 if(stage=="select")for(int l=first;l<last;++l)selectLayer(paths,config,l,command);
 if(stage=="merge")merge(paths,config,command);
 if(stage=="temperature-paths")migrateTemperaturePaths(paths,config,command);
 if(stage=="plots")for(int l=first;l<last;++l)plots(paths,config,l,limit,scope,command);
 return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
