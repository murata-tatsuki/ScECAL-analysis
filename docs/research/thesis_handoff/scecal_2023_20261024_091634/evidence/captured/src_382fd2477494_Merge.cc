#include "Common.hh"

namespace selected_gain {
void merge(const Paths& p,const Json& config,const std::vector<std::string>& command) {
 auto signature=configHash(p);Json layerChecks=Json::object();std::vector<fs::path> inputs;
 for(int layer=0;layer<30;++layer){auto dir=p.out/"layers"/layerName(layer);
  require(fs::exists(dir/"SELECTION_COMPLETE"),"Selection incomplete: "+layerName(layer));
  verifyComplete(dir,"SELECTION_COMPLETE",signature);layerChecks[layerName(layer)]=sha256(dir/"SELECTION_COMPLETE");inputs.push_back(dir/"selection.root");
 }
 if(fs::exists(p.out/"MERGE_COMPLETE")){
  auto saved=verifyComplete(p.out,"MERGE_COMPLETE",signature);require(saved.at("layers")==layerChecks,"Selection layers changed");
  std::cout<<"Already complete: merge"<<std::endl;return;
 }
 auto combined=openRoot(p.out/"selection.root","RECREATE");
 for(size_t k=0;k<OutputTrees.size();++k){TChain chain(OutputTrees[k].c_str());for(auto& path:inputs)require(chain.Add(path.c_str())==1,"Cannot add selection layer");combined->cd();
  if(chain.GetEntries()==0){RootRows empty(*combined,OutputTrees[k].c_str(),OutputSchemas[k]);empty.finish();}
  else {auto* t=chain.CloneTree(-1,"fast");require(t&&t->GetEntries()==chain.GetEntries(),"Selection ROOT merge failed");require(t->Write()>0,"Cannot write merged selection tree");t->SetDirectory(nullptr);delete t;}
 }
 metadata(*combined,metadataInfo(config,command));closeRoot(*combined);combined.reset();
 auto source=openRoot(p.out/"selection.root");auto* t=tree(*source,"Decisions");Reader reader(*t,DecisionSchema);
 // Keep only entry indices for the global inventory, not all per-spectrum evidence.
 std::map<std::string,std::map<int,Long64_t>> indices;std::map<std::string,Json> conditions;
 for(Long64_t i=0;i<t->GetEntries();++i){auto row=reader.at(i);int cell=physicalCell(row);auto c=conditionFields(row);auto key=str(c,"condition");
  require(key==str(row,"condition"),"Condition identifier mismatch");require(indices[key].emplace(cell,i).second,"Duplicate channel within a calibration condition");conditions[key]=c;
 }
 require(!conditions.empty(),"No LED conditions to export");Settings settings(config);Json manifest=Json::array();std::vector<std::string> outputs={"selection.root","conditions.tsv","summary.json"};long long totalSelected=0,totalRejected=0;
 for(auto& [key,c]:conditions){auto interval=temperatureInterval(num(c,"temperature_bin"),temperatureWidth(p,c));auto relativePath=conditionPath(c,interval)/"gain.root";auto path=p.out/relativePath;fs::create_directories(path.parent_path());auto file=openRoot(path,"RECREATE");
  Int_t cellIDs=0;Double_t gain=NaN,gainError=NaN,chiSquare=NaN,NDF=NaN,temp=NaN,tempError=NaN;
  std::string gainErrorKind="conservative_quality_scale_not_statistical",qualityStatistic="representative_fit_poisson_deviance";
  TTree calibration("T_Event","Selected LED single-photoelectron gain; see metadata_json for uncertainty and quality definitions");calibration.SetDirectory(file.get());
  calibration.Branch("cellIDs",&cellIDs,"cellIDs/I");calibration.Branch("gain",&gain,"gain/D");calibration.Branch("gainError",&gainError,"gainError/D");calibration.Branch("chiSquare",&chiSquare,"chiSquare/D");calibration.Branch("NDF",&NDF,"NDF/D");calibration.Branch("temp",&temp,"temp/D");calibration.Branch("tempError",&tempError,"tempError/D");
  calibration.Branch("gainErrorKind",&gainErrorKind);calibration.Branch("qualityStatistic",&qualityStatistic);
  RootRows selected(*file,"Selection",DecisionSchema),rejected(*file,"Rejected",DecisionSchema);
  for(int layer=0;layer<30;++layer)for(int chip=0;chip<6;++chip)for(int channel=0;channel<(chip==5?30:36);++channel){int cell=layer*100000+chip*10000+channel;auto found=indices[key].find(cell);
   auto row=found==indices[key].end()?rejectAbsent(c,cell,settings.dead.count(cell)):reader.at(found->second);
   if(num(row,"accepted")!=1){rejected.add(row);continue;}
   require(!settings.dead.count(cell),"Dead channel cannot enter calibration");
   gain=num(row,"gain");gainError=num(row,"gain_error");chiSquare=num(row,"representative_deviance");NDF=num(row,"representative_ndf");temp=num(row,"temperature_C");tempError=num(row,"temperature_error_C");
   require(std::isfinite(gain)&&gain>0&&std::isfinite(gainError)&&gainError>0&&std::isfinite(chiSquare)&&chiSquare>=0&&std::isfinite(NDF)&&NDF>0&&std::isfinite(temp)&&std::isfinite(tempError)&&tempError>=0,"Invalid selected calibration row");
   cellIDs=layer*1000000+chip*10000+channel;
   require(calibration.Fill()>0,"Cannot write calibration entry");selected.add(row);
  }
  file->cd();require(calibration.Write()>0,"Cannot save T_Event");calibration.SetDirectory(nullptr);selected.finish();rejected.finish();
  require(selected.size()+rejected.size()==6300&&selected.size()==calibration.GetEntries(),"Calibration inventory is incomplete");
  auto meta=metadataInfo(config,command);meta["condition"]=c;meta["selected_channels"]=selected.size();meta["rejected_channels"]=rejected.size();meta["cellIDs_encoding"]="layer*1000000+chip*10000+channel (legacy calibration)";meta["diagnostic_cellid_encoding"]="layer*100000+chip*10000+channel";metadata(*file,meta);closeRoot(*file);
  auto item=c;item["temperature_lower_C"]=interval["lower_C"];item["temperature_upper_C"]=interval["upper_C"];item["root"]=relativePath.string();item["selected_channels"]=selected.size();item["rejected_channels"]=rejected.size();item["status"]=selected.size()?"has_selected_channels":"no_reliable_gain";manifest.push_back(item);outputs.push_back(relativePath.string());totalSelected+=selected.size();totalRejected+=rejected.size();
 }
 writeTable(p.out/"conditions.tsv",manifest,{"condition","beam","period","temperature_lower_C","temperature_upper_C","temperature_bin","hit_tag","gain_tag","memory_cell","settings","selected_channels","rejected_channels","status","root"});
 Json summary={{"conditions",manifest.size()},{"selected_channel_conditions",totalSelected},{"rejected_channel_conditions",totalRejected},{"physical_channels_per_condition",6300},{"layers",layerChecks},{"command",commandLine(command)}};
 atomicJson(p.out/"summary.json",summary);seal(p.out,"MERGE_COMPLETE",signature,outputs,summary);
 std::cout<<"Completed merge: "<<manifest.size()<<" condition ROOT files, "<<totalSelected<<" selected channel-conditions"<<std::endl;
}
}
