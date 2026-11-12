#include "Common.hh"

namespace selected_gain {
double temperatureWidth(const Paths& p,const Json& c) {
 static std::map<std::string,double> cache;
 const auto identity=p.input.string()+"/"+str(c,"beam")+"/"+str(c,"period");
 if(cache.count(identity))return cache.at(identity);
 const auto comparison=readJson(p.input/"inputs/configuration.json");
 const auto& frozen=comparison.at("input_snapshot").at("collections");
 const auto prefix=str(c,"beam")+"_led_"+str(c,"period")+"_dac";double width=NaN;
 for(auto it=frozen.begin();it!=frozen.end();++it){
  if(it.key().rfind(prefix,0)!=0)continue;
  const auto directory=p.project/"analysis/result/sipm_gain"/p.label/"data/led"/token(it.key())/"collection";
  Lock lock(directory/".lock",false,true);
  require(it.value()==sha256(directory/"COLLECT_COMPLETE"),"Collection differs from compared input: "+it.key());
  auto record=readJson(directory/"inputs/collection.json");auto row=record.at("manifest_row");
  require(str(row,"kind")=="led"&&str(row,"beam")==str(c,"beam")&&str(row,"period")==str(c,"period"),"Temperature provenance mismatch");
  double w=num(record.at("collection"),"temperature_bin");require(std::isfinite(w)&&w>0,"Invalid collection temperature width");
  require(!std::isfinite(width)||width==w,"Mixed temperature widths within period");width=w;
 }
 require(std::isfinite(width),"No collection temperature provenance: "+identity);cache[identity]=width;return width;
}
namespace {
void atomicText(const fs::path& path,const std::string& text) {
 auto tmp=fs::path(path.string()+".tmp-"+std::to_string(getpid()));require(!fs::exists(tmp),"Temporary output exists: "+tmp.string());
 writeText(tmp,text);fs::rename(tmp,path);
}
std::string tableText(const Json& rows,const std::vector<std::string>& columns) {
 std::ostringstream s;for(size_t i=0;i<columns.size();++i)s<<(i?"\t":"")<<columns[i];s<<'\n';
 for(const auto& row:rows){for(size_t i=0;i<columns.size();++i)s<<(i?"\t":"")<<(row.contains(columns[i])?tsvValue(row[columns[i]]):"");s<<'\n';}return s.str();
}
}
void migrateTemperaturePaths(const Paths& p,const Json& config,const std::vector<std::string>& command) {
 (void)config;const auto signature=configHash(p);auto journal=p.out/"inputs/temperature_path_migration.json";Json plan;
 if(fs::exists(journal)&&readJson(journal).value("status","")=="pending"){
  plan=readJson(journal);require(plan.at("configuration_sha256")==signature,"Migration configuration changed");
 }else{
  plan={{"status","pending"},{"configuration_sha256",signature},{"command",commandLine(command)},{"moves",Json::object()},{"writes",Json::object()},{"markers",Json::array()}};
  std::vector<fs::path> markers;if(fs::exists(p.out/"MERGE_COMPLETE"))markers.push_back(p.out/"MERGE_COMPLETE");
  for(int l=0;l<30;++l){auto m=p.out/"plots/layers"/layerName(l)/"PLOTS_COMPLETE";if(fs::exists(m))markers.push_back(m);}
  require(!markers.empty(),"No completed ROOT/plot outputs to rename");
  auto rewrite=[&](const fs::path& base,const std::string& name){
   std::vector<std::string> parts;for(const auto& x:fs::path(name))parts.push_back(x.string());fs::path from,to;
   for(size_t i=0;i<parts.size();++i){const auto& x=parts[i];from/=x;
    if(x.rfind("temperature_bin_",0)==0){require(i>=3&&parts[i-3]=="data","Unexpected temperature path: "+name);size_t used=0;auto number=x.substr(16);double bin=std::stod(number,&used);require(used==number.size(),"Bad temperature index in path");
     Json c={{"beam",parts[i-2]},{"period",parts[i-1]}};to/=temperatureDirectory(temperatureInterval(bin,temperatureWidth(p,c)));
     auto old=(base/from).lexically_normal().lexically_relative(p.out).string(),dest=(base/to).lexically_normal().lexically_relative(p.out).string();
     require(old.rfind("../",0)!=0&&dest.rfind("../",0)!=0,"Output escapes selection directory");
     if(plan["moves"].contains(old))require(plan["moves"][old]==dest,"Conflicting temperature migration");else plan["moves"][old]=dest;
    }else to/=x;
   }return to.string();
  };
  auto addWrite=[&](const fs::path& path,const std::string& content,bool marker){auto old=readText(path);if(old==content)return;auto rel=path.lexically_relative(p.out).string();plan["writes"][rel]={{"before_sha256",digest(old)},{"after_sha256",digest(content)},{"text",content},{"marker",marker}};};
  for(const auto& marker:markers){const auto base=marker.parent_path();auto saved=verifyComplete(base,marker.filename().string().c_str(),signature);auto changed=saved;Json files=Json::object();
   const bool merged=marker.filename()=="MERGE_COMPLETE";const auto table=base/(merged?"conditions.tsv":"plot_index.tsv");
   auto rows=readTable(table);std::istringstream header(readText(table));std::string line;std::getline(header,line);auto columns=tsvLine(line);
   for(auto& row:rows){const auto field=merged?"root":"figure";row[field]=rewrite(p.out,str(row,field));
    if(merged){auto interval=temperatureInterval(num(row,"temperature_bin"),temperatureWidth(p,row));row["temperature_lower_C"]=interval["lower_C"];row["temperature_upper_C"]=interval["upper_C"];}
   }
   if(merged)for(const std::string key:{"temperature_lower_C","temperature_upper_C"})if(std::find(columns.begin(),columns.end(),key)==columns.end())columns.push_back(key);
   addWrite(table,tableText(rows,columns),false);
   for(auto it=saved.at("files").begin();it!=saved.at("files").end();++it){auto name=rewrite(base,it.key());require(!files.contains(name),"Duplicate migrated file");auto rel=(base/it.key()).lexically_normal().lexically_relative(p.out).string();files[name]=plan["writes"].contains(rel)?plan["writes"][rel]["after_sha256"]:it.value();}
   changed["files"]=files;addWrite(marker,changed.dump(2)+'\n',true);plan["markers"].push_back(marker.lexically_relative(p.out).string());
  }
  if(plan["moves"].empty()&&plan["writes"].empty()){std::cout<<"Temperature directories already use Celsius ranges"<<std::endl;return;}
  std::set<std::string> destinations;
  for(auto it=plan["moves"].begin();it!=plan["moves"].end();++it){auto dest=it.value().get<std::string>();require(destinations.insert(dest).second,"Temperature directory collision");require(fs::is_directory(p.out/it.key())&&!fs::exists(p.out/dest),"Migration source missing or destination already exists: "+it.key());}
  atomicJson(journal,plan);
 }
 for(auto it=plan["moves"].begin();it!=plan["moves"].end();++it){auto src=p.out/it.key(),dest=p.out/it.value().get<std::string>();if(fs::exists(src)){require(!fs::exists(dest),"Refusing to overwrite migration destination");fs::rename(src,dest);}else require(fs::is_directory(dest),"Both migration directories missing");}
 for(bool marker:{false,true})for(auto it=plan["writes"].begin();it!=plan["writes"].end();++it){auto& w=it.value();if(w.at("marker")!=marker)continue;auto path=p.out/it.key();const auto hash=sha256(path);require(hash==w.at("before_sha256")||hash==w.at("after_sha256"),"Output changed during migration: "+it.key());if(hash!=w.at("after_sha256"))atomicText(path,str(w,"text"));}
 for(auto& name:plan["markers"]){auto marker=p.out/name.get<std::string>();verifyComplete(marker.parent_path(),marker.filename().string().c_str(),signature);}
 plan["status"]="complete";plan["validation"]="All ROOT/PNG checksums unchanged; tables and completion manifests updated";atomicJson(journal,plan);
 std::cout<<"Renamed "<<plan["moves"].size()<<" temperature directories; ROOT/PNG contents preserved"<<std::endl;
}
}
