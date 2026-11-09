#pragma once
#include "GainIO.hh"
#include "LEDGroups.hh"
namespace gain {
inline const std::vector<std::string> StringFields={"sample","kind","beam","energy_GeV","period","dac","settings","method","spectrum_status","status","reasons","uncertainty_status","condition_note","source_sample","source_raw"};
inline const std::vector<std::string> NumericFields={"sid","run","cellid","hit_tag","gain_tag","temperature_bin","memory_cell","time_bin","entries","valid_temperature_entries","overflow","temperature_C","temperature_RMS_C","gain","stat_error","fit_error","stability_ADC","null_p","fit_status","covariance_status","deviance","ndf","alternative_delta_bic","bootstrap_requested","bootstrap_success","variants_success","variants_total","source_sid","led_group"};
inline fs::path samplePath(const Context& c,const std::string& label,const Json& row){return c.result/token(label)/"data"/str(row,"kind")/str(row,"sample");}
inline std::string partName(int group,int groups){std::ostringstream s;s<<"part_"<<std::setfill('0')<<std::setw(4)<<group<<"_of_"<<std::setw(4)<<groups;return s.str();}
inline int channelGroup(int cell,int groups){return ((cell/100000)*210+(cell%100000/10000)*36+cell%100)%groups;}
inline std::vector<std::string> selectedMethods(const std::string& name){if(name=="all")return {Methods.begin(),Methods.end()};require(std::find(Methods.begin(),Methods.end(),name)!=Methods.end(),"Unknown --method "+name);return {name};}
// The raw-file manifest remains unchanged; analysis units aggregate disjoint LED groups.
inline Json analysisUnits(const Json& rows){
    Json units=Json::array();std::map<std::string,size_t> indices;std::map<std::string,std::set<int>> seen;
    for(const auto& row:rows){
        auto kind=str(row,"kind");std::string name=str(row,"sample");
        if(kind=="led"){
            auto dac=str(row,"dac");require(std::regex_match(dac,std::regex("[0-9]+")),"LED DAC must be an internal integer setting");
            name=token(str(row,"beam")+"_led_"+str(row,"period")+"_dac"+dac);
            require(seen[name].insert(int(num(row,"led_group"))).second,"Duplicate LED group in date/DAC unit: "+name);
        }
        auto key=kind+"/"+name;
        if(!indices.count(key)){Json unit=row;unit["sample"]=name;unit["members"]=Json::array();if(kind=="led"){unit["led_group"]="-1";unit["raw"]="";}indices[key]=units.size();units.push_back(unit);}
        auto& unit=units.at(indices.at(key));
        for(auto field:{"kind","beam","period","dac","settings","energy_GeV"})require(str(unit,field)==str(row,field),"Mixed conditions in analysis unit: "+name);
        unit["members"].push_back(row);
    }
    for(auto& unit:units)if(str(unit,"kind")=="led")std::sort(unit["members"].begin(),unit["members"].end(),[](const Json& a,const Json& b){return num(a,"led_group")<num(b,"led_group");});
    return units;
}
inline void writeGainRoot(const fs::path& path,const Json& rows,const Json& metadata){
    auto file=openRoot(path,"RECREATE");TTree t("Gain","Per-method SPE gain; stat_error is NaN when not evaluated");
    std::map<std::string,std::string> strings;std::map<std::string,double> numbers;
    for(auto& name:StringFields){strings[name]="";t.Branch(name.c_str(),&strings[name]);}
    for(auto& name:NumericFields){numbers[name]=NaN;t.Branch(name.c_str(),&numbers[name]);}
    for(auto& row:rows){for(auto& name:StringFields)strings[name]=str(row,name);for(auto& name:NumericFields)numbers[name]=num(row,name);t.Fill();}
    t.Write();TObjString info(metadata.dump().c_str());info.Write("metadata_json");file->Close();
}
inline Json readGainRoot(const fs::path& path){
    auto file=openRoot(path);auto* t=tree(*file,"Gain");std::map<std::string,std::string*> strings;
    for(auto& name:StringFields){strings[name]=nullptr;if(t->GetBranch(name.c_str()))require(t->SetBranchAddress(name.c_str(),&strings[name])>=0,"Invalid string branch "+name);else require(name=="source_sample"||name=="source_raw","Missing string branch "+name);}
    Json rows=Json::array();for(Long64_t i=0;i<t->GetEntries();++i){require(t->GetEntry(i)>0,"Cannot read gain record");Json row;for(auto& name:StringFields)row[name]=strings[name]?*strings[name]:"";for(auto& name:NumericFields)row[name]=t->GetLeaf(name.c_str())?leaf(*t,name.c_str()):NaN;rows.push_back(row);}
    t->ResetBranchAddresses();for(auto& item:strings)delete item.second;file->Close();return rows;
}
inline void freezeInputs(const Context& c,const std::string& label,const Json& rows,bool refresh=false){
    auto inputs=c.result/token(label)/"inputs";fs::create_directories(inputs);Lock lock(inputs/".lock",true);
    auto units=analysisUnits(rows);auto unitPath=inputs/"analysis_units.json";
    require(fs::exists(unitPath)||!fs::exists(c.result/token(label)/"data"),"Old per-file layout: use regroup --from-label OLD --label NEW");
    if(fs::exists(unitPath))require(readJson(unitPath)==units,"Analysis units changed; use a new label");else writeJson(unitPath,units);
    auto path=inputs/"manifest.tsv";if(fs::exists(path))require(readTable(path)==rows,"Label has a different manifest; use a new label");else writeTable(path,rows,ManifestFields);
    auto sources=c.sources();for(auto it=sources.begin();it!=sources.end();++it){fs::path rel=it.key(),src;
        std::string prefix="analysis/sipm_gain/";
        if(it.key().rfind(prefix,0)==0)src=c.here/it.key().substr(prefix.size());else src=c.project/rel;
        auto dest=inputs/"source"/rel;
        if(fs::exists(dest)&&!refresh)require(sha256(dest)==it.value(),"Source changed for label; use --refit or a new label");
        else{fs::create_directories(dest.parent_path());fs::copy_file(src,dest,fs::copy_options::overwrite_existing);}
    }writeJson(inputs/"source_manifest.json",sources);
}
inline void sealPart(const fs::path& out,const Json& config){
    Json checks;for(auto name:{"gain.root","gain_estimates.tsv","fits.json","diagnostics.root"})checks[name]=sha256(out/name);
    writeJson(out/"FIT_COMPLETE",{{"files",checks},{"config",config}});
}
inline Json verifyPart(const fs::path& out){auto marker=readJson(out/"FIT_COMPLETE");for(auto it=marker.at("files").begin();it!=marker.at("files").end();++it)require(sha256(out/it.key())==it.value(),"Changed result: "+(out/it.key()).string());return marker;}
void regroup(const Options&,const Context&);
void collectUnit(const Json&,const Options&,const Context&,const fs::path&);
void submit(const Options&,const Context&);
void task(const Options&,const Context&);
void run(const Options&,const Context&);
}
