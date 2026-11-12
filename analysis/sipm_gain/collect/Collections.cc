#include "Pipeline.hh"
#include <TH1D.h>
#include <TNamed.h>
#include <ctime>
// Existing raw reader is also used for new date/DAC collections.
void collectOne(const gain::Json&,const gain::Options&,const gain::Context&,const gain::fs::path&);
namespace gain {
namespace {
const std::vector<std::string> IntFields={"sid","run","cellid","hit_tag","gain_tag","temperature_bin","memory_cell","time_bin"};
const std::vector<std::string> CountFields={"entries","valid_temperature_entries","overflow"};
const std::vector<std::string> TempFields={"temperature_C","temperature_RMS_C"};
struct Input {Json row,record;fs::path directory;std::string hash;};
void mergeRoot(const std::vector<Input>& sources,const fs::path& output,const Json& provenance){
    auto out=openRoot(output,"RECREATE");
    TTree meta("Spectra","LED date/DAC collection; disjoint illuminated channels; unique spectrum IDs");
    std::map<std::string,int> ints;std::map<std::string,ULong64_t> counts;std::map<std::string,double> temperatures;
    for(auto& k:IntFields){ints[k]=0;meta.Branch(k.c_str(),&ints[k]);}
    for(auto& k:CountFields){counts[k]=0;meta.Branch(k.c_str(),&counts[k]);}
    for(auto& k:TempFields){temperatures[k]=NaN;meta.Branch(k.c_str(),&temperatures[k]);}
    int sourceIndex=0,sourceSid=0,ledGroup=0;
    meta.Branch("source_index",&sourceIndex);meta.Branch("source_sid",&sourceSid);meta.Branch("led_group",&ledGroup);
    TTree diagnostics("MemoryDiagnostics","HitTag0 conditional moments; source group retained");
    std::map<std::string,int> di;std::map<std::string,double> dd;
    for(auto k:{"run","cellid","gain_tag","memory_cell","temperature_bin"}){di[k]=0;diagnostics.Branch(k,&di[k]);}
    for(auto k:{"entries","mean_ADC","RMS_ADC"}){dd[k]=0;diagnostics.Branch(k,&dd[k]);}
    diagnostics.Branch("source_index",&sourceIndex);diagnostics.Branch("led_group",&ledGroup);
    Json sampling=Json::array();std::set<std::string> seen;int nextSid=0;
    for(const auto& source:sources){
        auto in=openRoot(source.directory/"spectra.root");auto* t=tree(*in,"Spectra");ledGroup=int(num(source.row,"led_group"));
        require(ledGroup>=0&&ledGroup<14,"Invalid source LED group");std::set<int> sourceIds;
        for(Long64_t e=0;e<t->GetEntries();++e){
            require(t->GetEntry(e)>0,"Cannot read source spectrum metadata");
            std::string key;
            for(auto& k:IntFields){ints[k]=int(leaf(*t,k.c_str()));if(k!="sid")key+=std::to_string(ints[k])+"/";}
            require(ints["memory_cell"]==0,"Source collection contains other memory cells");
            int cell=ints["cellid"];require(cell>=0&&cell/100000<30&&cell%10000/100==0&&ledChannel(ledGroup,cell%100000/10000*36+cell%100),"Source spectrum is outside illuminated LED group");
            require(seen.insert(key).second,"Duplicate spectrum in grouped collection");
            sourceSid=ints["sid"];require(sourceIds.insert(sourceSid).second,"Duplicate source sid");
            for(auto& k:CountFields){auto* p=t->GetLeaf(k.c_str());require(p,"Missing count field");counts[k]=p->GetValueLong64();}
            for(auto& k:TempFields)temperatures[k]=leaf(*t,k.c_str());
            std::unique_ptr<TH1D> h(in->Get<TH1D>(("s"+std::to_string(sourceSid)).c_str()));require(h&&h->GetNbinsX()==4096,"Missing source ADC histogram");h->SetDirectory(nullptr);
            ints["sid"]=nextSid++;h->SetName(("s"+std::to_string(ints["sid"])).c_str());out->cd();require(h->Write()>0,"Cannot write merged spectrum");meta.Fill();
        }
        auto* md=tree(*in,"MemoryDiagnostics");
        for(Long64_t e=0;e<md->GetEntries();++e){require(md->GetEntry(e)>0,"Cannot read memory diagnostics");for(auto& [k,v]:di)v=leaf(*md,k.c_str());for(auto& [k,v]:dd)v=leaf(*md,k.c_str());
            require(di["memory_cell"]==0&&ledChannel(ledGroup,di["cellid"]%100000/10000*36+di["cellid"]%100),"Invalid source memory diagnostic");diagnostics.Fill();}
        auto* tag=in->Get<TNamed>("sampling");sampling.push_back({{"source_index",sourceIndex},{"led_group",ledGroup},{"sample",source.row.at("sample")},{"raw",source.row.at("raw")},{"sampling",tag?tag->GetTitle():""}});
        ++sourceIndex;
    }
    out->cd();meta.Write();diagnostics.Write();TNamed("bootstrap_storage","disabled; memory cell 0 only").Write();
    TObjString(provenance.dump().c_str()).Write("source_collections_json");TObjString(sampling.dump().c_str()).Write("sources_json");
    TNamed("sampling","LED groups combined by beam/date/DAC; see sources_json; DAC is an internal setting, not volts").Write();out->Close();
}
void assemble(const Json& unit,const Context& c,const fs::path& out,const std::vector<fs::path>& paths,bool rebuild,bool temporarySources){
    require(paths.size()==unit.at("members").size(),"Wrong source collection count");
    std::vector<Input> inputs;std::vector<std::unique_ptr<Lock>> locks;Json provenance=Json::array(),cfg;
    for(size_t i=0;i<paths.size();++i){
        const auto& path=paths[i];require(fs::exists(path/"COLLECT_COMPLETE"),"Collection incomplete: "+path.string());locks.emplace_back(std::make_unique<Lock>(path/".lock",false,true));
        auto record=readJson(path/"inputs/collection.json");require(record.at("manifest_row")==unit.at("members").at(i),"Source collection manifest changed");
        const auto& selection=record.at("collection");require(num(selection,"memory_cell_selection",-1)==0&&!selection.value("bootstrap_storage",true),"Source must be memory-cell-0 collection without bootstrap storage");
        if(i==0)cfg=selection;else require(cfg==selection,"Mixed collection settings in one date/DAC unit");
        auto hash=trim(readText(path/"COLLECT_COMPLETE"));require(hash==sha256(path/"spectra.root"),"Source collection checksum mismatch");
        inputs.push_back({unit.at("members").at(i),record,path,hash});
        provenance.push_back({{"directory",path.string()},{"temporary_source",temporarySources},{"spectra_sha256",hash},{"record",record}});
    }
    Json expected={{"merge_version",1},{"manifest_row",unit},{"collection",cfg},{"source_collections",provenance},{"source",c.sources()}};
    auto done=out/"COLLECT_COMPLETE",recordPath=out/"inputs/collection.json";
    if(fs::exists(done)&&!rebuild){auto saved=readJson(recordPath);saved.erase("created_unix");require(saved==expected,"Grouped collection settings changed; use a new label");require(trim(readText(done))==sha256(out/"spectra.root"),"Grouped collection checksum changed");std::cout<<"Reused "<<out<<std::endl;return;}
    fs::create_directories(out/"inputs");fs::remove(done);
    auto scratch=c.project/".scratch/sipm_gain_merge";fs::create_directories(scratch);auto temp=scratch/(std::to_string(getpid())+"_"+str(unit,"sample")+".root");require(!fs::exists(temp),"Unexpected merge temporary file");
    try{
        if(str(unit,"kind")=="led")mergeRoot(inputs,temp,provenance);
        else{require(inputs.size()==1,"Non-LED units must have one input");fs::copy_file(inputs.front().directory/"spectra.root",temp);}
        fs::rename(temp,out/"spectra.root");
    }catch(...){fs::remove(temp);throw;}
    expected["created_unix"]=std::time(nullptr);writeJson(recordPath,expected);writeText(out/"inputs/collection_command.sh",commandLine(c.command)+'\n');writeText(done,sha256(out/"spectra.root")+'\n');
    std::cout<<"Collected unit "<<str(unit,"sample")<<" from "<<inputs.size()<<" inputs: "<<out<<std::endl;
}
}
void collectUnit(const Json& unit,const Options& o,const Context& c,const fs::path& out){
    // A complete collection is immutable unless explicitly recollected.
    if(fs::exists(out/"COLLECT_COMPLETE")&&!o.has("recollect")){
        auto saved=readJson(out/"inputs/collection.json");require(saved.at("manifest_row")==unit,"Collection unit changed");
        auto cfg=saved.at("collection");require(str(cfg,"channels")==o.get("channels","all")&&num(cfg,"events_per_file")==o.integer("events-per-file",0)&&num(cfg,"temperature_bin")==o.real("temperature-bin",.5)&&num(cfg,"time_bin")==o.real("time-bin",0),"Collection selection changed");
        for(const auto& source:saved.at("source_collections"))require(fingerprint(str(source.at("record").at("manifest_row"),"raw"))==source.at("record").at("input"),"Raw input changed since collection");
        require(trim(readText(out/"COLLECT_COMPLETE"))==sha256(out/"spectra.root"),"Collection checksum mismatch");return;
    }
    auto work=c.project/".scratch/sipm_gain_collect"/token(o.required("label"))/str(unit,"sample");fs::create_directories(work);
    std::vector<fs::path> paths;
    for(const auto& row:unit.at("members")){auto path=work/str(row,"sample");fs::create_directories(path);collectOne(row,o,c,path);paths.push_back(path);}
    assemble(unit,c,out,paths,o.has("recollect"),true);
    // All input records and source sampling metadata are retained in the final collection.
    fs::remove_all(work);
}
void regroup(const Options& o,const Context& c){
    auto from=token(o.required("from-label")),label=token(o.required("label"));require(from!=label,"Regroup into a new label to preserve existing results");
    auto source=c.result/from;require(!fs::exists(source/"inputs/analysis_units.json"),"Source is already grouped; use its fit stage directly");
    auto rows=readTable(source/"inputs/manifest.tsv"),units=analysisUnits(rows);
    auto allUnits=units;
    if(o.has("row")){auto i=o.integer("row",-1);require(i>=0&&i<static_cast<long long>(units.size()),"Analysis unit index out of range");units=Json::array({units.at(i)});}
    if(o.has("dry-run")){std::cout<<Json({{"raw_inputs",rows.size()},{"analysis_units",allUnits.size()},{"selected_units",units.size()},{"from_label",from},{"label",label}}).dump(2)<<'\n';return;}
    freezeInputs(c,label,rows);
    Json origin={{"from_label",from},{"source_manifest_sha256",sha256(source/"inputs/manifest.tsv")},{"source_code_manifest",readJson(source/"inputs/source_manifest.json")}};
    auto originPath=c.result/label/"inputs/regroup_source.json";
    if(fs::exists(originPath))require(readJson(originPath)==origin,"Regroup source changed");else writeJson(originPath,origin);
    for(const auto& unit:units){auto out=samplePath(c,label,unit)/"collection";fs::create_directories(out);Lock lock(out/".lock");std::vector<fs::path> paths;
        for(const auto& row:unit.at("members"))paths.push_back(samplePath(c,from,row)/"collection");
        assemble(unit,c,out,paths,false,false);
    }
    std::cout<<"Regroup complete. No raw events reread and no fits submitted. Output: "<<c.result/label<<std::endl;
}
}
