#include "Pipeline.hh"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <ctime>

using namespace gain;
const std::string DefaultChannels="all";
const std::vector<std::string> EstimateFields={"sample","kind","beam","settings","period","dac","energy_GeV","sid","run","cellid","hit_tag","gain_tag","temperature_bin","memory_cell","time_bin","entries","temperature_C","temperature_RMS_C","spectrum_status","method","gain","stat_error","fit_error","q16","q84","stability_ADC","null_p","status","reasons","fit_status","covariance_status","deviance","ndf","alternative_delta_bic","bootstrap_success","bootstrap_requested","variants_success","variants_total"};
const std::vector<std::string> PairFields={"sample","sid","cellid","method_a","method_b","delta_ADC","paired_stat_error_ADC","covariance_ADC2","paired_replicates","pull","both_accepted","stability_allowance_ADC","disagreement"};
const std::vector<std::string> MemoryFields={"run","cellid","gain_tag","memory_cell","temperature_bin","entries","mean_ADC","RMS_ADC"};
Options parse(int argc,char** argv){
    Options o;o.action=argc>1?argv[1]:"help";
    const std::set<std::string> switches={"help","refit","recollect","dry-run","allow-partial","execute"};
    for(int i=2;i<argc;++i){std::string k=argv[i];require(k.rfind("--",0)==0,"Expected option, got "+k);k=k.substr(2);require(!o.has(k),"Duplicate --"+k);if(switches.count(k))o.flags.insert(k);else{require(i+1<argc,"Missing value for --"+k);o.values[k]=argv[++i];}}
    std::map<std::string,std::set<std::string>> allowed={
        {"discover",{"output","data-root","profile","beams","kinds","periods","energies"}},
        {"run",{"manifest","label","row","stage","channels","events-per-file","temperature-bin","time-bin","method","channel-group","channel-groups","hit-tags","plot-limit","refit","recollect","dry-run","gain-min","gain-max","adc-min","adc-max","max-span","min-entries","min-cycles","null-replicates","bootstrap","alpha","stability-fraction","seed"}},
        {"task",{"tasks","index"}},
        {"regroup",{"from-label","label","row","dry-run"}},
        {"compare",{"label","manifest","allow-partial","reference","reference-format","temperature-tolerance","reference-temperature","stage","layer","comparison-name"}},
        {"submit",{"manifest","label","channels","events-per-file","bootstrap","min-entries","partition","max-concurrent","memory","walltime","execute","stage","channel-groups","hit-tags","plot-limit","temperature-bin","null-replicates","max-span","gain-min","gain-max","adc-min","adc-max"}}
    };
    if(o.action=="help"||o.action=="--help"||o.action=="-h"||o.has("help"))return o;
    require(allowed.count(o.action),"Unknown command: "+o.action);
    for(auto& [k,v]:o.values)require(allowed[o.action].count(k),"Unknown --"+k+" for "+o.action);
    for(auto& k:o.flags)require(allowed[o.action].count(k),"Unknown --"+k+" for "+o.action);
    return o;
}
Config configuration(const Options& o){Config c;
#define REAL(field) {std::string key=#field;std::replace(key.begin(),key.end(),'_','-');c.field=o.real(key,c.field);}
#define INTEGER(field) {std::string key=#field;std::replace(key.begin(),key.end(),'_','-');auto v=o.integer(key,c.field);require(v>=0&&v<=2147483647LL,"Invalid --"+key);c.field=v;}
    REAL(gain_min);REAL(gain_max);REAL(adc_min);REAL(adc_max);REAL(min_cycles);REAL(alpha);REAL(stability_fraction);
    INTEGER(max_span);INTEGER(min_entries);INTEGER(null_replicates);INTEGER(bootstrap);INTEGER(seed);
#undef REAL
#undef INTEGER
    c.validate();return c;
}
void channelsValid(const std::string& selection){
    if(selection=="all")return;
    require(std::regex_match(selection,std::regex("[0-9]+(,[0-9]+)*")),"Invalid channel list");std::set<int> seen;
    for(auto& text:split(selection)){long long c=std::stoll(text);require(c>=0&&c/100000<30&&c%100000/10000<6&&c%100<36&&c%10000/100==0&&!(c%100000/10000==5&&c%100>=30),"Use connected physical CellID, layer*100000+chip*10000+channel");require(seen.insert(c).second,"Duplicate channel");}
}
std::vector<fs::path> children(const fs::path& root,bool directories){require(fs::is_directory(root),"Missing input directory "+root.string());std::vector<fs::path> p;for(auto& e:fs::directory_iterator(root))if(directories?e.is_directory():(e.is_regular_file()&&e.path().extension()==".root"))p.push_back(e.path());std::sort(p.begin(),p.end());return p;}
bool selected(const std::string& wanted,const std::string& value){if(wanted.empty())return true;auto list=split(wanted);return std::find(list.begin(),list.end(),value)!=list.end();}
void discover(const Options& o){
    fs::path output=o.required("output");require(!fs::exists(output),"Manifest already exists: "+output.string());
    auto profile=o.get("profile","pilot");require(profile=="pilot"||profile=="full","Invalid profile");Json rows=Json::array();
    fs::path data=o.get("data-root","/megraid01/users/data_beamtest/ECAL_data/analysed/2023");
    for(auto& beam:split(o.get("beams","sps,ps"))){require(beam=="ps"||beam=="sps","Invalid beam");for(auto& kind:split(o.get("kinds","led,electron,muon"))){
        require(kind=="led"||kind=="electron"||kind=="muon","Unknown kind");auto root=data/beam/"decode"/(kind=="led"?(beam=="sps"?"LEDCali":"lightCalibData"):kind=="electron"?"e-":"mu-");
        auto dirs=children(root,true);
        if(kind=="led"){
            if(!o.has("periods")&&profile=="pilot"){require(!dirs.empty(),"No LED periods");dirs={beam=="sps"?dirs.back():dirs.front()};}
            for(auto& period:dirs){if(!selected(o.get("periods"),period.filename().string()))continue;
                for(auto& file:children(period,false)){std::smatch match;auto name=file.filename().string();require(std::regex_match(name,match,std::regex("group([0-9]+)_a?dac([0-9]+)\\.root")),"Unrecognized LED filename: "+name);
                    int group=std::stoi(match[1]),dac=std::stoi(match[2]);if(profile=="pilot"&&((group!=6&&group!=7)||(dac!=2150&&dac!=2300)))continue;
                    rows.push_back({{"sample",beam+"_led_"+period.filename().string()+"_group"+std::to_string(group)+"_dac"+std::to_string(dac)},
                        {"kind",kind},{"beam",beam},{"energy_GeV",""},{"period",period.filename().string()},{"led_group",std::to_string(group)},
                        {"dac",std::to_string(dac)},{"settings","unconfirmed"},{"raw",fs::canonical(file).string()}});
                }
            }
        }else for(auto& energy:dirs){auto ename=energy.filename().string(),e=ename;if(e.size()>=3&&e.substr(e.size()-3)=="GeV")e.resize(e.size()-3);
            if(o.has("energies")){if(!selected(o.get("energies"),e))continue;}else if(profile=="pilot"&&ename!=(beam=="sps"?"100GeV":kind=="muon"?"10GeV":"5GeV"))continue;
            auto files=children(energy,false);if(profile=="pilot"&&!files.empty()){std::string preferred=kind=="muon"?(beam=="sps"?"Run44_":"Run11_"):"";auto chosen=files.front();if(!preferred.empty())for(auto& f:files)if(f.filename().string().find(preferred)!=std::string::npos){chosen=f;break;}files={chosen};}
            for(auto& file:files)rows.push_back({{"sample",beam+"_"+kind+"_"+ename+"_"+file.stem().string()},{"kind",kind},{"beam",beam},{"energy_GeV",e},
                {"period",file.stem().string()},{"led_group","-1"},{"dac",""},{"settings","unconfirmed"},{"raw",fs::canonical(file).string()}});
        }
    }}require(!rows.empty(),"No matching inputs");writeTable(output,rows,ManifestFields);std::cout<<rows.size()<<" inputs: "<<output<<'\n';
}
Json collectionConfig(const Options& o){
    auto channels=o.get("channels",DefaultChannels);channelsValid(channels);
    auto events=o.integer("events-per-file",0);double temperature=o.real("temperature-bin",.5),time=o.real("time-bin",0.);
    require(events>=0&&temperature>0&&time>=0,"Invalid collection configuration");
    return {{"channels",channels},{"events_per_file",events},{"temperature_bin",temperature},{"time_bin",time},{"memory_cell_selection",0},{"bootstrap_storage",false}};
}
void collectOne(const Json& row,const Options& o,const Context& context,const fs::path& out){
    auto cfg=collectionConfig(o);Json expected={{"input",fingerprint(str(row,"raw"))},{"manifest_row",row},{"collection",cfg},{"source",context.sources(true)}};
    auto done=out/"COLLECT_COMPLETE",record=out/"inputs/collection.json";
    if(fs::exists(done)&&!o.has("recollect")){auto saved=readJson(record);saved.erase("created_unix");require(saved==expected,"Cached collection differs; choose another label or original configuration");require(sha256(out/"spectra.root")==trim(readText(done)),"Collection checksum changed");return;}
    require(!fs::exists(out/"spectra.root")||o.has("recollect"),"Incomplete existing collection: "+out.string());fs::create_directories(out/"inputs");
    auto scratch=context.project/".scratch/sipm_gain";fs::create_directories(scratch);auto temp=scratch/(std::to_string(getpid())+"_"+str(row,"sample")+".root");require(!fs::exists(temp),"Unexpected temporary output");
    std::vector<std::string> cmd={(context.here/"Collect").string(),str(row,"raw"),temp.string(),str(cfg,"channels"),str(cfg,"events_per_file"),str(cfg,"temperature_bin"),str(row,"led_group"),str(cfg,"time_bin")};
    writeText(out/"inputs/collect_command.sh",commandLine(cmd)+'\n');execute(cmd,out/"collect.log");require(fingerprint(str(row,"raw"))==expected["input"],"Input changed while reading");
    fs::rename(temp,out/"spectra.root");fs::remove(out/"FIT_COMPLETE");expected["created_unix"]=std::time(nullptr);writeJson(record,expected);writeText(done,sha256(out/"spectra.root")+'\n');
}
void diagnostics(TFile& file,const std::string& stem,const Result& result,const fs::path& png){
    auto& d=result.detail;file.cd();auto dir=file.mkdir(stem.c_str());require(dir,"Cannot create diagnostic directory");dir->cd();
    TH1D h(("spectrum_"+stem).c_str(),";HG ADC;entries",d.x.size(),d.x.front()-.5,d.x.back()+.5);h.SetStats(false);
    for(size_t i=0;i<d.x.size();++i){h.SetBinContent(i+1,d.y[i]);h.SetBinError(i+1,std::sqrt(d.y[i]));}h.Write("spectrum");
    std::vector<TGraph> curves;curves.reserve(3);curves.push_back(graph(d.x,d.env,"env_"+stem));curves.back().Write("env");
    std::vector<std::string> curveNames={"env"};
    for(int i=0;i<2;++i){if(d.curves[i].empty())continue;curveNames.push_back(Methods[i]);curves.push_back(graph(d.x,d.curves[i],Methods[i]+"_"+stem));curves.back().Write(Methods[i].c_str());Vec residual;for(size_t k=0;k<d.x.size();++k)residual.push_back((d.y[k]-d.curves[i][k])/std::sqrt(std::max(d.curves[i][k],1.)));auto g=graph(d.x,residual,Methods[i]+"_residual_"+stem);g.Write((Methods[i]+"_residual").c_str());}
    auto fft=graph(d.freq,d.power,"fft_"+stem),acf=graph(d.lags,d.acf,"autocorrelation_"+stem);fft.Write("fft");acf.Write("autocorrelation");
    if(png.empty())return;
    TCanvas canvas(("c_"+stem).c_str(),stem.c_str(),1400,1000);canvas.Divide(2,2);canvas.cd(1);h.SetTitle((stem+";HG ADC;entries").c_str());h.Draw("E");TLegend legend(.60,.64,.89,.89);legend.AddEntry(&h,"data","lep");int colors[]={kGray+2,kBlue+1,kRed+1};
    for(size_t i=0;i<curves.size();++i){curves[i].SetLineColor(colors[i]);curves[i].SetLineWidth(2);curves[i].Draw("L SAME");legend.AddEntry(&curves[i],curveNames[i].c_str(),"l");}legend.Draw();
    canvas.cd(2);fft.SetTitle("FFT of whitened, tapered residual;frequency [1/ADC];power");fft.Draw("AL");canvas.cd(3);acf.SetTitle("Linear autocorrelation;lag [ADC];correlation");acf.Draw("AL");
    canvas.cd(4);TPaveText text(.06,.06,.94,.94,"NDC");text.SetFillColor(0);text.SetTextAlign(12);text.AddText(("Status: "+str(result.report,"status")).c_str());
    for(auto it=result.report["methods"].begin();it!=result.report["methods"].end();++it){const auto& r=it.value();std::ostringstream line;line<<it.key()<<": "<<std::fixed<<std::setprecision(3)<<num(r,"gain")<<" ADC/p.e.";text.AddText(line.str().c_str());text.AddText(("  "+str(r,"status")+" "+str(r,"reasons")).c_str());}
    text.AddText("Statistical uncertainty: not evaluated (bootstrap off)");
    text.Draw();canvas.SaveAs(png.c_str());canvas.Close();file.cd();
}
void fitPart(const Json& row,const Options& o,const Context& context,const fs::path& sample,const std::string& method,int group,int groups){
    auto collection=sample/"collection",out=sample/method/partName(group,groups);
    fs::create_directories(out);Lock partLock(out/".lock");Lock collectionLock(collection/".lock",false,true);
    auto collected=readJson(collection/"inputs/collection.json");require(collected["manifest_row"]==row,"Manifest changed since collection");
    require(num(collected["collection"],"memory_cell_selection",-1)==0,"Collection must explicitly select memory cell 0");
    require(fs::exists(collection/"COLLECT_COMPLETE"),"Collection incomplete");auto rootHash=trim(readText(collection/"COLLECT_COMPLETE"));
    require(sha256(collection/"spectra.root")==rootHash,"Collected spectrum checksum mismatch");
    Config cfg=configuration(o);int limit=o.integer("plot-limit",2);require(limit>=0,"Invalid plot limit");
    auto tags=split(o.get("hit-tags","1"));std::set<int> chosen;for(auto& tag:tags){require(tag=="0"||tag=="1","Invalid hit tag");chosen.insert(std::stoi(tag));}require(!chosen.empty(),"Empty hit tags");
    Json config={{"config",cfg.json()},{"method",method},{"group",group},{"groups",groups},{"hit_tags",o.get("hit-tags","1")},{"source",context.sources()},{"spectra_sha256",rootHash},{"plot_limit",limit},{"root_version",gROOT->GetVersion()}};
    auto done=out/"FIT_COMPLETE";
    if(fs::exists(done)&&!o.has("refit")){require(verifyPart(out)["config"]==config,"Part settings differ; use a new label or --refit");std::cout<<"Reused "<<out<<std::endl;return;}
    fs::remove(done);writeJson(out/"inputs/fit_configuration.json",config);writeText(out/"inputs/fit_command.sh",commandLine(context.command)+'\n');
    auto input=openRoot(collection/"spectra.root");auto* metadata=tree(*input,"Spectra");auto root=openRoot(out/"diagnostics.root","RECREATE");
    Json reports=Json::array(),flat=Json::array();int plotted=0;
    for(Long64_t entry=0;entry<metadata->GetEntries();++entry){
        require(metadata->GetEntry(entry)>0,"Cannot read metadata");
        int cell=leaf(*metadata,"cellid"),tag=leaf(*metadata,"hit_tag");
        if(!chosen.count(tag)||channelGroup(cell,groups)!=group)continue;
        require(leaf(*metadata,"memory_cell")==0,"Unexpected memory cell");
        Json r;for(auto key:{"sid","run","cellid","hit_tag","gain_tag","temperature_bin","memory_cell","time_bin","entries","valid_temperature_entries","overflow"})r[key]=static_cast<long long>(leaf(*metadata,key));
        r["source_sid"]=metadata->GetLeaf("source_sid")?leaf(*metadata,"source_sid"):leaf(*metadata,"sid");
        r["led_group"]=metadata->GetLeaf("led_group")?leaf(*metadata,"led_group"):num(row,"led_group");
        r["source_sample"]=metadata->GetLeaf("source_index")?str(row.at("members").at(int(leaf(*metadata,"source_index"))),"sample"):str(row,"sample");
        r["source_raw"]=metadata->GetLeaf("source_index")?str(row.at("members").at(int(leaf(*metadata,"source_index"))),"raw"):str(row,"raw");
        r["temperature_C"]=leaf(*metadata,"temperature_C");r["temperature_RMS_C"]=leaf(*metadata,"temperature_RMS_C");int sid=num(r,"sid");
        std::unique_ptr<TH1D> h(dynamic_cast<TH1D*>(input->Get(("s"+std::to_string(sid)).c_str())));require(h&&h->GetNbinsX()==4096,"Missing raw spectrum");h->SetDirectory(nullptr);
        Vec counts(4096);for(int k=0;k<4096;++k)counts[k]=h->GetBinContent(k+1);
        cfg.seed=configuration(o).seed+int(num(r,"source_sid"));auto result=analyze(counts,str(row,"kind"),cfg,method);
        result.report["spectrum"]=r;result.report["sample"]=row;reports.push_back(result.report);
        Json line=r;for(auto& key:{"sample","kind","beam","settings","period","dac","energy_GeV"})line[key]=row[key];
        line["method"]=method;line["spectrum_status"]="comparison_only";line["condition_note"]=num(r,"temperature_bin")==-1?"temperature_missing":"memory_cell_0";
        for(auto it=result.report["methods"][method].begin();it!=result.report["methods"][method].end();++it)if(!it.value().is_array()&&!it.value().is_object())line[it.key()]=it.value();
        flat.push_back(line);
        if(!result.detail.x.empty()){fs::path png;if(limit==0||plotted<limit){fs::create_directories(out/"figures");png=out/"figures"/("s"+std::to_string(sid)+".png");++plotted;}diagnostics(*root,"s"+std::to_string(sid),result,png);}
        std::cout<<str(row,"sample")<<' '<<method<<" sid="<<sid<<" cell="<<cell<<' '<<str(result.report,"status")<<std::endl;
    }
    root->Close();input->Close();writeJson(out/"fits.json",reports);
    auto fields=EstimateFields;fields.push_back("uncertainty_status");fields.push_back("condition_note");for(auto key:{"source_sid","led_group","source_sample","source_raw"})fields.push_back(key);writeTable(out/"gain_estimates.tsv",flat,fields);
    writeGainRoot(out/"gain.root",flat,config);sealPart(out,config);
}
namespace gain {
void run(const Options& o,const Context& context){
    auto rows=manifest(o.required("manifest"));auto label=token(o.required("label"));auto cfg=collectionConfig(o);auto fitConfig=configuration(o);
    auto stage=o.get("stage","collect");require(stage=="collect"||stage=="fit"||stage=="all","Invalid stage");
    int groups=o.integer("channel-groups",4);require(groups>0&&groups<=6300,"Invalid channel-groups");
    int first=0,last=groups;if(o.has("channel-group")){first=o.integer("channel-group",-1);require(first>=0&&first<groups,"Invalid channel-group");last=first+1;}
    auto methods=selectedMethods(o.get("method","all"));auto allRows=rows;rows=analysisUnits(rows);
    if(o.has("row")){auto i=o.integer("row",-1);require(i>=0&&i<static_cast<long long>(rows.size()),"Manifest row out of range");rows=Json::array({rows.at(i)});}
    if(o.has("dry-run")){std::cout<<Json({{"stage",stage},{"rows",rows.size()},{"methods",methods},{"groups",groups},{"selected_group_count",last-first},{"bootstrap",0},{"memory_cell",0},{"output",(context.result/label).string()}}).dump(2)<<'\n';return;}
    freezeInputs(context,label,allRows,o.has("refit"));
    if(stage!="collect"){
        auto plan=context.result/label/"inputs/analysis_plan.json";Lock lock(plan.parent_path()/".plan.lock",true);
        Json expected={{"channel_groups",groups},{"methods",Methods},{"fit_config",fitConfig.json()},{"hit_tags",o.get("hit-tags","1")},{"bootstrap",0}};
        if(fs::exists(plan)&&!o.has("refit"))require(readJson(plan)==expected,"Analysis plan differs; use a new label or --refit");else writeJson(plan,expected);
    }
    for(auto& row:rows){auto sample=samplePath(context,label,row);auto collected=sample/"collection";
        if(stage=="collect"||stage=="all"){fs::create_directories(collected);Lock lock(collected/".lock");collectUnit(row,o,context,collected);}
        if(stage=="fit"||stage=="all")for(auto& method:methods)for(int group=first;group<last;++group)fitPart(row,o,context,sample,method,group,groups);
    }
    std::cout<<"Output: "<<context.result/label<<'\n';
}
}
void help(){std::cout<<
    "Compiled ROOT/C++ SPE gain analysis: separate estimators, memory cell 0, bootstrap off\n"
    "  Analyze discover --output TSV [--profile pilot|full] [--beams sps,ps] [--kinds led,electron,muon]\n"
    "  Analyze run --manifest TSV --label NAME [--stage collect|fit|all] [--row UNIT_INDEX]\n"
    "  Analyze regroup --from-label OLD --label NEW [--row UNIT_INDEX] [--dry-run]\n"
    "    --method all|multi_gaussian|envelope_comb|fft|autocorrelation\n"
    "    --channel-groups 4 [--channel-group N] --channels all|ID,... --events-per-file 0\n"
    "    --temperature-bin 0.5 --time-bin 0 --hit-tags 1 --plot-limit 2 --bootstrap 0\n"
    "    --gain-min 6 --gain-max 40 --adc-min 250 --adc-max 1400 --max-span 400\n"
    "    --min-entries 120 --min-cycles 3 --null-replicates 199 --alpha .01\n"
    "    --stability-fraction .05 --seed 20261013 --refit --recollect --dry-run\n"
    "  Analyze compare --label NAME [--manifest TSV] [--allow-partial]\n"
    "    [--stage numeric|prepare|analyze|merge|plots] [--layer 0..29] [--comparison-name NAME]\n"
    "    [--reference ROOT --reference-format led|physical] [--temperature-tolerance .5]\n"
    "    Default compare stage: numeric only. Plots are a separate command.\n"
    "  Analyze submit --manifest TSV --label NAME [--stage all|collect|fit]\n"
    "    [--channel-groups 4 --max-concurrent 16 --memory 32G --walltime 24:00:00]\n"
    "    [--channels all --events-per-file 0 --min-entries 120 --bootstrap 0] [--execute]\n"
    "Default stage: collect only. LED units: beam/date/DAC. --row indexes analysis_units.json.\n"
    "Without --execute, submit only prepares frozen tasks. Positive bootstrap counts are rejected.\n";}
int main(int argc,char** argv)try{
    gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);Context context(argc,argv);auto options=parse(argc,argv);
    if(options.action=="help"||options.action=="--help"||options.action=="-h"||options.has("help")){help();return 0;}
    if(options.action=="discover")discover(options);else if(options.action=="run")gain::run(options,context);else if(options.action=="compare")compare(options,context);else if(options.action=="submit")gain::submit(options,context);else if(options.action=="task")gain::task(options,context);else if(options.action=="regroup")gain::regroup(options,context);
    return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
