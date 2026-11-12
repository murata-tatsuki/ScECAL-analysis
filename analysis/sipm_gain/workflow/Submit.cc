#include "Pipeline.hh"
#include <sys/stat.h>
gain::Config configuration(const gain::Options&);
void channelsValid(const std::string&);
namespace gain {
void task(const Options& o,const Context& c){
    auto plan=readJson(o.required("tasks"));auto index=o.integer("index",-1);
    require(index>=0&&index<static_cast<long long>(plan.at("tasks").size()),"Task index out of range");
    require(sha256(c.here/"Analyze")==plan.at("analyze_sha256")&&sha256(c.here/"Collect")==plan.at("collect_sha256"),"Frozen executable changed");
    require(sha256(plan.at("manifest").get<std::string>())==plan.at("manifest_sha256"),"Frozen manifest changed");
    Options next;next.action="run";next.values=plan.at("options").get<std::map<std::string,std::string>>();
    const auto& item=plan.at("tasks").at(index);next.values["stage"]=str(item,"stage");next.values["row"]=str(item,"row");
    if(str(item,"stage")=="fit"){next.values["method"]=str(item,"method");next.values["channel-group"]=str(item,"group");}
    run(next,c);
}
void submit(const Options& o,const Context& c){
    auto rawRows=manifest(o.required("manifest"));auto rows=analysisUnits(rawRows);auto label=token(o.required("label"));auto stage=o.get("stage","collect");
    require(stage=="all"||stage=="collect"||stage=="fit","Invalid submission stage");
    int groups=o.integer("channel-groups",4),concurrent=o.integer("max-concurrent",16);
    require(groups>0&&groups<=6300&&concurrent>0,"Invalid parallel settings");
    size_t totalRequested=rows.size()*((stage!="fit"?1:0)+(stage!="collect"?4*groups:0));
    require(totalRequested<=35000,"Too many pending tasks for one submission; reduce channel-groups or split the manifest (local MaxJobCount is 40000)");
    auto dest=c.log/label;fs::create_directories(dest);Lock lock(dest/".submit.lock");
    Json previousJobs=fs::exists(dest/"jobs.tsv")?readTable(dest/"jobs.tsv"):Json::array();
    if(!previousJobs.empty()){require(stage=="fit","This label already has a submission record");for(auto& job:previousJobs)require(str(job,"stage")=="collect","Fit jobs are already recorded for this label");}
    Options settings;settings.action="run";
    for(auto& [key,value]:o.values)if(key!="partition"&&key!="max-concurrent"&&key!="memory"&&key!="walltime"&&key!="stage")settings.values[key]=value;
    settings.values["channels"]=o.get("channels","all");settings.values["events-per-file"]=o.get("events-per-file","0");
    settings.values["bootstrap"]=o.get("bootstrap","0");settings.values["min-entries"]=o.get("min-entries","120");
    settings.values["channel-groups"]=std::to_string(groups);configuration(settings);channelsValid(settings.values["channels"]);
    if(stage=="fit")for(auto& row:rows)require(fs::exists(samplePath(c,label,row)/"collection/COLLECT_COMPLETE"),"Collect all inputs before submit --stage fit");
    freezeInputs(c,label,rawRows);
    auto frozen=dest/"manifest.tsv";writeTable(frozen,rawRows,ManifestFields);settings.values["manifest"]=frozen.string();
    auto runtime=dest/"runtime";fs::create_directories(runtime);
    for(auto& entry:fs::recursive_directory_iterator(c.here)){
        if(!entry.is_regular_file())continue;
        auto path=entry.path();auto ext=path.extension().string();
        if(!(ext==".cc"||ext==".hh"||ext==".sh"||path.filename()=="Makefile"||path.filename()=="Analyze"||path.filename()=="Collect"))continue;
        auto target=runtime/fs::relative(path,c.here);fs::create_directories(target.parent_path());if(fs::exists(target))require(sha256(path)==sha256(target),"Prepared runtime differs; use a new label");else fs::copy_file(path,target);
    }
    require(fs::is_regular_file(runtime/"Analyze")&&fs::is_regular_file(runtime/"Collect"),"Build executables before submission");
    Json base={{"manifest",frozen.string()},{"manifest_sha256",sha256(frozen)},{"options",settings.values},{"analyze_sha256",sha256(runtime/"Analyze")},{"collect_sha256",sha256(runtime/"Collect")}};
    std::vector<std::pair<std::string,Json>> stages;
    if(stage!="fit"){Json tasks=Json::array();for(size_t row=0;row<rows.size();++row)tasks.push_back({{"stage","collect"},{"row",row}});stages.emplace_back("collect",tasks);}
    if(stage!="collect"){Json tasks=Json::array();for(size_t row=0;row<rows.size();++row)for(int group=0;group<groups;++group)for(auto& method:Methods)tasks.push_back({{"stage","fit"},{"row",row},{"group",group},{"method",method}});stages.emplace_back("fit",tasks);}
    std::string previousJob;Json jobs=previousJobs;std::ostringstream script;script<<"#!/bin/bash\nset -euo pipefail\n";
    size_t taskCount=0;int wave=previousJobs.size(),firstWave=wave;
    for(auto& [name,tasks]:stages)for(size_t first=0;first<tasks.size();first+=1000){
        size_t n=std::min(size_t(1000),tasks.size()-first);Json wavePlan=base;wavePlan["tasks"]=Json::array();for(size_t i=0;i<n;++i)wavePlan["tasks"].push_back(tasks[first+i]);
        auto waveDir=dest/name/("wave_"+std::to_string(wave));fs::create_directories(waveDir);auto taskFile=waveDir/"tasks.json";writeJson(taskFile,wavePlan);
        std::vector<std::string> prefix={"sbatch","--parsable","--job-name","sipm_gain_"+name,"--partition",o.get("partition","all"),"--array","0-"+std::to_string(n-1)+"%"+std::to_string(concurrent),"--mem",o.get("memory","32G"),"--time",o.get("walltime","24:00:00"),"--cpus-per-task=1","--output",(waveDir/"slurm-%A_%a.out").string(),"--error",(waveDir/"slurm-%A_%a.err").string()};
        std::vector<std::string> suffix={(runtime/"slurm_array.sh").string(),runtime.string(),c.project.string(),taskFile.string(),sha256(taskFile)};
        script<<"job_"<<wave<<"=$("<<commandLine(prefix);if(wave>firstWave)script<<" --dependency=afterok:${job_"<<wave-1<<"}";script<<' '<<commandLine(suffix)<<")\nprintf '%s\\n' \"$job_"<<wave<<"\"\n";
        if(o.has("execute")){
            // Record every accepted wave immediately; a later failure cannot hide live jobs.
            auto cmd=prefix;if(!previousJob.empty()){cmd.push_back("--dependency");cmd.push_back("afterok:"+previousJob);}cmd.insert(cmd.end(),suffix.begin(),suffix.end());
            auto job=trim(execute(cmd));require(std::regex_match(job,std::regex("[0-9]+(;[A-Za-z0-9_.-]+)?")),"Unexpected sbatch response: "+job);previousJob=job.substr(0,job.find(';'));
            jobs.push_back({{"job",previousJob},{"stage",name},{"wave",wave},{"tasks",n},{"task_manifest",taskFile.string()},{"sha256",sha256(taskFile)}});writeTable(dest/"jobs.tsv",jobs);std::cout<<"Submitted "<<previousJob<<" "<<name<<" tasks="<<n<<std::endl;
        }
        taskCount+=n;++wave;
    }
    writeText(dest/"submit_command.sh",script.str());writeJson(dest/"submission.json",{{"stage",stage},{"inputs",rawRows.size()},{"analysis_units",rows.size()},{"channel_groups",groups},{"task_count",taskCount},{"waves",wave-firstWave},{"max_concurrent",concurrent},{"bootstrap",0},{"frozen_runtime",runtime.string()},{"command",commandLine(c.command)}});
    std::cout<<"Prepared "<<taskCount<<" tasks in "<<wave-firstWave<<" dependent waves; at most "<<concurrent<<" active tasks. "<<(o.has("execute")?"Submitted.":"No jobs submitted. Add --execute to submit.")<<'\n';
}
}
