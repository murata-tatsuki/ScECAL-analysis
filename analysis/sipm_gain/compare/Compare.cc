#include "CompareIO.hh"
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TChain.h>

namespace gain {
namespace {
using namespace comparison;
std::map<int,Json> reference(const fs::path& path,const std::string& format){
    require(format=="led"||format=="physical","Invalid reference ID format");auto file=openRoot(path);auto* t=tree(*file,"T_Event");
    std::string id=t->GetBranch("cellIDs")?"cellIDs":"CellID";require(t->GetBranch(id.c_str())&&t->GetBranch("gain"),"Unsupported gain reference");std::map<int,Json> result;
    for(Long64_t i=0;i<t->GetEntries();++i){require(t->GetEntry(i)>0,"Reference read failed");int n=leaf(*t,id.c_str());require(n>=0,"Negative reference channel");
        int cell=format=="led"?n/1000000*100000+n%1000000/10000*10000+n%100:n/100000*100000+n%100000/10000*10000+n%100;
        if(cell/100000>=30||(cell%100000/10000==5&&cell%100>=30))continue;
        require(cell%100000/10000<6&&cell%100<36,"Malformed reference channel");require(!result.count(cell),"Duplicate reference channel");
        result[cell]={{"reference_gain",leaf(*t,"gain")},{"reference_error",t->GetLeaf("gainError")?leaf(*t,"gainError"):NaN},{"reference_temperature",t->GetLeaf("temp")?leaf(*t,"temp"):NaN}};
    }return result;
}
std::string groupKey(const Json& r,std::initializer_list<const char*> fields){std::string s;for(auto key:fields){s+=str(r,key);s+='\t';}return s;}
void plot(const Json& rows,const fs::path& out,TFile& root,TsvRows& index){
    root.cd();std::set<int> cells;for(auto& r:rows)cells.insert(int(num(r,"cellid")));const int colors[]={kBlack,kBlue+1,kRed+1,kGreen+2};
    for(int cell:cells){std::vector<const Json*> valid;std::vector<std::pair<std::string,std::string>> spectra;
        for(auto& r:rows)if(int(num(r,"cellid"))==cell&&std::isfinite(num(r,"gain"))){valid.push_back(&r);auto key=std::make_pair(str(r,"sample"),str(r,"sid"));if(std::find(spectra.begin(),spectra.end(),key)==spectra.end())spectra.push_back(key);}
        for(size_t first=0;first<spectra.size();first+=24){size_t count=std::min(size_t(24),spectra.size()-first);std::map<std::pair<std::string,std::string>,int> slots;for(size_t i=0;i<count;++i)slots[spectra[first+i]]=i;
            std::vector<const Json*> points;double high=1;for(auto* r:valid)if(slots.count({str(*r,"sample"),str(*r,"sid")})){points.push_back(r);double error=num(*r,"stat_error");high=std::max(high,num(*r,"gain")+(std::isfinite(error)?error:0.));}
            std::string stem="channel_"+std::to_string(cell)+(spectra.size()>24?"_page"+std::to_string(first/24+1):"");
            TCanvas canvas(stem.c_str(),"Gain estimates",1500,1000);canvas.SetBottomMargin(.36);canvas.SetRightMargin(.25);
            TH1D axis(("axis_"+stem).c_str(),("CellID "+std::to_string(cell)+";;gain [ADC/p.e.]").c_str(),count,-.5,count-.5);axis.SetStats(false);axis.SetMinimum(0);axis.SetMaximum(high*1.1);axis.GetXaxis()->SetLabelSize(.022);
            for(auto& [key,slot]:slots){const Json* r=nullptr;for(auto* p:points)if(str(*p,"sample")==key.first&&str(*p,"sid")==key.second){r=p;break;}require(r,"Missing plot row");
                auto sample=str(*r,"sample");std::smatch match;std::string run=std::regex_search(sample,match,std::regex("Run[0-9]+"))?match[0].str():str(*r,"period");
                std::string tag=str(*r,"kind")=="led"?"LED "+str(*r,"period").substr(0,8)+" DAC"+str(*r,"dac"):str(*r,"kind")+" "+run;
                std::string label=str(*r,"beam")+" "+tag+" s"+str(*r,"sid");axis.GetXaxis()->SetBinLabel(slot+1,label.c_str());
                index.add({{"figure",stem},{"slot",slot},{"cellid",cell},{"sample",str(*r,"sample")},{"sid",str(*r,"sid")},{"temperature_C",num(*r,"temperature_C")}});
            }
            axis.GetXaxis()->LabelsOption("v");axis.Draw();TLegend legend(.76,.54,.995,.92);legend.SetTextSize(.019);std::vector<TGraphErrors> graphs;graphs.reserve(8);
            for(int j=0;j<4;++j)for(auto status:{"candidate","review"}){graphs.emplace_back();auto& g=graphs.back();auto name=Methods[j]+"_"+status+"_"+stem;g.SetName(name.c_str());g.SetMarkerStyle(std::string(status)=="candidate"?20+j:24+j);g.SetMarkerSize(1.2);g.SetMarkerColor(colors[j]);g.SetLineColor(colors[j]);
                for(auto* r:points)if(str(*r,"method")==Methods[j]&&str(*r,"status")==status){int n=g.GetN();g.SetPoint(n,slots.at({str(*r,"sample"),str(*r,"sid")})+(j-1.5)*.15,num(*r,"gain"));double error=num(*r,"stat_error");g.SetPointError(n,0,std::isfinite(error)?error:0);}
                if(g.GetN()){g.Draw("P SAME");legend.AddEntry(&g,(Methods[j]+" "+status).c_str(),"p");g.Write();}
            }legend.Draw();canvas.Write();canvas.SaveAs((out/(stem+".png")).c_str());canvas.Close();
        }
    }
}
Long64_t compareChannel(const Json& rows, double tolerance,
                        const std::function<void(const Json&)>& pair,
                        const std::function<void(const Json&)>& cross) {
    std::map<std::string,std::vector<const Json*>> bySpectrum;
    for(auto& row:rows)bySpectrum[groupKey(row,{"sample","sid"})].push_back(&row);
    for(auto& [key,values]:bySpectrum)for(size_t a=0;a<values.size();++a)for(size_t b=a+1;b<values.size();++b){
        auto& x=*values[a];auto& y=*values[b];pair({{"sample",x["sample"]},{"sid",x["sid"]},{"cellid",x["cellid"]},{"method_a",x["method"]},{"method_b",y["method"]},{"status_a",x["status"]},{"status_b",y["status"]},
            {"gain_a",num(x,"gain")},{"gain_b",num(y,"gain")},{"delta_ADC",num(x,"gain")-num(y,"gain")},{"ratio",num(y,"gain")>0?num(x,"gain")/num(y,"gain"):NaN},
            {"paired_stat_error_ADC",NaN},{"covariance_ADC2",NaN},{"pull",NaN},{"compatibility","not_assessed_no_statistical_uncertainty"}});
    }
    std::map<std::string,std::vector<const Json*>> matched;Long64_t candidates=0;
    for(auto& row:rows)if(str(row,"status")=="candidate"){++candidates;matched[groupKey(row,{"cellid","hit_tag","gain_tag","memory_cell","method"})].push_back(&row);}
    for(auto& [key,values]:matched)for(auto* a:values)for(auto kind:{"led","electron","muon"}){
        const Json* b=nullptr;double best=INFINITY;
        for(auto* test:values)if(str(*test,"kind")==kind&&str(*test,"sample")!=str(*a,"sample")){
            double distance=std::abs(num(*a,"temperature_C")-num(*test,"temperature_C"));if(std::isfinite(distance)&&distance<best){best=distance;b=test;}
        }
        if(!b||groupKey(*a,{"sample","sid"})>groupKey(*b,{"sample","sid"}))continue;
        std::string condition="matched_declared_conditions";
        if(str(*a,"settings")=="unconfirmed"||str(*b,"settings")=="unconfirmed")condition="settings_unconfirmed";
        else if(str(*a,"settings")!=str(*b,"settings"))condition="settings_differ";
        if(best>tolerance)condition+=";temperature_mismatch";
        cross({{"cellid",(*a)["cellid"]},{"method",(*a)["method"]},{"sample_a",(*a)["sample"]},{"sid_a",(*a)["sid"]},{"sample_b",(*b)["sample"]},{"sid_b",(*b)["sid"]},
            {"delta_temperature_C",num(*a,"temperature_C")-num(*b,"temperature_C")},{"delta_gain_ADC",num(*a,"gain")-num(*b,"gain")},{"ratio",num(*b,"gain")>0?num(*a,"gain")/num(*b,"gain"):NaN},
            {"condition_status",condition},{"compatibility","not_assessed_no_statistical_uncertainty"}});
    }

    return candidates;
}

struct FitInputs {
    Json plan, source, samples;
    int groups;
};
FitInputs fitInputs(const Options& o, const Context& c) {
    auto inputs = c.result / token(o.required("label")) / "inputs";
    require(fs::exists(inputs / "analysis_plan.json"), "No split-analysis plan; complete the fit stage first");
    FitInputs f;
    f.plan = readJson(inputs / "analysis_plan.json"); f.groups = f.plan.at("channel_groups");
    require(f.groups > 0 && f.groups <= 6300, "Invalid fit channel groups");
    f.source = readJson(inputs / "source_manifest.json"); f.samples = readTable(inputs / "manifest.tsv");
    if (o.has("manifest")) require(readTable(o.required("manifest")) == f.samples, "Comparison manifest differs from analysis inputs");
    if (fs::exists(inputs / "analysis_units.json")) {
        auto units = readJson(inputs / "analysis_units.json");
        require(units == analysisUnits(f.samples), "Analysis unit manifest changed"); f.samples = units;
    }
    return f;
}
// Small metadata only. Input contents are checksummed by verifyPart during prepare.
Json inputSnapshot(const Context& c, const std::string& label, const FitInputs& f) {
    auto base = c.result / label; Json snapshot = {{"inputs", Json::object()}, {"parts", Json::array()}, {"collections", Json::object()}};
    for (auto name : {"manifest.tsv", "analysis_plan.json", "source_manifest.json", "analysis_units.json"})
        if (fs::exists(base / "inputs" / name)) snapshot["inputs"][name] = sha256(base / "inputs" / name);
    for (auto& sample : f.samples) {
        auto sampleDir = samplePath(c, label, sample); auto collect = sampleDir / "collection/COLLECT_COMPLETE";
        snapshot["collections"][str(sample, "sample")] = fs::exists(collect) ? sha256(collect) : "missing";
        for (auto& method : Methods) for (int group = 0; group < f.groups; ++group) {
            auto path = sampleDir / method / partName(group, f.groups);
            Json entry = {{"path", fs::relative(path, base).string()}};
            if (fs::exists(path / "FIT_COMPLETE")) {
                entry["marker_sha256"] = sha256(path / "FIT_COMPLETE");
                for (auto name : {"gain.root", "gain_estimates.tsv", "fits.json", "diagnostics.root"}) entry["files"][name] = fingerprint(path / name);
            } else entry["missing"] = true;
            snapshot["parts"].push_back(std::move(entry));
        }
    }
    return snapshot;
}
std::string configurationHash(const fs::path& out) { return sha256(out / "inputs/configuration.json"); }
void copySources(const Context& c, const fs::path& out, const Json& sources) {
    const std::string prefix = "analysis/sipm_gain/";
    for (auto it = sources.begin(); it != sources.end(); ++it) {
        fs::path from = it.key().rfind(prefix, 0) == 0 ? c.here / it.key().substr(prefix.size()) : c.project / it.key();
        auto dest = out / "inputs/source" / it.key(); fs::create_directories(dest.parent_path());
        fs::copy_file(from, dest, fs::copy_options::overwrite_existing);
        require(sha256(dest) == it.value(), "Comparison source changed while snapshotting");
    }
}
Json prepare(const Options& o, const Context& c, const fs::path& out) {
    auto label = token(o.required("label")); auto f = fitInputs(o, c);
    double tolerance = o.real("temperature-tolerance", .5); require(tolerance >= 0, "Invalid temperature tolerance");
    require(!o.has("reference-temperature"), "Temperature correction is not implemented; omit --reference-temperature");
    require(!o.has("reference-format") || o.has("reference"), "--reference-format requires --reference");
    Json config = {{"format_version", 1}, {"label", label}, {"temperature_tolerance", tolerance},
        {"allow_partial", o.has("allow-partial")}, {"comparison_sources", c.sources()}, {"input_snapshot", inputSnapshot(c, label, f)}};
    if (o.has("reference")) {
        config["reference"] = fingerprint(o.required("reference")); config["reference_sha256"] = sha256(o.required("reference"));
        config["reference_format"] = o.get("reference-format", "led");
    }
    auto configPath = out / "inputs/configuration.json";
    if (fs::exists(configPath)) require(readJson(configPath) == config, "Inputs, source or options changed; use a new --comparison-name");
    else {
        require(!fs::exists(out / "all_estimates.root") && !fs::exists(out / "summary.json") && !fs::exists(out / "comparison.root"), "Legacy comparison exists; preserve it with a new --comparison-name");
        fs::create_directories(out / "inputs"); copySources(c, out, config.at("comparison_sources")); atomicJson(configPath, config);
    }
    auto signature = configurationHash(out);
    if (fs::exists(out / "PREPARE_COMPLETE")) {
        // Recheck original fit checksums as well as the frozen prepared data on explicit prepare reruns.
        for (auto& sample : f.samples) for (auto& method : Methods) for (int group = 0; group < f.groups; ++group) {
            auto path = samplePath(c, label, sample) / method / partName(group, f.groups);
            if (fs::exists(path / "FIT_COMPLETE")) { Lock lock(path / ".lock", false, true); verifyPart(path); }
        }
        verifyComplete(out, "PREPARE_COMPLETE", signature); std::cout << "Reuse prepared comparison: " << out << std::endl; return config;
    }
    std::map<int, Json> ref;
    if (o.has("reference")) ref = reference(o.required("reference"), o.get("reference-format", "led"));
    fs::create_directories(out / "prepared");
    std::vector<std::unique_ptr<TFile>> files; std::vector<std::unique_ptr<RootRows>> writers;
    for (int layer = 0; layer < 30; ++layer) {
        files.push_back(openRoot(out / "prepared" / (layerName(layer) + ".root"), "RECREATE"));
        writers.push_back(std::make_unique<RootRows>(*files.back(), "Gain", gainSchema(true)));
    }
    Json missing = Json::array(); Long64_t total = 0; size_t unit = 0;
    for (auto& sample : f.samples) {
        // Bounded to one input unit; sids remain unique across all methods and channel parts.
        std::set<std::string> seen; Json common; auto name = str(sample, "sample");
        for (auto& method : Methods) for (int group = 0; group < f.groups; ++group) {
            auto path = samplePath(c, label, sample) / method / partName(group, f.groups);
            if (!fs::exists(path / "FIT_COMPLETE")) { missing.push_back(fs::relative(path, c.result / label).string()); continue; }
            Lock partLock(path / ".lock", false, true); auto settings = verifyPart(path).at("config");
            require(settings.at("source") == f.source, "Part uses an older source snapshot");
            require(str(settings, "method") == method && num(settings, "group") == group && num(settings, "groups") == f.groups, "Part identity mismatch");
            require(settings.at("config") == f.plan.at("fit_config") && settings.at("hit_tags") == f.plan.at("hit_tags"), "Mixed fit settings in comparison");
            require(settings.at("spectra_sha256") == trim(readText(samplePath(c, label, sample) / "collection/COLLECT_COMPLETE")), "Collection changed after fitting; refit affected parts");
            Json identity = {{"source", settings.at("source")}, {"spectra_sha256", settings.at("spectra_sha256")}};
            if (common.is_null()) common = identity; else require(common == identity, "Methods used different source or collected data");
            visitRows(path / "gain.root", [&](Json& row) {
                require(str(row, "sample") == name && str(row, "method") == method && str(row, "kind") == str(sample, "kind"), "ROOT record identity differs");
                int cell = physicalCell(row);
                require(num(row, "memory_cell") == 0 && channelGroup(cell, f.groups) == group, "ROOT record in wrong channel group");
                require(seen.insert(groupKey(row, {"sample", "sid", "method"})).second, "Duplicate gain record");
                if (ref.count(cell)) row.update(ref.at(cell));
                double g = num(row, "reference_gain"); row["gain_over_reference"] = g > 0 ? num(row, "gain") / g : NaN;
                row["reference_condition"] = "settings_unconfirmed";
                // Reproduce original method ordering and nearest-temperature tie breaking.
                row["_ordinal"] = total++; writers.at(cell / 100000)->add(row);
            });
        }
        std::cout << "prepare unit " << ++unit << '/' << f.samples.size() << ": " << name << ", rows=" << total << std::endl;
    }
    require(missing.empty() || o.has("allow-partial"), "Incomplete method/channel parts; use --allow-partial for a marked partial comparison");
    require(total > 0, "No completed gain records");
    require(inputSnapshot(c, label, f) == config.at("input_snapshot"), "Inputs changed while preparing comparison; use a new --comparison-name");
    Json counts = Json::object(); std::vector<std::string> outputs;
    for (int layer = 0; layer < 30; ++layer) {
        auto name = layerName(layer); counts[name] = writers[layer]->size();
        writers[layer]->finish(true); writers[layer].reset();
        metadata(*files[layer], {{"configuration_sha256", signature}, {"layer", layer}, {"rows", counts[name]}});
        closeRoot(*files[layer]); files[layer].reset(); outputs.push_back("prepared/" + name + ".root");
    }
    seal(out, "PREPARE_COMPLETE", signature, outputs, {{"rows", total}, {"layers", counts}, {"missing", missing}, {"partial", !missing.empty()}, {"command", commandLine(c.command)}});
    return config;
}
Json loadConfiguration(const Options& o, const Context& c, const fs::path& out) {
    require(fs::exists(out / "PREPARE_COMPLETE"), "Run compare --stage prepare first");
    auto config = readJson(out / "inputs/configuration.json");
    require(config.at("format_version") == 1 && config.at("label") == o.required("label"), "Incompatible comparison configuration");
    require(config.at("comparison_sources") == c.sources(), "Comparison source changed; use a new --comparison-name and prepare");
    if (o.has("temperature-tolerance")) require(o.real("temperature-tolerance", .5) == num(config, "temperature_tolerance"), "Temperature tolerance differs from prepare");
    if (o.has("reference")) require(config.contains("reference_sha256") && sha256(o.required("reference")) == config.at("reference_sha256"), "Reference differs from prepare");
    if (o.has("reference-format")) require(o.get("reference-format") == str(config, "reference_format"), "Reference format differs from prepare");
    require(!o.has("reference-temperature"), "Temperature correction is not implemented");
    require(!o.has("manifest"), "--manifest is only used with prepare/numeric");
    if (o.has("allow-partial")) require(config.at("allow_partial") == true, "--allow-partial differs from prepare");
    return config;
}
Json preparedInfo(const fs::path& out, const std::string& signature, int layer) {
    auto info = readJson(out / "PREPARE_COMPLETE");
    require(str(info, "configuration_sha256") == signature, "Prepared configuration changed");
    auto path = "prepared/" + layerName(layer) + ".root";
    require(sha256(out / path) == info.at("files").at(path), "Prepared ROOT changed: " + path);
    return info;
}
const std::vector<std::string> NumericFiles = {"numeric.root", "all_estimates.tsv", "method_comparison.tsv", "cross_sample.tsv"};
void analyzeLayer(const Options& o, const Context& c, const fs::path& out, const Json& config, int layer) {
    auto signature = configurationHash(out); auto prepared = preparedInfo(out, signature, layer);
    auto name = layerName(layer); auto dest = out / "layers" / name; fs::create_directories(dest); Lock lock(dest / ".lock");
    auto preparedSHA = prepared.at("files").at("prepared/" + name + ".root");
    if (fs::exists(dest / "NUMERIC_COMPLETE")) {
        auto info = verifyComplete(dest, "NUMERIC_COMPLETE", signature);
        require(info.at("prepared_sha256") == preparedSHA, "Layer uses another prepared ROOT");
        std::cout << "Reuse numerical comparison: " << name << std::endl; return;
    }
    auto file = openRoot(dest / "numeric.root", "RECREATE");
    RootRows gain(*file, "Gain", gainSchema()), pairs(*file, "MethodComparison", PairSchema), cross(*file, "CrossSample", CrossSchema);
    auto columns = gainSchema().columns(); std::sort(columns.begin(), columns.end());
    TsvRows gt(dest / "all_estimates.tsv", columns), pt(dest / "method_comparison.tsv", PairColumns), ct(dest / "cross_sample.tsv", CrossColumns);
    Long64_t candidates = 0, channels = 0;
    visitChannels(out / "prepared" / (name + ".root"), [&](const Json& rows) {
        require(physicalCell(rows.at(0)) / 100000 == layer, "Wrong layer in prepared ROOT");
        for (auto& row : rows) { gain.add(row); gt.add(row); }
        candidates += compareChannel(rows, num(config, "temperature_tolerance"),
            [&](const Json& row) { pairs.add(row); pt.add(row); }, [&](const Json& row) { cross.add(row); ct.add(row); });
        if (++channels % 30 == 0) std::cout << name << ": channels=" << channels << ", rows=" << gain.size() << std::endl;
    });
    Json info = {{"layer", layer}, {"rows", gain.size()}, {"candidates", candidates}, {"method_pairs", pairs.size()}, {"cross_sample_pairs", cross.size()},
        {"prepared_sha256", preparedSHA}, {"command", commandLine(c.command)}};
    require(gain.size() == prepared.at("layers").at(name).get<Long64_t>(), "Layer row count differs");
    gain.finish(); pairs.finish(); cross.finish(); metadata(*file, info); closeRoot(*file); gt.finish(); pt.finish(); ct.finish();
    seal(dest, "NUMERIC_COMPLETE", signature, NumericFiles, info);
    std::cout << "Completed numerical comparison: " << name << ", rows=" << info["rows"] << std::endl;
    (void)o;
}
void mergeTsv(const fs::path& out, const std::string& name) {
    std::ofstream target(out / name); require(bool(target), "Cannot create merged TSV"); std::string expected;
    for (int layer = 0; layer < 30; ++layer) {
        std::ifstream input(out / "layers" / layerName(layer) / name); std::string line;
        require(bool(std::getline(input, line)), "Missing TSV header");
        if (layer == 0) { expected = line; target << line << '\n'; } else require(line == expected, "Different layer TSV schemas");
        while (std::getline(input, line)) target << line << '\n';
        require(input.eof() && bool(target), "Merged TSV IO error");
    }
    target.close(); require(bool(target), "Merged TSV close failed");
}
void mergeLayers(const Context& c, const fs::path& out, const Json& config) {
    Lock mergeLock(out / ".merge.lock"); auto signature = configurationHash(out); auto prepared = readJson(out / "PREPARE_COMPLETE");
    require(str(prepared, "configuration_sha256") == signature, "Prepared configuration changed");
    std::vector<std::unique_ptr<Lock>> locks; Json layerChecks = Json::object(), layerCounts = Json::object(); Long64_t rows = 0, candidates = 0, pairs = 0, cross = 0;
    for (int layer = 0; layer < 30; ++layer) {
        auto name = layerName(layer); auto dest = out / "layers" / name;
        require(fs::exists(dest / "NUMERIC_COMPLETE"), "Numerical layer unfinished: " + name);
        locks.push_back(std::make_unique<Lock>(dest / ".lock", false, true));
        auto info = verifyComplete(dest, "NUMERIC_COMPLETE", signature);
        require(info.at("prepared_sha256") == prepared.at("files").at("prepared/" + name + ".root"), "Prepared layer identity differs");
        layerChecks[name] = sha256(dest / "NUMERIC_COMPLETE");
        layerCounts[name] = {{"Gain", info.at("rows")}, {"MethodComparison", info.at("method_pairs")}, {"CrossSample", info.at("cross_sample_pairs")}};
        rows += info.at("rows").get<Long64_t>(); candidates += info.at("candidates").get<Long64_t>();
        pairs += info.at("method_pairs").get<Long64_t>(); cross += info.at("cross_sample_pairs").get<Long64_t>();
    }
    require(rows == prepared.at("rows").get<Long64_t>(), "Merged row count differs");
    if (fs::exists(out / "NUMERIC_COMPLETE")) {
        auto existing = verifyComplete(out, "NUMERIC_COMPLETE", signature);
        require(existing.at("layers") == layerChecks, "Completed layer comparison changed");
        std::cout << "Reuse merged comparison: " << out << std::endl; return;
    }
    Json summary = {{"rows", rows}, {"candidates", candidates}, {"method_pairs", pairs}, {"cross_sample_pairs", cross},
        {"missing", prepared.at("missing")}, {"partial", prepared.at("partial")}, {"bootstrap", 0}, {"statistical_compatibility", "not_evaluated"},
        {"no_automatic_calibration_export", true}, {"command", commandLine(c.command)}, {"configuration_sha256", signature}, {"plots", "separate_stage"},
        {"comparison_source_sha256", config.at("comparison_sources").at("analysis/sipm_gain/compare/Compare.cc")}};
    for (auto key : {"reference", "reference_sha256"}) if (config.contains(key)) summary[key] = config.at(key);
    auto output = openRoot(out / "all_estimates.root", "RECREATE");
    for (auto name : {"Gain", "MethodComparison", "CrossSample"}) {
        TChain chain(name);
        for (int layer = 0; layer < 30; ++layer) {
            auto count = layerCounts.at(layerName(layer)).at(name).get<Long64_t>();
            if (count) require(chain.Add((out / "layers" / layerName(layer) / "numeric.root").c_str(), count) == 1, "Cannot add layer TTree");
        }
        output->cd();
        // CloneTree on a chain containing only empty trees can crash in ROOT.
        if (!chain.GetEntries()) {
            Schema schema = std::string(name) == "Gain" ? gainSchema() : std::string(name) == "MethodComparison" ? PairSchema : CrossSchema;
            RootRows empty(*output, name, schema); empty.finish(); continue;
        }
        auto* merged = chain.CloneTree(-1, "fast");
        require(merged && merged->GetEntries() == chain.GetEntries(), "Merged TTree entry count differs");
        require(merged->Write() > 0, "Cannot save merged TTree"); merged->SetDirectory(nullptr); delete merged;
    }
    metadata(*output, summary); closeRoot(*output);
    for (auto name : {"all_estimates.tsv", "method_comparison.tsv", "cross_sample.tsv"}) mergeTsv(out, name);
    writeJson(out / "summary.json", summary);
    seal(out, "NUMERIC_COMPLETE", signature, {"all_estimates.root", "all_estimates.tsv", "method_comparison.tsv", "cross_sample.tsv", "summary.json"}, {{"layers", layerChecks}, {"rows", rows}});
    std::cout << "Completed merged comparison: " << out << std::endl;
}
void plotLayer(const Context& c, const fs::path& out, int layer) {
    auto signature = configurationHash(out); auto prepared = preparedInfo(out, signature, layer); auto name = layerName(layer);
    auto numbers = out / "layers" / name; require(fs::exists(numbers / "NUMERIC_COMPLETE"), "Run numerical comparison for " + name + " first");
    Lock numericLock(numbers / ".lock", false, true); auto numeric = verifyComplete(numbers, "NUMERIC_COMPLETE", signature);
    require(numeric.at("prepared_sha256") == prepared.at("files").at("prepared/" + name + ".root"), "Numerical and prepared layers differ");
    auto dest = out / "plots" / name; fs::create_directories(dest); Lock lock(dest / ".lock");
    if (fs::exists(dest / "PLOTS_COMPLETE")) { verifyComplete(dest, "PLOTS_COMPLETE", signature); std::cout << "Reuse plots: " << name << std::endl; return; }
    auto root = openRoot(dest / "comparison.root", "RECREATE"); TsvRows index(dest / "plot_spectra.tsv", PlotColumns);
    visitChannels(out / "prepared" / (name + ".root"), [&](const Json& rows) { plot(rows, dest, *root, index); });
    metadata(*root, {{"configuration_sha256", signature}, {"layer", layer}, {"command", commandLine(c.command)}}); closeRoot(*root); index.finish();
    std::vector<std::string> files = {"comparison.root", "plot_spectra.tsv"};
    for (auto& entry : fs::directory_iterator(dest)) if (entry.is_regular_file() && entry.path().extension() == ".png") files.push_back(entry.path().filename().string());
    seal(dest, "PLOTS_COMPLETE", signature, files, {{"layer", layer}, {"command", commandLine(c.command)}});
    std::cout << "Completed plots: " << name << std::endl;
}
}
void compare(const Options& o, const Context& c) {
    auto stage = o.get("stage", "numeric");
    require(stage == "numeric" || stage == "prepare" || stage == "analyze" || stage == "merge" || stage == "plots", "Invalid comparison stage");
    int first = 0, last = 30;
    if (o.has("layer")) {
        require(stage == "analyze" || stage == "plots", "--layer is only valid for analyze or plots");
        auto layer = o.integer("layer", -1); require(layer >= 0 && layer < 30, "--layer must be 0..29"); first = layer; last = first + 1;
    }
    auto out = c.result / token(o.required("label")) / token(o.get("comparison-name", "comparison"));
    // Keep result directories disjoint from fit inputs and per-method output.
    require(out.filename() != "inputs" && out.filename() != "data", "Reserved --comparison-name");
    fs::create_directories(out); Lock lock(out / ".lock", false, stage != "numeric" && stage != "prepare");
    Json config = (stage == "numeric" || stage == "prepare") ? prepare(o, c, out) : loadConfiguration(o, c, out);
    if (stage == "numeric" || stage == "analyze") for (int layer = first; layer < last; ++layer) analyzeLayer(o, c, out, config, layer);
    if (stage == "numeric" || stage == "merge") mergeLayers(c, out, config);
    if (stage == "plots") for (int layer = first; layer < last; ++layer) plotLayer(c, out, layer);
}
}
