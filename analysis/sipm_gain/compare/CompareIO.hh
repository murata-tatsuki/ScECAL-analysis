#pragma once
#include "Pipeline.hh"
#include <TTreeIndex.h>
#include <functional>
#include <limits>

namespace gain::comparison {
struct Schema {
    std::vector<std::string> strings, numbers;
    std::vector<std::string> columns() const {
        auto result = strings;
        result.insert(result.end(), numbers.begin(), numbers.end());
        return result;
    }
};
inline Schema gainSchema(bool prepared = false) {
    Schema s{StringFields, NumericFields};
    s.strings.push_back("reference_condition");
    for (auto key : {"reference_gain", "reference_error", "reference_temperature", "gain_over_reference"}) s.numbers.push_back(key);
    if (prepared) s.numbers.push_back("_ordinal");
    return s;
}
inline const Schema PairSchema{
    {"sample", "method_a", "method_b", "status_a", "status_b", "compatibility"},
    {"sid", "cellid", "gain_a", "gain_b", "delta_ADC", "ratio", "paired_stat_error_ADC", "covariance_ADC2", "pull"}};
inline const Schema CrossSchema{
    {"method", "sample_a", "sample_b", "condition_status", "compatibility"},
    {"cellid", "sid_a", "sid_b", "delta_temperature_C", "delta_gain_ADC", "ratio"}};
inline const std::vector<std::string> PairColumns = {"sample", "sid", "cellid", "method_a", "method_b", "status_a", "status_b", "gain_a", "gain_b", "delta_ADC", "ratio", "paired_stat_error_ADC", "covariance_ADC2", "pull", "compatibility"};
inline const std::vector<std::string> CrossColumns = {"cellid", "method", "sample_a", "sid_a", "sample_b", "sid_b", "delta_temperature_C", "delta_gain_ADC", "ratio", "condition_status", "compatibility"};
inline const std::vector<std::string> PlotColumns = {"figure", "slot", "cellid", "sample", "sid", "temperature_C"};
inline std::string layerName(int layer) {
    std::ostringstream s; s << "layer_" << std::setfill('0') << std::setw(2) << layer; return s.str();
}
inline int physicalCell(const Json& row) {
    double value = num(row, "cellid");
    require(std::isfinite(value) && value >= 0 && value < 3000000 && std::floor(value) == value, "Invalid comparison CellID");
    int cell = value, chip = cell % 100000 / 10000, ch = cell % 100;
    require(chip < 6 && ch < 36 && cell % 10000 / 100 == 0 && !(chip == 5 && ch >= 30), "Disconnected or malformed comparison CellID");
    return cell;
}
// Fixed branch buffers; never hold the dataset or serialized tables in memory.
class RootRows {
    TFile& file_;
    TTree tree_;
    Schema schema_;
    std::map<std::string, std::string> strings_;
    std::map<std::string, double> numbers_;
public:
    RootRows(TFile& file, const char* name, Schema schema): file_(file), tree_(name, name), schema_(std::move(schema)) {
        tree_.SetDirectory(&file_); tree_.SetAutoFlush(4096); tree_.SetAutoSave(0);
        TTree::SetMaxTreeSize(std::numeric_limits<Long64_t>::max());
        for (auto& key : schema_.strings) tree_.Branch(key.c_str(), &strings_[key]);
        for (auto& key : schema_.numbers) tree_.Branch(key.c_str(), &numbers_[key]);
    }
    void add(const Json& row) {
        for (auto& key : schema_.strings) strings_[key] = str(row, key);
        for (auto& key : schema_.numbers) numbers_[key] = num(row, key);
        require(tree_.Fill() >= 0 && !file_.TestBit(TFile::kWriteError), "ROOT comparison write failed");
    }
    Long64_t size() const { return tree_.GetEntries(); }
    void finish(bool index = false) {
        file_.cd();
        if (index && size()) {
            require(tree_.BuildIndex("cellid", "_ordinal") >= 0, "Cannot index comparison channels");
            auto* idx = dynamic_cast<TTreeIndex*>(tree_.GetTreeIndex());
            require(idx && idx->GetN() == size(), "Incomplete comparison channel index");
        }
        require(tree_.Write() > 0 && !file_.TestBit(TFile::kWriteError), "Cannot save comparison TTree");
        tree_.SetDirectory(nullptr);
    }
};
class TsvRows {
    std::ofstream file_;
    std::vector<std::string> fields_;
public:
    TsvRows(const fs::path& path, std::vector<std::string> fields): file_(path), fields_(std::move(fields)) {
        require(bool(file_), "Cannot write " + path.string());
        for (size_t i = 0; i < fields_.size(); ++i) file_ << (i ? "\t" : "") << fields_[i];
        file_ << '\n';
    }
    void add(const Json& row) {
        for (size_t i = 0; i < fields_.size(); ++i) file_ << (i ? "\t" : "") << (row.contains(fields_[i]) ? tsvValue(row[fields_[i]]) : "");
        file_ << '\n'; require(bool(file_), "Comparison TSV write failed");
    }
    void finish() { file_.close(); require(bool(file_), "Comparison TSV close failed"); }
};
class Reader {
    TTree& tree_;
    Schema schema_;
    std::map<std::string, std::string*> strings_;
public:
    Reader(TTree& t, Schema schema): tree_(t), schema_(std::move(schema)) {
        for (auto& key : schema_.strings) {
            strings_[key] = nullptr;
            if (tree_.GetBranch(key.c_str())) require(tree_.SetBranchAddress(key.c_str(), &strings_[key]) >= 0, "Invalid string branch " + key);
            else require(key == "source_sample" || key == "source_raw" || key == "reference_condition", "Missing string branch " + key);
        }
    }
    ~Reader() { tree_.ResetBranchAddresses(); for (auto& item : strings_) delete item.second; }
    Json at(Long64_t entry) {
        require(tree_.GetEntry(entry) > 0, "Cannot read comparison record"); Json row;
        for (auto& key : schema_.strings) row[key] = strings_[key] ? *strings_[key] : "";
        for (auto& key : schema_.numbers) row[key] = tree_.GetLeaf(key.c_str()) ? leaf(tree_, key.c_str()) : NaN;
        return row;
    }
};
inline void visitRows(const fs::path& path, const std::function<void(Json&)>& consume) {
    auto file = openRoot(path); auto* t = tree(*file, "Gain"); Reader reader(*t, {StringFields, NumericFields});
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { auto row = reader.at(i); consume(row); }
}
inline void visitChannels(const fs::path& path, const std::function<void(const Json&)>& consume) {
    auto file = openRoot(path); auto* t = tree(*file, "Gain");
    if (!t->GetEntries()) return;
    auto* index = dynamic_cast<TTreeIndex*>(t->GetTreeIndex());
    require(index && index->GetN() == t->GetEntries(), "Missing prepared channel index");
    Reader reader(*t, gainSchema(true)); Json rows = Json::array(); int previous = -1; double ordinal = -1;
    for (Long64_t i = 0; i < index->GetN(); ++i) {
        auto row = reader.at(index->GetIndex()[i]); int cell = physicalCell(row); double current = num(row, "_ordinal");
        require(cell >= previous && std::isfinite(current) && (cell != previous || current > ordinal), "Invalid prepared record order");
        if (cell != previous && !rows.empty()) { consume(rows); rows = Json::array(); }
        previous = cell; ordinal = current; rows.push_back(std::move(row));
    }
    if (!rows.empty()) consume(rows);
}
inline void atomicJson(const fs::path& path, const Json& value) {
    auto temp = path; temp += ".tmp"; writeJson(temp, value); fs::rename(temp, path);
}
inline Json fileChecks(const fs::path& out, const std::vector<std::string>& names) {
    Json checks = Json::object(); for (auto& name : names) checks[name] = sha256(out / name); return checks;
}
inline Json verifyComplete(const fs::path& out, const char* marker, const std::string& configuration) {
    auto info = readJson(out / marker); require(str(info, "configuration_sha256") == configuration, "Comparison configuration changed: " + out.string());
    for (auto it = info.at("files").begin(); it != info.at("files").end(); ++it)
        require(sha256(out / it.key()) == it.value(), "Changed comparison output: " + (out / it.key()).string());
    return info;
}
inline void seal(const fs::path& out, const char* marker, const std::string& configuration, const std::vector<std::string>& files, Json info = Json::object()) {
    info["configuration_sha256"] = configuration; info["files"] = fileChecks(out, files); atomicJson(out / marker, info);
}
inline void metadata(TFile& file, const Json& info) {
    file.cd(); TObjString text(info.dump().c_str());
    require(text.Write("metadata_json") > 0, "Cannot save comparison metadata");
}
inline void closeRoot(TFile& file) {
    file.Close(); require(!file.TestBit(TFile::kWriteError), "Comparison ROOT close failed");
}
}
