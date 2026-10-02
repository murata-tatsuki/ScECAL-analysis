#pragma once
// Event/CellID/pedestal helpers adapted from adc_energy/ChannelResponse.cc.
#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using EventKey = std::tuple<int, int, int>;
using FilePtr = std::unique_ptr<TFile>;

FilePtr openFile(const std::string& path, const char* mode = "READ") {
    FilePtr f(TFile::Open(path.c_str(), mode));
    if (!f || f->IsZombie()) throw std::runtime_error("Cannot open ROOT file: " + path);
    return f;
}
TTree* tree(TFile& f, const char* name) {
    auto* t = dynamic_cast<TTree*>(f.Get(name));
    if (!t) throw std::runtime_error(std::string("Missing tree ") + name + " in " + f.GetName());
    t->SetBranchStatus("*", false);
    return t;
}
template<class T> void bind(TTree* t, const char* name, T* value) {
    if (!t->GetBranch(name)) throw std::runtime_error(std::string("Missing branch: ") + name);
    t->SetBranchStatus(name, true);
    if (t->SetBranchAddress(name, value) < 0)
        throw std::runtime_error(std::string("Incompatible branch type: ") + name);
}
int physicalID(int c) {
    int layer = c / 100000, chip = c % 100000 / 10000, channel = c % 100;
    if (c < 0 || layer >= 32 || chip >= 6 || channel >= 36)
        throw std::runtime_error("Invalid CellID: " + std::to_string(c));
    return layer * 100000 + chip * 10000 + channel;
}
bool connected(int c) { return c % 100000 / 10000 != 5 || c % 100 < 30; }
int channelSlot(int physical) {
    return (physical / 100000 * 6 + physical % 100000 / 10000) * 36 + physical % 100;
}
size_t hitSlot(int full) {
    const int p = physicalID(full);
    return ((p / 100000 * 6 + p % 100000 / 10000) * 100 + full % 10000 / 100) * 36 + p % 100;
}
struct HitIndex {
    uint64_t calibrated = 0, raw = 0, signal = 0;
    size_t index = 0;
};
void cacheActiveBranches(TTree* t) {
    t->SetCacheSize(0);
    t->SetCacheSize(32 * 1024 * 1024);
    TIter next(t->GetListOfBranches());
    while (auto* branch = next()) if (t->GetBranchStatus(branch->GetName()))
        t->AddBranchToCache(branch->GetName(), true);
    t->StopCacheLearningPhase();
}

std::set<int> integerList(const std::string& s, bool ranges) {
    std::set<int> result;
    std::istringstream input(s);
    for (std::string item; std::getline(input, item, ',');) {
        auto dash = item.find('-');
        int lo = std::stoi(item), hi = lo;
        if (dash != std::string::npos) {
            if (!ranges) throw std::runtime_error("Channel ranges are not supported");
            hi = std::stoi(item.substr(dash + 1));
        }
        if (lo < 0 || hi < lo || (ranges && hi >= 32)) throw std::runtime_error("Invalid selection: " + s);
        for (int i = lo; i <= hi; ++i) result.insert(i);
    }
    if (result.empty()) throw std::runtime_error("Empty selection");
    return result;
}
struct Pedestal { double hg, lg; };
template<class Real> std::map<int, Pedestal> readPedestalTree(TTree* t) {
    std::vector<int>* ids = nullptr;
    std::vector<Real>* h = nullptr;
    std::vector<Real>* l = nullptr;
    bind(t, "CellID", &ids); bind(t, "PedHighMean", &h); bind(t, "PedLowMean", &l);
    std::map<int, Pedestal> result;
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
        if (t->GetEntry(e) <= 0) throw std::runtime_error("Cannot read pedestal entry");
        if (ids->size() != h->size() || ids->size() != l->size()) throw std::runtime_error("Pedestal vector length mismatch");
        for (size_t j = 0; j < ids->size(); ++j) {
            int id = physicalID(ids->at(j));
            Pedestal value{double(h->at(j)), double(l->at(j))};
            if (!std::isfinite(value.hg) || !std::isfinite(value.lg)) throw std::runtime_error("Nonfinite pedestal");
            auto [it, fresh] = result.emplace(id, value);
            if (!fresh && (std::abs(it->second.hg - value.hg) > 1e-6 || std::abs(it->second.lg - value.lg) > 1e-6))
                throw std::runtime_error("Conflicting pedestal values for channel " + std::to_string(id));
        }
    }
    t->ResetBranchAddresses();
    delete ids; delete h; delete l;
    return result;
}
std::map<int, Pedestal> readPedestal(const std::string& name) {
    auto f = openFile(name); auto* t = tree(*f, "ChnLevel");
    auto* branch = t->GetBranch("PedHighMean");
    if (!branch) throw std::runtime_error("Missing PedHighMean");
    std::string type = branch->GetClassName();
    if (type.find("double") != std::string::npos) return readPedestalTree<double>(t);
    if (type.find("float") != std::string::npos) return readPedestalTree<float>(t);
    throw std::runtime_error("Unsupported pedestal branch type: " + type);
}

