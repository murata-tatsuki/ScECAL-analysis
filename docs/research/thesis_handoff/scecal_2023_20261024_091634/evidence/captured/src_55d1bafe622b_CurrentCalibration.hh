#pragma once
#include "MuonMIPCommon.hh"
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace mmr {

struct Parameter {
    double mip = 1., ph = 0., pl = 0., gain = 1., offset = 0.;
    double rawMip = std::numeric_limits<double>::quiet_NaN();
    double chi2ndf = std::numeric_limits<double>::quiet_NaN();
    double threshold = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    bool mipFallback = true, pedPresent = false, hlFallback = true;
    bool thresholdFallback = true;
};

namespace current_calibration_detail {
inline TTree* getTree(TFile& file, const char* name) {
    auto* t = dynamic_cast<TTree*>(file.Get(name));
    require(t != nullptr, std::string("Missing ") + name + " in " + file.GetName());
    t->SetBranchStatus("*", 0);
    return t;
}
template<class T> inline void setBranch(TTree* t, const char* name, T* ptr) {
    t->SetBranchStatus(name, 1);
    require(t->SetBranchAddress(name, ptr) >= 0, std::string("Invalid branch ") + name);
}
inline int productionPhysical(int id) {
    if (id < 0) return -1;
    const int l = id / 100000, p = id % 100000 / 10000, c = id % 100;
    if (l >= 32 || p >= 6 || c >= 36 || (p == 5 && c > 29)) return -1;
    return l * 100000 + p * 10000 + c;
}
inline std::map<int, Parameter> emptyParameters() {
    std::map<int, Parameter> result;
    for (int l = 0; l < 30; ++l)
        for (int p = 0; p < 6; ++p)
            for (int c = 0; c < (p == 5 ? 30 : 36); ++c)
                result.emplace(l * 100000 + p * 10000 + c, Parameter{});
    return result;
}
} // namespace current_calibration_detail

// Current (2026-10-10) production MIP policy. Invalid fits never contribute to
// individual constants or fallback means. Each accepted physical channel is
// counted once; duplicate accepted records use the last accepted value.
inline std::map<int, Parameter> loadMip(const std::string& path, bool allow = false) {
    using namespace current_calibration_detail;
    auto result = emptyParameters();
    auto file = openFile(path, allow);
    auto* t = getTree(*file, "MIP_Fit");
    int id = 0, ndf = 0;
    double mpv = 0., chi = 0.;
    setBranch(t, "CellID", &id);
    setBranch(t, "LandauMPV", &mpv);
    setBranch(t, "ChiSquare", &chi);
    setBranch(t, "NDF", &ndf);
    std::set<int> accepted;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        require(t->GetEntry(i) > 0, "Cannot read MIP parameter entry");
        const int key = productionPhysical(id);
        if (!result.count(key)) continue;
        auto& p = result.at(key);
        p.rawMip = mpv;
        p.chi2ndf = ndf > 0 && std::isfinite(chi)
            ? chi / ndf : std::numeric_limits<double>::quiet_NaN();
        const double limit = group(key) == 1 ? 2. : 1.7;
        if (!std::isfinite(mpv) || mpv <= 0. || ndf <= 0 ||
            !std::isfinite(chi) || chi < 0. || chi / ndf > limit) continue;
        p.mip = mpv;
        accepted.insert(key);
    }
    t->ResetBranchAddresses();
    long double sum[2] = {};
    int count[2] = {};
    for (const auto& [id, p] : result) if (accepted.count(id)) {
        sum[group(id)] += p.mip;
        ++count[group(id)];
    }
    double mean[2];
    for (int g = 0; g < 2; ++g) {
        require(count[g] > 0, "No valid MIP donors in a SiPM group");
        mean[g] = static_cast<double>(sum[g] / count[g]);
        require(std::isfinite(mean[g]) && mean[g] > 0., "Invalid fallback MIP");
    }
    for (auto& [id, p] : result) {
        p.mipFallback = accepted.count(id) == 0;
        if (p.mipFallback) p.mip = mean[group(id)];
    }
    return result;
}

inline std::map<int, Parameter> loadParameters(const std::string& pedestal,
                                               const std::string& mip,
                                               const std::string& hl,
                                               const std::string& threshold,
                                               bool allow = false) {
    using namespace current_calibration_detail;
    auto result = loadMip(mip, allow);
    {
        auto file = openFile(hl, allow);
        auto* t = getTree(*file, "InterCalib");
        int id = 0;
        double slope = 0., intercept = 0.;
        setBranch(t, "CellID", &id);
        setBranch(t, "Slope", &slope);
        setBranch(t, "Intercept", &intercept);
        std::set<int> accepted;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            require(t->GetEntry(i) > 0, "Cannot read HG/LG parameter entry");
            const int key = productionPhysical(id);
            if (!result.count(key) || !std::isfinite(slope) ||
                slope < .02 || slope > .05 || !std::isfinite(intercept)) continue;
            const double gain = 1. / slope, offset = -intercept / slope;
            if (!std::isfinite(gain) || !std::isfinite(offset)) continue;
            auto& p = result.at(key);
            p.gain = gain; p.offset = offset;
            accepted.insert(key);
        }
        t->ResetBranchAddresses();
        long double gainSum[2] = {}, offsetSum[2] = {};
        int count[2] = {};
        for (const auto& [id, p] : result) if (accepted.count(id)) {
            gainSum[group(id)] += p.gain;
            offsetSum[group(id)] += p.offset;
            ++count[group(id)];
        }
        double gainMean[2], offsetMean[2];
        for (int g = 0; g < 2; ++g) {
            require(count[g] > 0, "No valid HG/LG donors in a SiPM group");
            gainMean[g] = static_cast<double>(gainSum[g] / count[g]);
            offsetMean[g] = static_cast<double>(offsetSum[g] / count[g]);
            require(std::isfinite(gainMean[g]) && gainMean[g] > 0. &&
                    std::isfinite(offsetMean[g]), "Invalid fallback HG/LG constants");
        }
        for (auto& [id, p] : result) {
            p.hlFallback = accepted.count(id) == 0;
            if (p.hlFallback) {
                p.gain = gainMean[group(id)]; p.offset = offsetMean[group(id)];
            }
        }
    }
    {
        auto file = openFile(pedestal, allow);
        auto* t = getTree(*file, "ChnLevel");
        std::vector<int>* ids = nullptr;
        std::vector<float> *hf = nullptr, *lf = nullptr;
        std::vector<double> *hd = nullptr, *ld = nullptr;
        setBranch(t, "CellID", &ids);
        auto* highBranch = t->GetBranch("PedHighMean");
        auto* lowBranch = t->GetBranch("PedLowMean");
        require(highBranch && lowBranch, "Missing pedestal mean branches");
        const std::string highType = highBranch->GetClassName();
        const std::string lowType = lowBranch->GetClassName();
        const bool highFloat = highType.find("vector<float") != std::string::npos;
        const bool lowFloat = lowType.find("vector<float") != std::string::npos;
        require(highFloat || highType.find("vector<double") != std::string::npos,
                "Unsupported HG pedestal vector type");
        require(lowFloat || lowType.find("vector<double") != std::string::npos,
                "Unsupported LG pedestal vector type");
        if (highFloat) setBranch(t, "PedHighMean", &hf);
        else setBranch(t, "PedHighMean", &hd);
        if (lowFloat) setBranch(t, "PedLowMean", &lf);
        else setBranch(t, "PedLowMean", &ld);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            require(t->GetEntry(i) > 0, "Cannot read pedestal parameter entry");
            require(ids && (highFloat ? hf != nullptr : hd != nullptr) &&
                    (lowFloat ? lf != nullptr : ld != nullptr), "Null pedestal vector");
            require(ids->size() == (highFloat ? hf->size() : hd->size()) &&
                    ids->size() == (lowFloat ? lf->size() : ld->size()),
                    "Pedestal vector size mismatch");
            for (size_t j = 0; j < ids->size(); ++j) {
                const int key = productionPhysical(ids->at(j));
                if (!result.count(key)) continue;
                auto& p = result.at(key);
                // Production Calibration reads vector<float> even when the
                // stored branch is double precision. Preserve that rounding.
                p.ph = static_cast<float>(highFloat ? hf->at(j) : hd->at(j));
                p.pl = static_cast<float>(lowFloat ? lf->at(j) : ld->at(j));
                p.pedPresent = true;
            }
        }
        t->ResetBranchAddresses();
        delete ids; delete hf; delete lf; delete hd; delete ld;
        // Missing entries retain production's zero, and pedPresent=false.
    }
    // Thresholds are diagnostics of the trigger plateau, not extra hit cuts.
    // Unknown channel/chip thresholds stay NaN rather than becoming zero.
    if (threshold.empty() || threshold == "-") return result;
    auto file = openFile(threshold, allow);
    auto* t = dynamic_cast<TTree*>(file->Get("expErfThre"));
    if (!t) return result;
    t->SetBranchStatus("*", 0);
    int id = 0;
    double value = 0., width = 0.;
    setBranch(t, "CellID", &id);
    setBranch(t, "threshold", &value);
    setBranch(t, "sigma", &width);
    double chipSum[30][6] = {}, chipSigmaSum[30][6] = {};
    int chipCount[30][6] = {};
    std::set<int> present;
    for (auto& [key, p] : result) { p.threshold = 1.; p.sigma = 0.; }
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        require(t->GetEntry(i) > 0, "Cannot read threshold parameter entry");
        if (id < 0) continue;
        const int l = id / 100000, c = id % 100000 / 10000, ch = id % 100;
        if (l >= 30 || c >= 6 || ch >= 36) continue;
        // Production accumulates source rows per chip (also spare channels),
        // rather than averaging only the individually accepted MIP channels.
        chipSum[l][c] += value;
        chipSigmaSum[l][c] += width;
        ++chipCount[l][c];
        const int key = l * 100000 + c * 10000 + ch;
        if (!result.count(key)) continue;
        auto& p = result.at(key);
        p.threshold = value; p.sigma = width;
        present.insert(key);
    }
    t->ResetBranchAddresses();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (auto& [key, p] : result) {
        const int l = key / 100000, c = key % 100000 / 10000;
        const bool replaceThreshold = p.threshold == 1.;
        const bool replaceSigma = p.sigma == 0.;
        p.thresholdFallback = !present.count(key) || replaceThreshold || replaceSigma;
        if (replaceThreshold)
            p.threshold = chipCount[l][c] ? chipSum[l][c] / chipCount[l][c] : nan;
        if (replaceSigma)
            p.sigma = chipCount[l][c] ? chipSigmaSum[l][c] / chipCount[l][c] : nan;
        if (!std::isfinite(p.threshold) || !std::isfinite(p.sigma) || p.sigma <= 0.) {
            p.threshold = nan; p.sigma = nan; p.thresholdFallback = true;
        }
    }
    return result;
}

} // namespace mmr
