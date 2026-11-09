#pragma once
#include "../../analysis/impact_studies/StudyCommon.hh"
#include <TBranch.h>
#include <set>

// Reproduce the effective constants of the production Calibration.cxx and
// Extract.cxx revised on 2026-10-10. Keep historical parameters() unchanged.
// The returned map covers the 6,300 connected physical channels in layers 0-29.
// Production accepts layer 30/31 records, but they do not contribute to the
// fallback means and cannot affect the channels returned here.
inline std::map<int, Parameter> currentParameters(const std::string& hl,
                                                 const std::string& mip,
                                                 const std::string& ped) {
    std::map<int, Parameter> a;
    for (int l = 0; l < 30; ++l)
        for (int p = 0; p < 6; ++p)
            for (int c = 0; c < (p == 5 ? 30 : 36); ++c)
                a.emplace(l * 100000 + p * 10000 + c, Parameter{});
    std::set<int> validMip, validHL;
    auto productionPhysical = [](int id) {
        if (id < 0) return -1;
        const int l = id / 100000, p = id % 100000 / 10000, c = id % 100;
        if (l >= 32 || p >= 6 || c >= 36 || (p == 5 && c > 29)) return -1;
        return l * 100000 + p * 10000 + c;
    };
    {
        auto f = openFile(mip);
        auto* t = object<TTree>(*f, "MIP_Fit");
        int id = 0, ndf = 0;
        double mpv = 0., chi = 0.;
        bind(t, "CellID", &id);
        bind(t, "LandauMPV", &mpv);
        bind(t, "ChiSquare", &chi);
        bind(t, "NDF", &ndf);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            require(t->GetEntry(i) > 0, "Cannot read MIP factor entry");
            const int key = productionPhysical(id);
            if (key < 0) continue;
            const double limit = group(key) == 1 ? 2. : 1.7;
            if (!std::isfinite(mpv) || mpv <= 0. || ndf <= 0 ||
                !std::isfinite(chi) || chi < 0. || chi / ndf > limit) continue;
            if (a.count(key)) {
                a.at(key).mip = mpv;
                validMip.insert(key);
            }
        }
        t->ResetBranchAddresses();
    }
    {
        auto f = openFile(hl);
        auto* t = object<TTree>(*f, "InterCalib");
        int id = 0;
        double slope = 0., intercept = 0.;
        bind(t, "CellID", &id);
        bind(t, "Slope", &slope);
        bind(t, "Intercept", &intercept);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            require(t->GetEntry(i) > 0, "Cannot read HL factor entry");
            const int key = productionPhysical(id);
            if (key < 0 || !std::isfinite(slope) || slope < .02 || slope > .05 ||
                !std::isfinite(intercept)) continue;
            const double gain = 1. / slope, offset = -intercept / slope;
            if (!std::isfinite(gain) || !std::isfinite(offset)) continue;
            if (a.count(key)) {
                auto& p = a.at(key);
                p.gain = gain;
                p.intercept = offset;
                validHL.insert(key);
            }
        }
        t->ResetBranchAddresses();
    }
    long double mipSum[2] = {}, gainSum[2] = {}, offsetSum[2] = {};
    int mipN[2] = {}, hlN[2] = {};
    // Unique physical channels and ascending layer/chip/channel order reproduce
    // production's long-double sums, including the treatment of duplicates.
    for (const auto& [id, p] : a) {
        const int g = group(id);
        if (validMip.count(id)) { mipSum[g] += p.mip; ++mipN[g]; }
        if (validHL.count(id)) {
            gainSum[g] += p.gain;
            offsetSum[g] += p.intercept;
            ++hlN[g];
        }
    }
    double mipMean[2], gainMean[2], offsetMean[2];
    for (int g = 0; g < 2; ++g) {
        require(mipN[g] > 0 && hlN[g] > 0, "No valid constants in a SiPM group");
        mipMean[g] = static_cast<double>(mipSum[g] / mipN[g]);
        gainMean[g] = static_cast<double>(gainSum[g] / hlN[g]);
        offsetMean[g] = static_cast<double>(offsetSum[g] / hlN[g]);
        require(std::isfinite(mipMean[g]) && mipMean[g] > 0. &&
                std::isfinite(gainMean[g]) && gainMean[g] > 0. &&
                std::isfinite(offsetMean[g]), "Invalid fallback constants");
    }
    for (auto& [id, p] : a) {
        const int g = group(id);
        // A valid MPV=1 is a measured value, not the old sentinel convention.
        p.mipFallback = validMip.count(id) == 0;
        p.hlFallback = validHL.count(id) == 0;
        if (p.mipFallback) p.mip = mipMean[g];
        if (p.hlFallback) { p.gain = gainMean[g]; p.intercept = offsetMean[g]; }
    }
    {
        auto f = openFile(ped);
        auto* t = object<TTree>(*f, "ChnLevel");
        std::vector<int>* ids = nullptr;
        std::vector<float>* hf = nullptr;
        std::vector<float>* lf = nullptr;
        std::vector<double>* hd = nullptr;
        std::vector<double>* ld = nullptr;
        bind(t, "CellID", &ids);
        auto* hb = t->GetBranch("PedHighMean");
        auto* lb = t->GetBranch("PedLowMean");
        require(hb && lb, "Missing pedestal mean branches");
        const std::string highType = hb->GetClassName(), lowType = lb->GetClassName();
        const bool highFloat = highType.find("vector<float") != std::string::npos;
        const bool lowFloat = lowType.find("vector<float") != std::string::npos;
        require(highFloat || highType.find("vector<double") != std::string::npos,
                "Unsupported HG pedestal branch type");
        require(lowFloat || lowType.find("vector<double") != std::string::npos,
                "Unsupported LG pedestal branch type");
        if (highFloat) bind(t, "PedHighMean", &hf); else bind(t, "PedHighMean", &hd);
        if (lowFloat) bind(t, "PedLowMean", &lf); else bind(t, "PedLowMean", &ld);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            require(t->GetEntry(i) > 0, "Cannot read pedestal entry");
            require(ids && (highFloat ? hf != nullptr : hd != nullptr) &&
                    (lowFloat ? lf != nullptr : ld != nullptr), "Null pedestal vector");
            require(ids->size() == (highFloat ? hf->size() : hd->size()) &&
                    ids->size() == (lowFloat ? lf->size() : ld->size()),
                    "Pedestal size mismatch");
            for (size_t j = 0; j < ids->size(); ++j) {
                const int key = productionPhysical(ids->at(j));
                if (!a.count(key)) continue;
                auto& p = a.at(key);
                // Calibration.cxx reads vector<float>, including float rounding
                // when the stored parameter vector is double precision.
                p.ph = static_cast<float>(highFloat ? hf->at(j) : hd->at(j));
                p.pl = static_cast<float>(lowFloat ? lf->at(j) : ld->at(j));
                p.ped = true;
            }
        }
        t->ResetBranchAddresses();
        delete ids; delete hf; delete lf; delete hd; delete ld;
        // Missing pedestal entries retain zero, matching current production;
        // p.ped remains false so the caller can report them explicitly.
    }
    return a;
}
