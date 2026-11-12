#pragma once
// Exact, process-local reuse of the current muon histogram-fit result.
// No comparison tolerance, bin merging, response scaling, or fit approximation.
#include <TH1D.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <list>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace mmr {
struct FitCacheStats {
    std::uint64_t requests = 0;
    std::uint64_t hits = 0;
    std::uint64_t misses = 0;
    std::uint64_t actualFits = 0;
    std::uint64_t evictions = 0;
    std::uint64_t oversizeSkips = 0;
    std::uint64_t bufferBypasses = 0;
};

namespace exact_fit_cache_detail {
template<class T> inline void append(std::string& key, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    key.append(reinterpret_cast<const char*>(&value), sizeof(value));
}
inline void appendString(std::string& key, const char* value) {
    const std::size_t size = std::strlen(value);
    append(key, size);
    key.append(value, size);
}

// This canonicalization mirrors fitHistogram's existing order exactly.
// Low-statistics and missing-threshold failures retain their raw bounds.
// For later outcomes, the reported and fitted bounds have already been clamped.
inline void effectiveBounds(const TH1D* h, double& low, double& high, double minHits) {
    if (!h || !std::isfinite(low) || !std::isfinite(high)) return;
    double sumw = 0, sumw2 = 0;
    for (int bin = 0; bin <= h->GetNbinsX() + 1; ++bin) {
        sumw += h->GetBinContent(bin);
        sumw2 += std::pow(h->GetBinError(bin), 2);
    }
    const double neff = sumw2 > 0 ? sumw * sumw / sumw2 : 0;
    if (neff < minHits) return;
    const double width = h->GetBinWidth(1);
    low = std::max(width * .5, low);
    high = std::min(high, h->GetXaxis()->GetXmax() - width * .5);
}

inline std::string fingerprint(const TH1D* h, double low, double high,
                               double scale, double normalization,
                               double minHits, double maxChi2) {
    std::string key;
    if (h) key.reserve(std::size_t(h->GetNcells()) * 80 + 256);
    const bool present = h != nullptr;
    append(key, present);
    effectiveBounds(h, low, high, minHits);
    for (double value : {low, high, scale, normalization, minHits, maxChi2})
        append(key, value);
    if (!h) return key;

    // The kernel sees these exact native coordinates, bins, weights, and errors.
    // Function names and histogram titles do not affect its custom fit.
    appendString(key, h->ClassName());
    append(key, h->GetDimension());
    append(key, h->GetNbinsX());
    append(key, h->GetNcells());
    append(key, h->GetEntries());
    append(key, h->GetNormFactor());
    append(key, h->GetBinErrorOption());
    append(key, h->GetStatOverflows());
    append(key, h->TestBit(TH1::kIsAverage));

    const auto* axis = h->GetXaxis();
    append(key, axis->GetXmin());
    append(key, axis->GetXmax());
    append(key, axis->GetFirst());
    append(key, axis->GetLast());
    append(key, axis->TestBit(TAxis::kAxisRange));
    append(key, axis->CanExtend());
    const auto* edges = axis->GetXbins();
    append(key, edges->GetSize());
    for (int i = 0; i < edges->GetSize(); ++i) append(key, edges->At(i));
    // Include actual coordinate calculations, including both flow bins.
    for (int bin = 0; bin <= h->GetNbinsX() + 1; ++bin) {
        append(key, axis->GetBinLowEdge(bin));
        append(key, axis->GetBinUpEdge(bin));
        append(key, axis->GetBinCenter(bin));
        append(key, axis->GetBinWidth(bin));
    }
    for (int bin = 0; bin < h->GetNcells(); ++bin) {
        append(key, h->GetBinContent(bin));
        append(key, h->GetBinError(bin));
        append(key, h->GetBinErrorLow(bin));
        append(key, h->GetBinErrorUp(bin));
    }
    const auto* errors = h->GetSumw2();
    append(key, errors->GetSize());
    for (int bin = 0; bin < errors->GetSize(); ++bin) append(key, errors->At(bin));
    return key;
}
} // namespace exact_fit_cache_detail

// Result must be copyable and expose rootStatus; the current Fit uses -999 for
// returns that did not enter ROOT's minimizer. Each cache must use the same
// fitting kernel and ROOT minimizer configuration throughout its lifetime.
// The native fit signature is included in full, including normalization.
//
// get(..., callback) executes callback with the ORIGINAL raw arguments on a
// miss; key-only bound canonicalization cannot change failure metadata.
// Instances are owned by individual processes; this class is not thread-safe.
template<class Result> class ExactFitCache {
    struct Entry {
        std::string key;
        Result result;
        std::size_t bytes;
    };
    using Entries = std::list<Entry>;
    using Iterator = typename Entries::iterator;
    // string_view points into stable list storage. Hash collisions are followed
    // by full byte equality; a hash value alone never identifies a fit.
    using Index = std::unordered_map<std::string_view, Iterator>;
    std::size_t budget_, used_ = 0;
    Entries entries_;
    Index index_;
    FitCacheStats stats_;

    template<class Fitter> Result calculate(Fitter&& fitter) {
        Result result = std::forward<Fitter>(fitter)();
        if (result.rootStatus != -999) ++stats_.actualFits;
        return result;
    }
    void evictOne() {
        auto last = std::prev(entries_.end());
        index_.erase(std::string_view(last->key));
        used_ -= last->bytes;
        entries_.erase(last);
        ++stats_.evictions;
    }
public:
    explicit ExactFitCache(std::size_t budgetBytes = 32u * 1024u * 1024u)
        : budget_(budgetBytes) {}
    ExactFitCache(const ExactFitCache&) = delete;
    ExactFitCache& operator=(const ExactFitCache&) = delete;
    ExactFitCache(ExactFitCache&&) = delete;
    ExactFitCache& operator=(ExactFitCache&&) = delete;

    template<class Fitter>
    Result get(const TH1D* h, double low, double high, double scale,
               double normalization, double minHits, double maxChi2,
               Fitter&& fitter) {
        ++stats_.requests;
        if (!budget_ || (h && h->GetBufferLength() != 0)) {
            ++stats_.misses;
            if (h && h->GetBufferLength() != 0) ++stats_.bufferBypasses;
            return calculate(std::forward<Fitter>(fitter));
        }
        std::string key = exact_fit_cache_detail::fingerprint(
            h, low, high, scale, normalization, minHits, maxChi2);
        auto found = index_.find(std::string_view(key));
        if (found != index_.end()) {
            ++stats_.hits;
            entries_.splice(entries_.begin(), entries_, found->second);
            return found->second->result;
        }
        ++stats_.misses;
        Result result = calculate(std::forward<Fitter>(fitter));
        // Charge fixed overhead conservatively, in addition to the stored key
        // and result. The budget applies to retained entries, not an input key.
        const std::size_t bytes = sizeof(Entry) + sizeof(typename Index::value_type)
                                + key.capacity() + 1 + 256;
        if (bytes > budget_) {
            ++stats_.oversizeSkips;
            return result;
        }
        while (used_ > budget_ - bytes) evictOne();
        entries_.push_front(Entry{std::move(key), result, bytes});
        try {
            index_.emplace(std::string_view(entries_.front().key), entries_.begin());
        } catch (...) {
            entries_.pop_front();
            throw;
        }
        used_ += bytes;
        return result;
    }
    const FitCacheStats& stats() const { return stats_; }
    std::size_t size() const { return entries_.size(); }
    std::size_t bytes() const { return used_; }
    std::size_t budgetBytes() const { return budget_; }
    void clear() {
        index_.clear();
        index_.rehash(0);
        entries_.clear();
        used_ = 0;
    }
    void reset() {
        clear();
        stats_ = {};
    }
};
} // namespace mmr

