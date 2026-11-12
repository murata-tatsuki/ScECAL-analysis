#pragma once
#include <nlohmann/json.hpp>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace gain {
using Json = nlohmann::json;
using Vec = std::vector<double>;
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline const std::array<std::string,4> Methods = {"multi_gaussian","envelope_comb","fft","autocorrelation"};
struct Config {
    double gain_min=6., gain_max=40., adc_min=250., adc_max=1400.;
    int max_span=400, min_entries=120;
    double min_cycles=3.;
    int null_replicates=199, bootstrap=0;
    double alpha=.01, stability_fraction=.05;
    unsigned seed=20261013;
    Json json() const;
    void validate() const;
};
struct Detail {
    Vec x,y,env,freq,power,lags,acf;
    std::array<Vec,2> curves;
};
struct Result { Json report; Detail detail; };
// One independent estimator; no bootstrap/refitting for statistical errors.
Result analyze(const Vec& counts,const std::string& kind,const Config& config,const std::string& method);
}
