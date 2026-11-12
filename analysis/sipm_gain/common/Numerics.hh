#pragma once
#include "GainMethods.hh"
#include <functional>
#include <utility>
namespace gain::detail {
using Bounds=std::vector<std::pair<double,double>>;
struct Minimum { Vec p,error;bool converged=false;int status=-1,covariance=0; };
struct Spectrum {double fftGain=NaN,acfGain=NaN,fftScore=0,acfScore=0;Vec candidates,freq,power,lags,acf,env;};
struct Fit {Json values;Vec mu;};
double sum(const Vec&); double square(double); Vec slice(const Vec&,int,int); int maximum(const Vec&);
double number(const Json&,const std::string&,double fallback=NaN); std::string join(const std::vector<std::string>&);
Vec smooth(const Vec&,double); double cdf(double); double deviance(const Vec&,const Vec&);
Minimum minimize(const std::function<double(const Vec&)>&,Vec,const Bounds&,bool,int);
double interpolate(double,const Vec&,const Vec&);
std::pair<Vec,Json> envelope(const Vec&,const Vec&,const std::string&,const Config&,bool);
Spectrum spectral(const Vec&,const Vec&,const Vec&,const Config&);
void autocorrelation(Spectrum&,const Vec&,const Vec&,double,const Config&);
Fit comb(const Vec&,const Vec&,const Vec&,const Config&,const Vec&,bool,bool);
Fit multiGaussian(const Vec&,const Vec&,const Vec&,const Config&,const Vec&,bool);
Fit envelopeComb(const Vec&,const Vec&,const Vec&,const Config&,const Vec&,bool);
}
