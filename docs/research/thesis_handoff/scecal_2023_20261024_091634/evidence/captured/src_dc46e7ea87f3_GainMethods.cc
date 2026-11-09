#include "Numerics.hh"
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace gain {
using namespace detail;
Json Config::json()const{return {{"gain_min",gain_min},{"gain_max",gain_max},{"adc_min",adc_min},{"adc_max",adc_max},{"max_span",max_span},
    {"min_entries",min_entries},{"min_cycles",min_cycles},{"null_replicates",null_replicates},{"bootstrap",bootstrap},
    {"alpha",alpha},{"stability_fraction",stability_fraction},{"seed",seed}};}
void Config::validate()const{
    if(!(0<gain_min&&gain_min<gain_max&&0<=adc_min&&adc_min<adc_max&&adc_max<=4095&&max_span>=3*gain_min&&min_entries>0&&min_cycles>0&&
        null_replicates>0&&bootstrap==0&&alpha>0&&alpha<1&&stability_fraction>0&&seed>0))throw std::runtime_error("Invalid fit configuration (bootstrap must be 0 in this workflow)");
}
namespace {
struct Estimate { Fit fit; Spectrum spectral; Json envelopeInfo; };
Estimate estimate(const Vec& x,const Vec& y,const std::string& kind,const Config& c,int method,bool fast=false,double seed=NaN){
    auto [env,info]=envelope(x,y,kind,c,fast);
    auto sp=spectral(x,y,env,c);
    if(kind=="led"){
        env=smooth(y,std::max(2.,sp.fftGain*.6/(x[1]-x[0])));
        for(double& v:env)v=std::max(v,.01);
        sp=spectral(x,y,env,c);
    }
    Fit fit;
    Vec candidates=std::isfinite(seed)?Vec{seed}:sp.candidates;
    if(method==0)fit=multiGaussian(x,y,env,c,candidates,fast);
    else if(method==1)fit=envelopeComb(x,y,env,c,candidates,fast);
    else fit.values={{"gain",method==2?sp.fftGain:sp.acfGain},{"score",method==2?sp.fftScore:sp.acfScore}};
    return {fit,sp,info};
}
}
Result analyze(const Vec& counts,const std::string& kind,const Config& c,const std::string& method){
    c.validate();
    if(kind!="led"&&kind!="electron"&&kind!="muon")throw std::runtime_error("Unknown spectrum kind");
    auto found=std::find(Methods.begin(),Methods.end(),method);
    if(found==Methods.end())throw std::runtime_error("Unknown estimator: "+method);
    int index=found-Methods.begin();
    if(counts.size()!=4096)throw std::runtime_error("Expected 4096 raw ADC bins");
    int low=std::ceil(c.adc_min),high=std::floor(c.adc_max);double total=0;
    for(int i=low;i<=high;++i)total+=counts[i];
    int lo=-1,hi=-1;double cumulative=0;
    for(int i=low;i<=high&&total>0;++i){cumulative+=counts[i];if(lo<0&&cumulative/total>=.005)lo=i;if(hi<0&&cumulative/total>=.995)hi=i;}
    std::vector<int> roi;
    if(lo>=0){if(hi<lo)hi=high;hi=std::min(hi,lo+c.max_span-1);for(int i=lo;i<=hi;++i)roi.push_back(i);}
    Result output;auto& report=output.report;
    report={{"config",c.json()},{"roi",roi},{"entries",sum(counts)},{"method",method},{"uncertainty_status","not_evaluated"},{"bootstrap_requested",0}};
    Json result={{"gain",NaN},{"stat_error",NaN},{"fit_error",NaN},{"bootstrap_requested",0},{"bootstrap_success",0},{"uncertainty_status","not_evaluated"}};
    Vec x,y;for(int i:roi){x.push_back(i);y.push_back(counts[i]);}
    if(x.size()<c.min_cycles*c.gain_min||sum(y)<c.min_entries){
        result["status"]="insufficient_data";result["reasons"]="insufficient_entries_or_adc_span";
        report["status"]="insufficient_data";report["methods"][method]=result;return output;
    }
    auto nominal=estimate(x,y,kind,c,index);result.update(nominal.fit.values);
    double g=number(result,"gain");TRandom3 rng(c.seed);int exceed=0;
    // This lightweight smooth-null diagnostic performs FFT/ACF only; no gain refits.
    for(int iter=0;iter<c.null_replicates;++iter){
        Vec fake;for(double e:nominal.spectral.env)fake.push_back(rng.PoissonD(std::max(e,0.)));
        auto null=spectral(x,fake,nominal.spectral.env,c);
        exceed+=index==3?null.acfScore>=nominal.spectral.acfScore:null.fftScore>=nominal.spectral.fftScore;
    }
    auto& d=output.detail;d.x=x;d.y=y;d.env=nominal.spectral.env;
    d.freq=nominal.spectral.freq;d.power=nominal.spectral.power;d.lags=nominal.spectral.lags;d.acf=nominal.spectral.acf;
    if(index<2)d.curves[index]=nominal.fit.mu;
    Vec variations;
    auto vary=[&](const Vec& xx,const Vec& yy){auto v=estimate(xx,yy,kind,c,index,true,g);variations.push_back(v.fit.values.value("converged",true)?number(v.fit.values,"gain"):NaN);};
    if(x.size()>c.min_cycles*c.gain_max){int trim=std::max(1,int(.08*x.size()));vary(slice(x,trim,x.size()),slice(y,trim,y.size()));vary(slice(x,0,x.size()-trim),slice(y,0,y.size()-trim));}
    if(c.gain_min>=5&&x.size()>30){Vec xx,yy;for(size_t i=0;i+1<x.size();i+=2){xx.push_back((x[i]+x[i+1])/2);yy.push_back(y[i]+y[i+1]);}vary(xx,yy);}
    if(index>=2){auto env=smooth(y,std::max(2.,c.gain_max*.6));for(double& v:env)v=std::max(v,.01);auto sp=spectral(x,y,env,c);variations.push_back(index==2?sp.fftGain:sp.acfGain);}
    double stability=NaN;int success=0;
    for(double v:variations)if(std::isfinite(v)){stability=success?std::max(stability,std::abs(v-g)):std::abs(v-g);++success;}
    double p=(exceed+1.)/(c.null_replicates+1.);
    result.update({{"null_p",p},{"stability_ADC",stability},{"variants_success",success},{"variants_total",variations.size()}});
    std::vector<std::string> reasons;
    if(p>c.alpha)reasons.push_back("no_significant_periodicity");
    if(!(c.gain_min<g&&g<c.gain_max))reasons.push_back("gain_boundary");
    if((x.back()-x.front()+1)/g<c.min_cycles)reasons.push_back("too_few_periods");
    if(index<2){
        if(!result.value("converged",false))reasons.push_back("fit_failed");
        if(number(result,"covariance_status",0)<2)reasons.push_back("poor_covariance");
        if(result.value("boundary_gain",true))reasons.push_back("fit_boundary");
        if(number(result,"alternative_delta_bic",INFINITY)<6)reasons.push_back("ambiguous_period");
        if(number(result,"ndf",0)<=0||number(result,"deviance")/std::max(1.,number(result,"ndf",0))>3)reasons.push_back("poor_shape_fit");
        if(number(result,"width_fraction",1)>.38)reasons.push_back("unresolved_peaks");
        if(index==0&&number(result,"consecutive_peaks",0)<3)reasons.push_back("missing_adjacent_peaks");
    }
    if(!std::isfinite(stability)||success<int(variations.size()))reasons.push_back("stability_incomplete");
    else if(stability/g>c.stability_fraction)reasons.push_back("range_or_binning_unstable");
    result["status"]=reasons.empty()?"candidate":"review";result["reasons"]=join(reasons);
    report["methods"][method]=result;report["status"]=result["status"];report["envelope"]=nominal.envelopeInfo;report["roi_entries"]=sum(y);
    report["note"]="Bootstrap disabled. Candidate means shape/periodicity/stability checks passed; statistical uncertainty and method compatibility are not evaluated. fit_error is conditional Minuit error only.";
    return output;
}
}
