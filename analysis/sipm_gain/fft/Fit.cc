#include "Numerics.hh"
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVirtualFFT.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace gain::detail {
Spectrum spectral(const Vec& x,const Vec& y,const Vec& env,const Config& c){
    int n=x.size(),nfft=16;while(nfft<8*n)nfft*=2;double dx=x[1]-x[0];
    Vec z(nfft,0),window(n);double weight=0;
    for(int i=0;i<n;++i){window[i]=.5-.5*std::cos(2*TMath::Pi()*i/(n-1));weight+=square(window[i]);z[i]=(y[i]-env[i])/std::sqrt(std::max(env[i],1.))*window[i];}
    // Keep private plans, rather than ROOT's shared current transform.
    static std::map<int,std::unique_ptr<TVirtualFFT>> plans;
    auto& fft=plans[nfft];if(!fft)fft.reset(TVirtualFFT::FFT(1,&nfft,"R2C ES K"));
    if(!fft)throw std::runtime_error("ROOT FFTW backend is unavailable");
    fft->SetPoints(z.data());fft->Transform();Vec power(nfft/2+1);std::vector<int> band,peaks;
    Spectrum r;r.env=env;
    for(int k=0;k<=nfft/2;++k){double re,im;fft->GetPointComplex(k,re,im);power[k]=(re*re+im*im)/std::max(weight,1.);
        double freq=k/(nfft*dx);if(freq>=1/c.gain_max&&freq<=1/c.gain_min){band.push_back(k);r.freq.push_back(freq);r.power.push_back(power[k]);}}
    if(band.empty())throw std::runtime_error("Empty FFT search band");
    for(int k:band)if(power[k]>=power[std::max(0,k-1)]&&power[k]>=power[std::min(nfft/2,k+1)])peaks.push_back(k);
    if(peaks.empty())peaks.push_back(*std::max_element(band.begin(),band.end(),[&](int a,int b){return power[a]<power[b];}));
    std::stable_sort(peaks.begin(),peaks.end(),[&](int a,int b){return power[a]>power[b];});
    for(int k:peaks){double f=k/(nfft*dx);if(k>0&&k<nfft/2){double a=std::log(std::max(power[k-1],1e-30)),b=std::log(std::max(power[k],1e-30)),d=std::log(std::max(power[k+1],1e-30));
        double den=a-2*b+d;if(den<0)f+=std::clamp(.5*(a-d)/den,-.5,.5)/(nfft*dx);}
        double g=1/f;bool distinct=true;for(double v:r.candidates)distinct &= std::abs(g/v-1)>.025;
        if(distinct)r.candidates.push_back(g);
        if(r.candidates.size()==3)break;
    }
    r.fftGain=r.candidates.front();r.fftScore=*std::max_element(r.power.begin(),r.power.end());
    autocorrelation(r,z,window,dx,c);
    return r;
}
}
