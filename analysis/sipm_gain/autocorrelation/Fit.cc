#include "Numerics.hh"
#include <algorithm>
#include <cmath>
namespace gain::detail {
void autocorrelation(Spectrum& r,const Vec& z,const Vec& window,double dx,const Config& c){
    int n=window.size();
    // Direct evaluation of the linear correlation, restricted to needed lags.
    // This equals zero-padded FFT autocorrelation, without circular wraparound.
    int last=std::min(n-1,int(std::ceil(c.gain_max/dx))+1);Vec ac(last+1);
    for(int lag=0;lag<=last;++lag){double a=0,w=0;for(int i=0;i+lag<n;++i){a+=z[i]*z[i+lag];w+=window[i]*window[i+lag];}ac[lag]=a/std::max(w,1e-12);}
    double norm=std::max(ac[0],1e-12);for(double& a:ac)a/=norm;
    int best=-1;for(int lag=0;lag<=last;++lag)if(lag*dx>=c.gain_min&&lag*dx<=c.gain_max&&lag<n/2.){
        r.lags.push_back(lag*dx);r.acf.push_back(ac[lag]);if(best<0||ac[lag]>ac[best])best=lag;}
    if(best<0)throw std::runtime_error("Empty autocorrelation search band");
    r.acfGain=best*dx;r.acfScore=ac[best];if(best>0&&best<last){double den=ac[best-1]-2*ac[best]+ac[best+1];if(den<0)r.acfGain+=std::clamp(.5*(ac[best-1]-ac[best+1])/den,-.5,.5)*dx;}
}
}
