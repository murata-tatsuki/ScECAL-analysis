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
double sum(const Vec& v){return std::accumulate(v.begin(),v.end(),0.);}
double square(double x){return x*x;}
Vec slice(const Vec& v,int a,int b){return Vec(v.begin()+a,v.begin()+b);}
int maximum(const Vec& v){return int(std::max_element(v.begin(),v.end())-v.begin());}
double deviation(const Vec& v){if(v.size()<2)return NaN;double m=sum(v)/v.size(),s=0;for(double a:v)s+=square(a-m);return std::sqrt(s/(v.size()-1));}
double quantile(Vec v,double p){if(v.empty())return NaN;std::sort(v.begin(),v.end());double x=p*(v.size()-1);size_t i=size_t(x);return v[i]+(x-i)*(v[std::min(i+1,v.size()-1)]-v[i]);}
double number(const Json& j,const std::string& k,double fallback){return j.contains(k)&&j[k].is_number()?j[k].get<double>():fallback;}
std::string join(const std::vector<std::string>& v){std::string s;for(auto& x:v){if(!s.empty())s+=';';s+=x;}return s;}
Vec smooth(const Vec& y,double sigma){
    int radius=std::max(2,int(std::ceil(4*sigma))),n=y.size();Vec kernel(2*radius+1),out(n);
    for(int j=-radius;j<=radius;++j)kernel[j+radius]=std::exp(-.5*square(j/sigma));
    double z=sum(kernel);for(double& a:kernel)a/=z;
    for(int i=0;i<n;++i)for(int j=-radius;j<=radius;++j)out[i]+=kernel[j+radius]*y[std::clamp(i-j,0,n-1)];
    return out;
}
double cdf(double x){
    // A&S 7.1.26, retained to preserve the original estimator's bin integrals.
    double z=std::abs(x)/std::sqrt(2.),t=1/(1+.3275911*z);
    double erf=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*std::exp(-z*z);
    return .5*(1+(x>0?1:x<0?-1:0)*erf);
}
double deviance(const Vec& y,const Vec& mu){
    double d=0;for(size_t i=0;i<y.size();++i){double m=std::max(mu[i],1e-12);d+=m-y[i];if(y[i]>0)d+=y[i]*std::log(y[i]/m);}return 2*d;
}
Minimum minimize(const std::function<double(const Vec&)>& fun,Vec initial,const Bounds& bounds,bool hesse,int maxCalls){
    std::unique_ptr<ROOT::Math::Minimizer> engine(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
    if(!engine)throw std::runtime_error("ROOT Minuit2 is unavailable");
    engine->SetPrintLevel(-1);engine->SetErrorDef(.5);engine->SetTolerance(.02);
    engine->SetMaxFunctionCalls(maxCalls);engine->SetStrategy(1);
    auto callback=[&](const double* pp){Vec p(pp,pp+initial.size());double v=fun(p);return std::isfinite(v)?v:1e100;};
    ROOT::Math::Functor function(callback,initial.size());engine->SetFunction(function);
    for(size_t i=0;i<initial.size();++i){auto [lo,hi]=bounds[i];double v=std::clamp(initial[i],lo+1e-6*(hi-lo),hi-1e-6*(hi-lo));
        double step=std::max(1e-4,std::min((hi-lo)*.002,std::max(std::abs(v)*.01,.002)));
        engine->SetLimitedVariable(i,"p"+std::to_string(i),v,step,lo,hi);
    }
    bool ok=engine->Minimize();if(hesse&&ok)engine->Hesse();
    Minimum r;r.p.assign(engine->X(),engine->X()+initial.size());r.error.assign(engine->Errors(),engine->Errors()+initial.size());
    r.status=engine->Status();r.covariance=engine->CovMatrixStatus();r.converged=ok&&r.status==0;return r;
}
double interpolate(double q,const Vec& x,const Vec& y){
    if(q<=x.front())return y.front();
    if(q>=x.back())return y.back();
    size_t k=std::upper_bound(x.begin(),x.end(),q)-x.begin();return y[k-1]+(y[k]-y[k-1])*(q-x[k-1])/(x[k]-x[k-1]);
}
std::pair<Vec,Json> envelope(const Vec& x,const Vec& y,const std::string& kind,const Config& c,bool fast){
    double dx=x[1]-x[0],total=std::max(sum(y),1.),center=0,var=0;for(size_t i=0;i<x.size();++i)center+=x[i]*y[i]/total;
    for(size_t i=0;i<x.size();++i)var+=y[i]*square(x[i]-center)/total;
    double scale=std::max(c.gain_max*.8,std::sqrt(var)),width=x.back()-x.front();
    Vec baseline=smooth(y,std::max(2.,c.gain_max/dx));for(double& v:baseline)v=std::max(v,.01);
    if(kind=="led")return {baseline,{{"model","kernel"},{"status",0},{"parameters",Json::array()}}};
    auto predict=[&](const Vec& p){Vec shape(x.size());for(size_t i=0;i<x.size();++i){
        if(kind=="muon"){
            double landau=0;for(int j=0;j<32;++j){double z=-5.+(j+.5)*10./32;
                landau+=TMath::Landau(x[i]-p[2]*z,p[0]+.22278298*p[1],p[1],true)*std::exp(-.5*z*z)*10./32/std::sqrt(2*TMath::Pi());}
            shape[i]=landau*cdf((x[i]-p[3])/p[4])+p[5]/width;
        }else shape[i]=(std::exp(-p[0]*(x[i]-x.front()))+p[1])*cdf((x[i]-p[2])/p[3]);
        shape[i]=std::max(shape[i],1e-15);
    }double norm=sum(shape);for(double& v:shape)v*=total/norm;return shape;};
    Vec initial;Bounds bounds;
    if(kind=="muon"){
        initial={x[maximum(baseline)],scale/3,scale/4,x.front(),std::max(2.,c.gain_max/3),.02};
        bounds={{x.front()-.5*width,x.back()},{2.,width},{1.,width*.6},{x.front()-2*c.gain_max,x.front()+width*.4},{1.,2*c.gain_max},{0.,.5}};
    }else{
        initial={1/std::max(scale,20.),.02,x.front(),std::max(2.,c.gain_max/3)};
        bounds={{0.,.1},{0.,1.},{x.front()-2*c.gain_max,x.front()+width*.5},{1.,2*c.gain_max}};
    }
    auto fit=minimize([&](const Vec& p){return deviance(y,predict(p))/2;},initial,bounds,!fast,3000);
    auto mu=predict(fit.p);return {mu,{{"model",kind=="muon"?"langaus_erf":"exp_erf"},{"status",fit.status},{"converged",fit.converged},
        {"deviance",deviance(y,mu)},{"ndf",int(y.size())-int(initial.size())-1},{"parameters",fit.p}}};
}
}
