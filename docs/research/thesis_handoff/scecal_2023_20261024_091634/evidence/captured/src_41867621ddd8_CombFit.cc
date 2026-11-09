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
Fit comb(const Vec& x,const Vec& y,const Vec& env,const Config& c,const Vec& seeds,bool structured,bool fast){
    double dx=x[1]-x[0],total=std::max(sum(y),1.);std::vector<Fit> fits;
    Vec expanded(seeds.begin(),seeds.begin()+std::min(seeds.size(),size_t(fast?1:3)));
    if(!fast){auto original=expanded;for(double g:original)for(double v:{g/2,g*2}){bool distinct=true;for(double a:expanded)distinct &= std::abs(v/a-1)>.03;if(c.gain_min<v&&v<c.gain_max&&distinct)expanded.push_back(v);}}
    for(double g:expanded){
        if(!(c.gain_min<=g&&g<=c.gain_max))continue;
        Vec residual(y.size());for(size_t i=0;i<y.size();++i)residual[i]=(y[i]-env[i])/std::sqrt(std::max(env[i],1.));
        double phase=x.front()+std::fmod(x[maximum(residual)]-x.front(),g);
        int peaks=int(std::ceil((x.back()-x.front())/g))+3,n=x.size(),cols=peaks+1;
        Vec matrix(size_t(n)*cols),sums(cols),heights(cols),mu(n),ratios(n),fresh(cols);
        auto predict=[&](const Vec& p){
            std::fill(sums.begin(),sums.end(),0.);
            for(int j=0;j<peaks;++j){double center=p[1]+(j-1)*p[0],sigma=p[0]*p[2];
                heights[j]=interpolate(center,x,env)*p[0]/dx;
                for(int i=0;i<n;++i){double a=cdf((x[i]+dx/2-center)/sigma)-cdf((x[i]-dx/2-center)/sigma);matrix[size_t(i)*cols+j]=a;sums[j]+=a;}
            }
            if(structured){for(int i=0;i<n;++i){double value=0;for(int j=0;j<peaks;++j)value+=matrix[size_t(i)*cols+j]*heights[j];mu[i]=std::max(p[4]*(p[3]*value+(1-p[3])*env[i]),1e-10);}return mu;}
            heights.back()=.01*total;sums.back()=1.;for(int i=0;i<n;++i)matrix[size_t(i)*cols+peaks]=1./n;
            for(int j=0;j<cols;++j)if(sums[j]<1e-10)heights[j]=0;
            // Fixed iteration count keeps the profiled objective smooth for Minuit.
            for(int iter=0;iter<80;++iter){
                for(int i=0;i<n;++i){double v=0;for(int j=0;j<cols;++j)v+=matrix[size_t(i)*cols+j]*heights[j];ratios[i]=y[i]/std::max(v,1e-10);}
                std::fill(fresh.begin(),fresh.end(),0.);
                for(int i=0;i<n;++i)for(int j=0;j<cols;++j)fresh[j]+=matrix[size_t(i)*cols+j]*ratios[i];
                for(int j=0;j<cols;++j)heights[j]*=fresh[j]/std::max(sums[j],1e-10);
            }
            for(int i=0;i<n;++i){double v=0;for(int j=0;j<cols;++j)v+=matrix[size_t(i)*cols+j]*heights[j];mu[i]=std::max(v,1e-10);}return mu;
        };
        Vec initial={g,phase,.16};Bounds bounds={{std::max(c.gain_min,g*.88),std::min(c.gain_max,g*1.12)},{phase-g*.55,phase+g*.55},{.065,.42}};
        if(structured){initial.insert(initial.end(),{.8,1.});bounds.insert(bounds.end(),{{0.,1.},{.3,3.}});}
        auto fit=minimize([&](const Vec& p){return deviance(y,predict(p))/2;},initial,bounds,!fast,fast?2500:6000);
        Vec predicted=predict(fit.p);double d=deviance(y,predicted);int npars=initial.size(),consecutive=0,current=0;Vec areas;
        if(!structured){areas.assign(heights.begin(),heights.begin()+peaks);double threshold=10.;
            for(int j=0;j<peaks;++j){double center=fit.p[1]+(j-1)*fit.p[0];if(center>=x.front()&&center<=x.back()){++npars;threshold=std::max(threshold,.01*heights[j]*sums[j]);}}
            ++npars;for(int j=0;j<peaks;++j){double center=fit.p[1]+(j-1)*fit.p[0];current=(center>=x.front()&&center<=x.back()&&heights[j]*sums[j]>threshold)?current+1:0;consecutive=std::max(consecutive,current);}
        }
        fits.push_back({{{"gain",fit.p[0]},{"fit_error",fit.error[0]},{"fit_status",fit.status},{"covariance_status",fit.covariance},
            {"converged",fit.converged},{"deviance",d},{"ndf",n-npars},{"bic",d+npars*std::log(n)},
            {"parameters",fit.p},{"peak_areas",areas},{"consecutive_peaks",consecutive},
            {"boundary_gain",std::min(fit.p[0]-bounds[0].first,bounds[0].second-fit.p[0])<.002*g},{"width_fraction",fit.p[2]}},predicted});
    }
    if(fits.empty())return {{{"gain",NaN},{"fit_status",-1},{"converged",false}},Vec(y.size(),0)};
    int best=-1;bool hasConverged=false;for(auto& f:fits)hasConverged|=f.values["converged"].get<bool>();
    for(size_t i=0;i<fits.size();++i)if(!hasConverged||fits[i].values["converged"].get<bool>())if(best<0||number(fits[i].values,"bic")<number(fits[best].values,"bic"))best=i;
    auto result=fits[best];double alternative=std::numeric_limits<double>::infinity();Json trials=Json::array();
    for(auto& f:fits){auto& v=f.values;if(std::abs(number(v,"gain")/number(result.values,"gain")-1)>.05&&v["converged"].get<bool>())alternative=std::min(alternative,number(v,"bic")-number(result.values,"bic"));
        Json t;for(auto k:{"gain","fit_status","converged","deviance","bic","consecutive_peaks","width_fraction"})t[k]=v[k];trials.push_back(t);}
    result.values["alternative_delta_bic"]=alternative;result.values["trials"]=trials;return result;
}
}
