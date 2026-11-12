#include "Numerics.hh"
namespace gain::detail {
Fit multiGaussian(const Vec& x,const Vec& y,const Vec& env,const Config& c,const Vec& seeds,bool fast){
    return comb(x,y,env,c,seeds,false,fast);
}
}
