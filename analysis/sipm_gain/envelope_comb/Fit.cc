#include "Numerics.hh"
namespace gain::detail {
Fit envelopeComb(const Vec& x,const Vec& y,const Vec& env,const Config& c,const Vec& seeds,bool fast){
    return comb(x,y,env,c,seeds,true,fast);
}
}
