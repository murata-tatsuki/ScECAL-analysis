#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

inline double missing() { return std::numeric_limits<double>::quiet_NaN(); }
struct Moments {
    long long n=0;
    double mean=0, m2=0;
    void add(double x) { ++n; double d=x-mean; mean+=d/n; m2+=d*(x-mean); }
    double sem() const { return n>1 ? std::sqrt(std::max(0.,m2)/(n*(n-1.))) : missing(); }
};
struct Distribution {
    Moments m;
    std::vector<double> values;
    void add(double x) { m.add(x); values.push_back(x); }
    void sort() { std::sort(values.begin(),values.end()); }
    // Exact empirical quantiles: linear interpolation at (N-1)*p (type 7).
    double quantile(double p) const {
        if(values.empty()) return missing();
        double rank=(values.size()-1)*p; size_t k=static_cast<size_t>(rank);
        return values[k]+(values[std::min(k+1,values.size()-1)]-values[k])*(rank-k);
    }
};
struct Regression {
    long long n=0;
    double mx=0, my=0, xx=0, xy=0, yy=0;
    void add(double x,double y) {
        ++n; double dx=x-mx,dy=y-my; mx+=dx/n; my+=dy/n;
        xx+=dx*(x-mx); xy+=dx*(y-my); yy+=dy*(y-my);
    }
    bool valid() const { return n>=3 && xx>1e-12; }
    double slope() const { return valid()?xy/xx:missing(); }
    double at(double x) const { return valid()?my+slope()*(x-mx):missing(); }
    double errorAt(double x) const {
        return valid()?std::sqrt(std::max(0.,yy-xy*xy/xx)/(n-2)*(1./n+(x-mx)*(x-mx)/xx)):missing();
    }
};
struct ResidualBin {
    Moments x;
    Distribution delta;
    void add(double a,double d) { x.add(a); delta.add(d); }
};
struct SwitchBin {
    Moments x;
    Distribution hg,lg,selected;
    void add(double a,double h,double l,double s) { x.add(a); hg.add(h); lg.add(l); selected.add(s); }
};
