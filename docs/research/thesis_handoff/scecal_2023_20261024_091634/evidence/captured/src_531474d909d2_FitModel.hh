#pragma once
#include "Common.hh"
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
namespace mip {
// Switch integration coordinates for very narrow Landau widths. A capped
// Gaussian-coordinate grid undersamples the Landau peak when sigma/width is large.
class Convolution {
 double width,center,sigma;bool landauGrid;
 std::vector<double>offset,weight;
public:
 Convolution(double w,double mpv,double sig,double xmax,int precision=80):width(w),center(mpv+.22278298*w),sigma(sig),landauGrid(sig>8*w){
  if(!(w>0&&sig>0))return;
  if(!landauGrid){
   int steps=std::max(int(std::ceil(1.6*precision)),int(std::ceil(40*sig/w)));double dt=16./steps;
   for(int k=0;k<steps;++k){double t=-8+(k+.5)*dt;offset.push_back(sig*t);weight.push_back(std::exp(-.5*t*t)*dt*.3989422804014327);}
  }else{
   const int n=precision>80?32:16;
   static const auto quadratures=[](){std::array<std::pair<std::vector<double>,std::vector<double>>,2>q;for(int i=0;i<2;++i){int n=i?32:16;q[i].first.resize(n);q[i].second.resize(n);TF1::CalcGaussLegendreSamplingPoints(n,q[i].first.data(),q[i].second.data(),1e-14);}return q;}();
   const auto&nodes=quadratures[n==32];
   const double upper=(xmax-center+8*sig)/w,maximumStep=4*sig/w;
   for(double low=-12;low<upper;){
    double high=std::min(upper,low<12?low+3:low+std::min(low,maximumStep));
    double half=(high-low)/2,mid=(high+low)/2;
    for(int j=0;j<n;++j){double z=mid+half*nodes.first[j];offset.push_back(center+w*z);weight.push_back(TMath::Landau(z,0,1,true)*half*nodes.second[j]*.3989422804014327/sig);}
    low=high;
   }
  }
 }
 double operator()(double x)const{
  if(!(width>0&&sigma>=0))return 0;
  if(sigma==0)return TMath::Landau(x,center,width,true);
  double sum=0;
  for(size_t k=0;k<offset.size();++k){if(landauGrid){double z=(x-offset[k])/sigma;if(std::abs(z)<8)sum+=weight[k]*std::exp(-.5*z*z);}else sum+=TMath::Landau(x+offset[k],center,width,true)*weight[k];}
  return sum;
 }
};
// Landau MPV corrected for ROOT's -0.22278298 mode displacement.
inline double density(double x,double w,double mpv,double sigma,int precision=80){return Convolution(w,mpv,sigma,x,precision)(x);}
inline double langau(double*x,double*p){return p[3]*density(x[0],p[0],p[1],p[2]);}
struct InputCondition {Condition c;std::unique_ptr<TH1D>h;};
struct FitSpec {
 int variant=0,rebin=1,precision=80;
 double plateauSigmas=-1,upperScale=1,seedScale=1,thresholdShift=0,sigmaScale=1;
 double rawLow=0,mpvLower=.1;
 bool noGaussian=false,hesse=true;
 std::array<double,3> initial={nan,nan,nan};
};
struct BinnedCondition{Condition c;std::vector<double>counts,eff;double total=0;};
class Likelihood {
public:
 std::vector<double>x;std::vector<BinnedCondition>conditions;double width=0,low=0,high=0;int precision=80;
 std::vector<double> shape(const double*p)const {
  std::vector<double>a(x.size()*3);
  Convolution convolution(p[0],p[1],p[2],x.empty()?high:x.back()+width/2,precision);
  for(size_t i=0;i<x.size();++i)for(int j=0;j<3;++j)a[3*i+j]=convolution(x[i]+(j-1)*width*.3872983346207417);
  return a;
 }
 std::vector<double> expected(const BinnedCondition&c,const std::vector<double>&shape)const {
  std::vector<double>v(x.size());double sum=0;
  for(size_t i=0;i<x.size();++i){v[i]=std::max(1e-300,(shape[3*i]*c.eff[3*i]*5+shape[3*i+1]*c.eff[3*i+1]*8+shape[3*i+2]*c.eff[3*i+2]*5)/18);sum+=v[i];}
  for(double&n:v)n=std::max(1e-300,n*c.total/sum);
  return v;
 }
 double operator()(const double*p)const{auto s=shape(p);double dev=0;for(auto&c:conditions){auto v=expected(c,s);for(size_t i=0;i<v.size();++i){double n=c.counts[i];dev+=2*(v[i]-n+(n>0?n*std::log(n/v[i]):0));}}return std::isfinite(dev)?dev:1e100;}
};
inline Likelihood prepare(const std::vector<InputCondition>&input,const Meta&m,const FitSpec&spec,double maxADC){
 Likelihood q;q.precision=spec.precision;int base=mmr::group(m.CellID)?2:4;int rebin=base*spec.rebin;q.width=2*rebin;
 double limit=maxADC,plateau=0,low=0;
 for(auto&i:input)if(i.h->Integral()>0){low=std::max(low,spec.rawLow*i.c.factor);limit=std::min({limit,(2600-m.pedestal)*i.c.factor,i.h->GetXaxis()->GetXmax()});plateau=std::max(plateau,(m.threshold+spec.plateauSigmas*m.threshold_sigma)*i.c.factor);}
 q.low=spec.plateauSigmas<0?low:std::max(low,plateau);q.high=limit*spec.upperScale;
 if(!(q.high>q.low))return q;
 // Align to existing rebinned edges (-100 origin), keeping whole bins only.
 q.low=-100+std::ceil((q.low+100)/q.width)*q.width;q.high=-100+std::floor((q.high+100)/q.width)*q.width;
 for(double x=q.low+.5*q.width;x<q.high;x+=q.width)q.x.push_back(x);
 for(auto&i:input){BinnedCondition b;b.c=i.c;
  for(double x:q.x){double n=0;int first=i.h->FindFixBin(x-.5*q.width+1e-7);for(int k=0;k<rebin;++k)n+=i.h->GetBinContent(first+k);b.counts.push_back(n);b.total+=n;for(int j=0;j<3;++j){double coordinate=x+(j-1)*q.width*.3872983346207417;double eff=spec.plateauSigmas<0?.5*std::erfc(-((coordinate/i.c.factor-m.threshold-spec.thresholdShift*m.threshold_sigma)/(m.threshold_sigma*spec.sigmaScale))/std::sqrt(2.)):1.;b.eff.push_back(eff);}}
  if(b.total>0)q.conditions.push_back(std::move(b));
 }
 return q;
}
struct Fitted {Result r;Likelihood model;std::vector<double>expected;};
inline Fitted fit(const std::vector<InputCondition>&input,const Meta&m,const mmr::Config&cfg,const FitSpec&spec,double fixedMPV=nan){
 Fitted out;auto&r=out.r;r.CellID=m.CellID;r.variant=spec.variant;r.reason="insufficient_statistics";r.model=spec.noGaussian?1:0;r.raw_fit_low=spec.rawLow;
 double maxADC=cfg.number(mmr::group(m.CellID)?"fit-max-10":"fit-max-15",mmr::group(m.CellID)?1000:2000);
 out.model=prepare(input,m,spec,maxADC);auto&q=out.model;r.fit_low=q.low;r.fit_high=q.high;r.bin_width=q.width;r.conditions=q.conditions.size();
 double total=0;for(auto&c:q.conditions)total+=c.total;r.chnEntries=int(std::min(total,double(INT_MAX)));
 for(auto&i:input)if(i.h->Integral()>0)r.plateau=std::isfinite(r.plateau)?std::max(r.plateau,(m.threshold+4*m.threshold_sigma)*i.c.factor):(m.threshold+4*m.threshold_sigma)*i.c.factor;
 if(total<cfg.integer("min-fit-hits",200)||q.x.size()<8||q.conditions.empty())return out;
 std::vector<double>sum(q.x.size());for(auto&c:q.conditions)for(size_t i=0;i<sum.size();++i)sum[i]+=c.counts[i];
 // ADC data alone determine the seed. No old/MC MPV enters the fit.
 double modal=q.x[std::max_element(sum.begin(),sum.end())-sum.begin()];double seed=std::clamp(modal*spec.seedScale,2.,q.high*.8);
 double lower[3]={.3,spec.mpvLower,.3},upper[3]={q.high*.4,q.high*.95,q.high*.5};double start[3]={std::max(2.,seed*.12),seed,std::max(2.,seed*.2)};
 for(int i=0;i<3;++i)if(std::isfinite(spec.initial[i]))start[i]=spec.initial[i];
 auto minimizer=std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));require(bool(minimizer),"Cannot create Minuit2");minimizer->SetMaxFunctionCalls(3500);minimizer->SetMaxIterations(1500);minimizer->SetTolerance(.02);minimizer->SetPrintLevel(-1);minimizer->SetErrorDef(1);
 ROOT::Math::Functor fn([&](const double*p){return q(p);},3);minimizer->SetFunction(fn);
 for(int i=0;i<3;++i){
  const char* name=i==0?"LandauWidth":i==1?"LandauMPV":"GauSigma";
  // Define fixed variables before any limits: ROOT otherwise retains the old
  // bound and silently clips a requested zero sigma to the previous 0.3 floor.
  if(i==1&&std::isfinite(fixedMPV))minimizer->SetFixedVariable(i,name,fixedMPV);
  else if(i==2&&spec.noGaussian)minimizer->SetFixedVariable(i,name,0);
  else minimizer->SetLimitedVariable(i,name,std::clamp(start[i],lower[i]+.01,upper[i]-.01),std::max(.2,std::abs(start[i])*.05),lower[i],upper[i]);
 }
 bool success=minimizer->Minimize();if(success&&spec.hesse)minimizer->Hesse();r.status=minimizer->Status();r.covariance=minimizer->CovMatrixStatus();r.edm=minimizer->Edm();r.objective=minimizer->MinValue();auto*p=minimizer->X();auto*e=minimizer->Errors();
 r.LandauWidth=p[0];r.LandauMPV=p[1];r.GauSigma=p[2];r.LandauWidth_error=e[0];r.LandauMPV_error=e[1];r.GauSigma_error=e[2];r.relative_error=e[1]/p[1];r.correlation=0;
 for(int i=0;i<3;++i){
  if(!(i==1&&std::isfinite(fixedMPV))&&!(i==2&&spec.noGaussian)){
   double tolerance=std::max(.02,(upper[i]-lower[i])*.001);
   if(p[i]-lower[i]<tolerance)r.boundary_mask|=1<<i;
   if(upper[i]-p[i]<tolerance)r.boundary_mask|=1<<(i+3);
  }
  for(int j=0;j<3;++j)r.cov[3*i+j]=minimizer->CovMatrix(i,j);
  if(i!=1&&e[i]>0&&e[1]>0)r.correlation=std::max(r.correlation,std::abs(r.cov[3*i+1]/(e[i]*e[1])));
 }
 r.boundary=r.boundary_mask!=0;
 auto shape=q.shape(p);out.expected.assign(q.x.size(),0);double chi=0;int groups=0;
 for(auto&c:q.conditions){auto expected=q.expected(c,shape);double n=0,mu=0,lastN=0,lastMu=0;bool have=false;
  for(size_t i=0;i<expected.size();++i)if(expected[i]<.001)r.model_outlier_hits+=int(c.counts[i]);
  // Pearson chi2 with adjacent bins pooled to expected >=5, separately per stratum.
  for(size_t i=0;i<expected.size();++i){out.expected[i]+=expected[i];n+=c.counts[i];mu+=expected[i];if(mu>=5){if(have){chi+=(lastN-lastMu)*(lastN-lastMu)/lastMu;++groups;}lastN=n;lastMu=mu;n=mu=0;have=true;}}
  if(have){lastN+=n;lastMu+=mu;chi+=(lastN-lastMu)*(lastN-lastMu)/lastMu;++groups;}
  else if(mu>0){chi+=(n-mu)*(n-mu)/mu;++groups;}
 }
 r.ChiSquare=chi;r.NDF=groups-int(q.conditions.size())-(3-int(std::isfinite(fixedMPV))-int(spec.noGaussian));
 TF1 convolution("convolution",langau,0,q.high,4);convolution.SetParameters(p[0],p[1],p[2],1);r.convolution_peak=convolution.GetMaximumX(0,q.high);
 r.reason="candidate";r.decision=0;
 auto fail=[&](const std::string&reason,int decision=0){r.reason=reason;r.decision=decision;};
 if(!success||r.status!=0)fail("fit_failed",2);
 else if(r.covariance<3)fail("covariance_unreliable");
 else if(r.NDF<=0)fail("insufficient_chi2_groups");
 else if(r.boundary)fail("parameter_boundary");
 else if(!std::isfinite(r.relative_error)||r.relative_error<=0||r.relative_error>cfg.number("max-relative-error",.10))fail("MPV_error_too_large");
 else if(r.correlation>cfg.number("max-correlation",.98))fail("MPV_width_degeneracy");
 else if(r.model_outlier_hits>=std::max(3.,r.chnEntries*.005))fail("threshold_model_incompatible",2);
 else if(r.ChiSquare/r.NDF>(mmr::group(m.CellID)?2.:1.7))fail("chi2_rejected",2);
 return out;
}
// Independent data-only starts; choose the lowest converged deviance. The
// previous calibration and MC parameters never enter a seed or constraint.
inline Fitted fitBest(const std::vector<InputCondition>&input,const Meta&m,const mmr::Config&cfg,FitSpec spec,double fixedMPV=nan){
 auto best=fit(input,m,cfg,spec,fixedMPV);
 if(!std::isfinite(best.r.objective))return best;
 auto q=prepare(input,m,spec,cfg.number(mmr::group(m.CellID)?"fit-max-10":"fit-max-15",mmr::group(m.CellID)?1000:2000));
 std::vector<double>counts(q.x.size());for(auto&c:q.conditions)for(size_t i=0;i<counts.size();++i)counts[i]+=c.counts[i];
 double peak=q.x[std::max_element(counts.begin(),counts.end())-counts.begin()];
 for(auto seed:std::vector<std::array<double,3>>{{std::max(2.,peak*.08),peak*.8,std::max(2.,peak*.55)},{std::max(2.,peak*.15),peak*1.2,std::max(2.,peak*.05)}}){
  spec.initial=seed;auto alt=fit(input,m,cfg,spec,fixedMPV);
  if(std::isfinite(alt.r.objective)&&alt.r.status==0&&(best.r.status!=0||alt.r.objective<best.r.objective))best=std::move(alt);
 }
 return best;
}
inline void saveFit(TDirectory&d,const Fitted&f,const std::string&name){
 const auto&r=f.r;if(!std::isfinite(r.LandauMPV))return;d.cd();
 TF1 fun((name+"_langau").c_str(),langau,0,r.fit_high,4);fun.SetParameters(r.LandauWidth,r.LandauMPV,r.GauSigma,1);fun.SetParNames("LandauWidth","LandauMPV","GauSigma","density_scale");fun.SetNpx(1600);fun.Save(0,r.fit_high,0,0,0,0);fun.Write();
 TGraph graph(f.model.x.size(),f.model.x.data(),f.expected.data());graph.SetName((name+"_expected").c_str());graph.Write();
 auto h=hist(name+"_used",int(f.model.x.size()),r.fit_low,r.fit_high);for(size_t i=0;i<f.model.x.size();++i){double n=0;for(auto&c:f.model.conditions)n+=c.counts[i];h->SetBinContent(i+1,n);h->SetBinError(i+1,std::sqrt(n));}h->SetEntries(r.chnEntries);h->Write();
}
}
