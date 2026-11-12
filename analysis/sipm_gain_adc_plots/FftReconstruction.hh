#pragma once
#include "../sipm_gain/common/GainIO.hh"
#include <TH1.h>
#include <TVirtualFFT.h>
#include <cmath>

namespace adc_plots {
// A display-only band-pass inverse transform, with phase from the original ADC.
// This neither fits a density nor changes the saved FFT gain/status.
struct FftReconstruction {
 std::unique_ptr<TGraph> curve, envelope, filter, windowedResidual;
 gain::Json info={{"status","unavailable"}};
};
inline FftReconstruction reconstructFft(const TH1& roi,const TGraph& envelope,
                                      const TGraph& savedPower,const gain::Json& estimate){
 using namespace gain;
 FftReconstruction out;
 double g=num(estimate,"gain"),dx=roi.GetBinWidth(1);
 int n=roi.GetNbinsX();
 if(!std::isfinite(g)||g<=2*dx||n<3){out.info["reason"]="missing_or_invalid_fft_gain";return out;}
 require(envelope.GetN()==n,"FFT envelope/ADC dimension mismatch");
 const double pi=std::acos(-1.);
 int nfft=16;while(nfft<8*n)nfft*=2;
 Vec z(nfft,0),window(n),x(n),env(n);double weight=0;
 for(int i=0;i<n;++i){
  x[i]=roi.GetBinCenter(i+1);env[i]=envelope.GetPointY(i);
  require(std::isfinite(env[i])&&env[i]>=0&&std::abs(envelope.GetPointX(i)-x[i])<1e-9&&
          std::abs(roi.GetBinWidth(i+1)-dx)<1e-9,"Invalid FFT envelope or ADC spacing");
  window[i]=.5-.5*std::cos(2*pi*i/(n-1));weight+=window[i]*window[i];
  z[i]=(roi.GetBinContent(i+1)-env[i])/std::sqrt(std::max(env[i],1.))*window[i];
 }
 std::unique_ptr<TVirtualFFT> forward(TVirtualFFT::FFT(1,&nfft,"R2C ES K"));
 require(bool(forward),"ROOT FFTW backend is unavailable");
 forward->SetPoints(z.data());forward->Transform();
 Vec re(nfft/2+1),im(re.size()),frequencies(re.size()),mask(re.size());
 for(int k=0;k<=nfft/2;++k){forward->GetPointComplex(k,re[k],im[k]);frequencies[k]=k/(nfft*dx);}
 // Reproduce the saved spectrum before filtering; its power alone lacks phase.
 double maxPowerError=0;
 for(int j=0;j<savedPower.GetN();++j){
  double fk=savedPower.GetPointX(j)*nfft*dx;int k=std::lround(fk);
  require(k>=0&&k<=nfft/2&&std::abs(k-fk)<1e-7,"Saved FFT frequency grid mismatch");
  double power=(re[k]*re[k]+im[k]*im[k])/std::max(weight,1.);
  double error=std::abs(power-savedPower.GetPointY(j))/std::max(1.,std::abs(savedPower.GetPointY(j)));
  maxPowerError=std::max(maxPowerError,error);
 }
 require(savedPower.GetN()>0&&maxPowerError<1e-7,"Recomputed FFT does not reproduce saved power");
 const double f0=1/g,halfWidth=std::min(2/(n*dx),.45*f0),windowFloor=.25;
 int retained=0;
 for(int k=0;k<=nfft/2;++k){
  double distance=std::abs(frequencies[k]-f0)/halfWidth;
  // Flat central band and a cosine edge avoid a sharp-cutoff ringing display.
  mask[k]=distance<=.5?1.:distance<1?.5*(1+std::cos(pi*(2*distance-1))):0.;
  re[k]*=mask[k];im[k]*=mask[k];retained+=mask[k]>0;
 }
 require(retained>0,"No FFT bins in the selected gain band");
 std::unique_ptr<TVirtualFFT> inverse(TVirtualFFT::FFT(1,&nfft,"C2R ES K"));
 require(bool(inverse),"ROOT inverse FFT is unavailable");
 inverse->SetPointsComplex(re.data(),im.data());inverse->Transform();
 Vec filtered(nfft);inverse->GetPoints(filtered.data());for(double& v:filtered)v/=nfft;
 Vec selectedX,counts;
 for(int i=0;i<n;++i)if(window[i]>=windowFloor){
  selectedX.push_back(x[i]);counts.push_back(env[i]+filtered[i]/window[i]*std::sqrt(std::max(env[i],1.)));
  require(std::isfinite(counts.back()),"Nonfinite FFT reconstruction");
 }
 require(!counts.empty(),"No stable inverse-window region");
 out.curve=std::make_unique<TGraph>(selectedX.size(),selectedX.data(),counts.data());
 out.curve->SetName("fft_reconstruction");out.curve->SetTitle("FFT band-pass reconstruction (not a fit);HG ADC;Entries / ADC");
 out.envelope.reset(static_cast<TGraph*>(envelope.Clone("fft_envelope")));
 out.filter=std::make_unique<TGraph>(frequencies.size(),frequencies.data(),mask.data());
 out.filter->SetName("fft_filter");out.filter->SetTitle("Gain-selected frequency filter;frequency [1/ADC];weight");
 out.windowedResidual=std::make_unique<TGraph>(n,x.data(),filtered.data());
 out.windowedResidual->SetName("fft_filtered_windowed_residual");
 out.windowedResidual->SetTitle("Inverse FFT before undoing Hann/whitening;HG ADC;windowed whitened residual");
 out.info={{"status","drawn"},{"gain_ADC",g},{"gain_status",str(estimate,"status")},
  {"algorithm","gain-selected fundamental band-pass inverse FFT; saved envelope + dewhitened, untapered residual"},
  {"fundamental_frequency_per_ADC",f0},{"band_half_width_per_ADC",halfWidth},{"harmonics",1},
  {"filter","unity within half the half-width, cosine roll-off to zero at full half-width"},
  {"window","Hann; same as original FFT"},{"minimum_window_for_display",windowFloor},
  {"fft_size",nfft},{"retained_positive_frequency_bins",retained},{"inverse_normalization",1./nfft},
  {"saved_power_max_relative_error",maxPowerError},{"roi_adc_min",x.front()},{"roi_adc_max",x.back()},
  {"display_adc_min",selectedX.front()},{"display_adc_max",selectedX.back()},
  {"negative_points",std::count_if(counts.begin(),counts.end(),[](double v){return v<0;})},
  {"note","Display only. Original complex phase/amplitude retained in the selected band; no peak-height fitting, clipping, area rescaling, or gain update. ROI edges omitted to avoid dividing by a small Hann window. Higher harmonics are not retained."}};
 return out;
}
}
