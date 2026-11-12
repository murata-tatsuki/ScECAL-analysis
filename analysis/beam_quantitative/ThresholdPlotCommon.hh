#pragma once
#include "PlotCommon.hh"
#include <TGraphAsymmErrors.h>

inline std::unique_ptr<TGraphAsymmErrors> curve(TFile& f,const std::string& prefix,int color,int marker,const std::string& center="_mean") {
  auto a=clone(f,prefix+center,false),lo=clone(f,prefix+"_low",false),hi=clone(f,prefix+"_high",false);
  auto graph=std::make_unique<TGraphAsymmErrors>();
  for(int i=1;i<=a->GetNbinsX();++i) {
    double v=a->GetBinContent(i),l=lo->GetBinContent(i),h=hi->GetBinContent(i);
    if(!std::isfinite(v)||!std::isfinite(l)||!std::isfinite(h))continue;
    int n=graph->GetN();graph->SetPoint(n,a->GetBinCenter(i),v);
    // Percentile intervals need not contain the original estimate.
    graph->SetPointError(n,0,0,std::max(0.,v-l),std::max(0.,h-v));
  }
  graph->SetLineColor(color);graph->SetMarkerColor(color);graph->SetMarkerStyle(marker);graph->SetLineWidth(2);
  return graph;
}
inline std::string label(const std::string& key) {
  if(key=="uncut")return "No additional cut";
  std::string v=key.substr(4);std::replace(v.begin(),v.end(),'p','.');return "Signal >= "+v+" MIP";
}
inline void frame(TH1& h,const std::string& title,double low,double high) {
  pad();const auto separator=title.find(';');
  const auto heading=title.substr(0,separator);
  h.SetTitle((separator==std::string::npos?std::string(""):title.substr(separator)).c_str());
  h.SetMinimum(low);h.SetMaximum(high);h.GetYaxis()->SetTitleOffset(1.5);h.Draw();
  if(!heading.empty()){TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.05);text.DrawLatex(.15,.94,heading.c_str());}
}

inline double upperBound(const TGraphAsymmErrors& graph) {
  double high=0.;for(int i=0;i<graph.GetN();++i)high=std::max(high,graph.GetPointY(i)+graph.GetErrorYhigh(i));return high;
}

inline std::unique_ptr<TH1> normalizedSignal(TFile& file,const std::string& key,int layer=-1) {
  auto signal=clone(file,key,false),counts=clone(file,"pre/uncut/count_mean",false);
  const double total=layer<0?counts->Integral():counts->GetBinContent(layer+1);
  // Include strips outside the displayed/stored signal range in normalization.
  if(total>0)signal->Scale(1./total);
  return signal;
}
