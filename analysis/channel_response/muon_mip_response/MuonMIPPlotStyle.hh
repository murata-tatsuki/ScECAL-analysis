#ifndef MUON_MIP_PLOT_STYLE_HH
#define MUON_MIP_PLOT_STYLE_HH

#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "../adc_energy/FastPng.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mmr_plot {

inline constexpr int canvasWidth=2250;
inline constexpr int canvasHeight=1500;

// Landau width is the Landau scale parameter, not a standard deviation.
struct FitInfo {
 std::string status="missing_histogram";
 bool valid=false;
 double mpv=std::numeric_limits<double>::quiet_NaN();
 double mpvError=std::numeric_limits<double>::quiet_NaN();
 double landau=std::numeric_limits<double>::quiet_NaN();
 double gaussian=std::numeric_limits<double>::quiet_NaN();
 double chi2ndf=std::numeric_limits<double>::quiet_NaN();
};

inline std::string fitNumber(double value,int precision=4){
 if(!std::isfinite(value))return "n/a";
 std::ostringstream text;
 text<<std::setprecision(precision)<<value;
 return text.str();
}

// Curves must have the original (unnormalized) histogram units. Sampling the
// TF1 rather than changing its area also supports TF1 objects read from ROOT:
// their saved sampled values can be evaluated without the original callback.
// Histogram clones stay alive in lifetime until after the canvas is cleared;
// graph/legend DrawClone objects belong to the pad.
inline void drawResponsePanel(TCanvas& canvas,int pad,
 const TH1D* first,const TH1D* second,
 const std::string& title,const std::string& unit,
 std::vector<std::unique_ptr<TH1D>>& lifetime,
 const FitInfo* firstFit=nullptr,const FitInfo* secondFit=nullptr,
 const TF1* firstCurve=nullptr,const TF1* secondCurve=nullptr,
 const std::string& firstLabel="data",const std::string& secondLabel="MC"){
 canvas.cd(pad);
 const std::array<const TH1D*,2> sources={first,second};
 const std::array<const FitInfo*,2> fits={firstFit,secondFit};
 const std::array<const TF1*,2> curves={firstCurve,secondCurve};
 const std::array<std::string,2> labels={firstLabel,secondLabel};
 const std::array<int,2> histogramColors={kBlack,kRed+1};
 const std::array<int,2> fitColors={TColor::GetColor("#707070"),TColor::GetColor("#E87979")};
 std::array<TH1D*,2> histograms{};
 std::array<TGraph*,2> drawnCurves{};
 std::array<double,2> totals={1.,1.};
 double maximum=0;
 for(int sample=0;sample<2;++sample){
  if(!sources[sample])continue;
  lifetime.emplace_back(static_cast<TH1D*>(sources[sample]->Clone(
   ("response_plot_"+std::to_string(pad)+"_"+std::to_string(sample)).c_str())));
  auto* histogram=lifetime.back().get();
  histograms[sample]=histogram;
  histogram->SetDirectory(nullptr);
  const double total=histogram->Integral(0,histogram->GetNbinsX()+1);
  if(total>0){totals[sample]=total;histogram->Scale(1./total);}
  histogram->SetStats(false);
  histogram->SetLineColor(histogramColors[sample]);
  histogram->SetLineWidth(2);
  histogram->SetLineStyle(1);
  maximum=std::max(maximum,histogram->GetMaximum());
 }
 auto* frame=histograms[0]?histograms[0]:histograms[1];
 if(!frame){
  TLatex text;
  text.SetNDC();text.SetTextSize(.030);
  text.DrawLatex(.13,.5,"No usable spectrum");
  return;
 }
 frame->SetTitle((title+";"+unit+";Fraction / bin").c_str());
 frame->SetMaximum(std::max(.01,1.35*maximum));
 frame->Draw("HIST");
 if(histograms[0]&&histograms[1])histograms[1]->Draw("HIST SAME");
 for(int sample=0;sample<2;++sample){
  if(!histograms[sample]||!curves[sample]||(fits[sample]&&!fits[sample]->valid))continue;
  const auto* curve=curves[sample];
  const int points=std::max(2,curve->GetNpx()+1);
  TGraph display(points);
  display.SetName(("response_fit_"+std::to_string(pad)+"_"+std::to_string(sample)).c_str());
  for(int point=0;point<points;++point){
   const double x=curve->GetXmin()+(curve->GetXmax()-curve->GetXmin())*point/(points-1.);
   display.SetPoint(point,x,curve->Eval(x)/totals[sample]);
  }
  display.SetLineColor(fitColors[sample]);
  display.SetLineWidth(2);display.SetLineStyle(2);
  drawnCurves[sample]=static_cast<TGraph*>(display.DrawClone("L SAME"));
 }
 int rows=0;
 for(int sample=0;sample<2;++sample){
  if(!histograms[sample])continue;
  ++rows;
  if(fits[sample])rows+=fits[sample]->valid?3:1;
 }
 const double top=.88,bottom=std::max(.42,top-.039*rows);
 TLegend legend(.43,bottom,.89,top);
 legend.SetMargin(.14);
 legend.SetBorderSize(0);legend.SetFillColor(kWhite);
 legend.SetFillStyle(1001);legend.SetTextFont(42);legend.SetTextSize(.027);
 for(int sample=0;sample<2;++sample){
  if(!histograms[sample])continue;
  legend.AddEntry(histograms[sample],
   (labels[sample]+(fits[sample]?" histogram":"")).c_str(),"l");
  if(!fits[sample])continue;
  const auto& fit=*fits[sample];
  if(!fit.valid){
   const std::string status="fit: "+(fit.status.empty()?"unavailable":fit.status);
   auto* entry=legend.AddEntry(static_cast<TObject*>(nullptr),status.c_str(),"");
   // Precision 1 makes status identifiers literal text (underscores included).
   entry->SetTextFont(41);
   continue;
  }
  const std::string result=labels[sample]+" fit: MPV="+fitNumber(fit.mpv)+
   " #pm "+fitNumber(fit.mpvError)+" "+unit;
  legend.AddEntry(drawnCurves[sample],result.c_str(),drawnCurves[sample]?"l":"");
  const std::string widths="w_{L}="+fitNumber(fit.landau)+
   ", #sigma_{G}="+fitNumber(fit.gaussian)+" "+unit;
  legend.AddEntry(static_cast<TObject*>(nullptr),widths.c_str(),"");
  const std::string quality="#chi^{2}/ndf="+fitNumber(fit.chi2ndf,3);
  legend.AddEntry(static_cast<TObject*>(nullptr),quality.c_str(),"");
 }
 legend.DrawClone();
}


inline void drawMpvProfiles(TCanvas& canvas,int cid,bool matched,
 const std::array<std::array<const TGraphAsymmErrors*,2>,4>& graphs,
 const std::string& outputPng){
 canvas.Clear();canvas.SetCanvasSize(canvasWidth,canvasHeight);canvas.Divide(2,2);
 std::vector<std::unique_ptr<TGraphAsymmErrors>> owned;
 // Clear pad references before locally owned graph clones are destroyed.
 struct ClearCanvas {TCanvas& canvas;~ClearCanvas(){canvas.Clear();}} clear{canvas};
 const std::array<std::string,4> axes={
  "Strip longitudinal position [mm]","Strip transverse position [mm]",
  "Temperature [C]","Run"};
 for(int axis=0;axis<4;++axis){
  canvas.cd(axis+1);
  double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  bool any=false;
  for(const auto* graph:graphs[axis]){
   if(!graph||!graph->GetN())continue;
   any=true;
   for(int point=0;point<graph->GetN();++point){
    double x=0,y=0;graph->GetPoint(point,x,y);
    xmin=std::min(xmin,x-graph->GetErrorXlow(point));
    xmax=std::max(xmax,x+graph->GetErrorXhigh(point));
    ymin=std::min(ymin,y-graph->GetErrorYlow(point));
    ymax=std::max(ymax,y+graph->GetErrorYhigh(point));
   }
  }
  if(!any){
   TLatex text;text.SetNDC();text.SetTextSize(.030);
   text.DrawLatex(.13,.55,"No valid binwise MPV fit");
   text.DrawLatex(.13,.48,axes[axis].c_str());
   continue;
  }
  const double dx=std::max(axis==3?1.:.25,.1*(xmax-xmin));
  const double dy=std::max(.025,.15*(ymax-ymin));
  const double topDy=std::max(.025,.45*(ymax-ymin));
  TGraphAsymmErrors* frame=nullptr;
  std::array<TGraphAsymmErrors*,2> displayed{};
  for(int sample=0;sample<2;++sample){
   const auto* source=graphs[axis][sample];
   if(!source||!source->GetN())continue;
   owned.emplace_back(static_cast<TGraphAsymmErrors*>(source->Clone(
    ("mpv_plot_"+std::to_string(axis)+"_"+std::to_string(sample)).c_str())));
   auto* graph=owned.back().get();displayed[sample]=graph;
   graph->SetLineColor(sample?kRed+1:kBlack);
   graph->SetMarkerColor(sample?kRed+1:kBlack);
   graph->SetMarkerStyle(sample?24:20);
   graph->SetTitle(("Cell "+std::to_string(cid)+", "+(matched?"matched":"unmatched")+
    ";"+axes[axis]+";Landau MPV / current MIP").c_str());
   if(!frame){
    graph->SetMinimum(ymin-dy);graph->SetMaximum(ymax+topDy);
    graph->Draw("AP");graph->GetXaxis()->SetLimits(xmin-dx,xmax+dx);
    frame=graph;
   }else graph->Draw("P SAME");
  }
  TLegend legend(.74,.73,.89,.87);
  legend.SetBorderSize(0);
  for(int sample=0;sample<2;++sample)
   if(displayed[sample])legend.AddEntry(displayed[sample],sample?"MC":"data","lp");
  legend.DrawClone();
  TLatex text;text.SetNDC();text.SetTextSize(.030);
  text.DrawLatex(.13,.91,"Binwise Landau x Gaussian fits");
 }
 fastPng(canvas,outputPng.c_str());
}

} // namespace mmr_plot

#endif
