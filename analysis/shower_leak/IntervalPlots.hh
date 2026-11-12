#pragma once
#include "Common.hh"
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TLine.h>

// Plot coordinates: two categorical slots at negative u; u=log10(f/0.001)
// for the physical logarithmic axis. No event selection or fit lives here.
struct IntervalValue {
  int bin, sample;
  bool lowStatistics;
  double value, error;
};
using IntervalValues = std::array<std::vector<IntervalValue>, 4>;
inline const std::array<double, 9> intervalEdges =
    {0, .001, .003, .01, .03, .05, .1, .2, 1};
inline std::pair<double, double> intervalDisplayBounds(int bin) {
  if (bin == 0) return {-1.0, -1.0};
  if (bin == 1) return {-.72, 0};
  return {std::log10(intervalEdges.at(bin-1)/.001),
          std::log10(intervalEdges.at(bin)/.001)};
}
inline double intervalDisplayCenter(int bin) {
  auto [lo, hi] = intervalDisplayBounds(bin);
  return (lo + hi)/2;
}

// Recover full-precision y values/errors from saved canvases, including legacy
// categorical canvases. This allows display-only updates without refitting.
inline IntervalValues readIntervalCanvas(TCanvas& canvas) {
  IntervalValues values;
  const bool hybrid = canvas.GetListOfPrimitives()->FindObject("interval_axis") != nullptr;
  for (int q=0; q<4; ++q) {
    auto* pad = canvas.GetPad(q+1);
    require(pad != nullptr, "Missing interval pad");
    TIter next(pad->GetListOfPrimitives());
    while (auto* item = next()) {
      auto* graph = dynamic_cast<TGraph*>(item);
      if (!graph) continue;
      int sample = graph->GetMarkerColor() == kBlue+1 ? 0 :
                   graph->GetMarkerColor() == kRed+1 ? 1 : -1;
      require(sample >= 0, "Unexpected interval marker color");
      for (int j=0; j<graph->GetN(); ++j) {
        const double x=graph->GetPointX(j);
        int bin=std::lround(x);
        if (hybrid && x<0) {
          bin=x<-.5?0:1;
        } else if (hybrid) {
          double distance=1e9;
          for (int b=0;b<9;++b) if (std::abs(x-intervalDisplayCenter(b))<distance) {
            distance=std::abs(x-intervalDisplayCenter(b)); bin=b;
          }
          require(distance<1e-8, "Unknown hybrid interval position");
        }
        require(bin>=0 && bin<9, "Unknown interval bin");
        double error=0;
        if (auto* g=dynamic_cast<TGraphErrors*>(graph)) error=g->GetErrorY(j);
        else if (auto* g=dynamic_cast<TGraphAsymmErrors*>(graph)) error=g->GetErrorYhigh(j);
        else require(false, "Missing interval errors");
        values[q].push_back({bin,sample,graph->GetMarkerStyle()==24,
                             graph->GetPointY(j),error});
      }
    }
  }
  return values;
}

inline void drawIntervalCanvas(TCanvas& canvas, const IntervalValues& values,
                               const std::string& caption) {
  canvas.SetEditable(true);
  canvas.Clear();
  canvas.cd();
  canvas.Divide(2,2);
  auto* metadata=new TNamed("interval_axis","hybrid-log-v1: categorical 0,(0,0.001]; u=log10(f/0.001) above 0.001");
  metadata->SetBit(TObject::kCanDelete);
  canvas.GetListOfPrimitives()->Add(metadata);
  const char* titles[]={"Mean calibrated signal [GeV]", "Gaussian #sigma / mean",
                        "(q84-q16) / (2 median)", "RMS / mean"};
  for (int q=0;q<4;++q) {
    auto* pad=canvas.cd(q+1);
    pad->SetLeftMargin(.19); pad->SetRightMargin(.05);
    pad->SetBottomMargin(.23); pad->SetTopMargin(.24);
    double ymax=0;
    for (const auto& v:values[q]) ymax=std::max(ymax,v.value+v.error);
    ymax=std::max(.01,ymax*1.15);
    auto* frame=new TH1D(Form("interval_frame_%d",q),Form(";;%s",titles[q]),1,-1.30,3.10);
    frame->SetDirectory(nullptr); frame->SetBit(TObject::kCanDelete);
    frame->SetMinimum(0); frame->SetMaximum(ymax);
    frame->GetXaxis()->SetLabelSize(0); frame->GetXaxis()->SetTickLength(0);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetYaxis()->SetLabelSize(.055);frame->GetYaxis()->SetTitleSize(.060);frame->GetYaxis()->CenterTitle(true);
    frame->Draw();
    // Physical logarithmic tick labels start at 0.001; the two slots on the
    // left have deliberately categorical widths (about one third of the total axis width).
    auto* axis=new TGaxis(0,0,3,0,.001,1,510,"G");
    axis->SetLabelFont(42); axis->SetLabelSize(.055);
    axis->SetLabelOffset(.012); axis->SetBit(TObject::kCanDelete); axis->Draw();
    TLatex text; text.SetTextFont(42); text.SetTextSize(.050); text.SetTextAlign(23);
    text.DrawLatex(-1.0,-.035*ymax,"0");
    text.SetTextSize(.040); text.DrawLatex(-.36,-.18*ymax,"(0, 0.001]");
    text.SetNDC(); text.SetTextSize(.046);
    text.DrawLatex(.55,.045,"Indicator signal fraction (log scale above 0.001)");
    auto* divider=new TLine(0,0,0,ymax);
    divider->SetLineStyle(3); divider->SetLineColor(kGray+1);
    divider->SetBit(TObject::kCanDelete); divider->Draw();
    auto* legend=new TLegend(.15,.79,.95,.89);
    legend->SetNColumns(4); legend->SetBorderSize(0);
    legend->SetFillStyle(0); legend->SetTextSize(.030);
    legend->SetBit(TObject::kCanDelete);
    for (int low=0;low<2;++low) for (int sample=0;sample<2;++sample) {
      auto* graph=new TGraphAsymmErrors;
      graph->SetName(Form("interval_%d_%d_%d",q,sample,low));
      graph->SetMarkerColor(sample?kRed+1:kBlue+1);
      graph->SetLineColor(sample?kRed+1:kBlue+1);
      graph->SetMarkerStyle(low?24:20); graph->SetMarkerSize(1.1);
      graph->SetBit(TObject::kCanDelete);
      for (const auto& v:values[q]) if (v.sample==sample && v.lowStatistics==bool(low)) {
        auto [lo,hi]=intervalDisplayBounds(v.bin);
        const double x=intervalDisplayCenter(v.bin);
        int j=graph->GetN(); graph->SetPoint(j,x,v.value);
        graph->SetPointError(j,x-lo,hi-x,v.error,v.error);
      }
      graph->Draw("P SAME");
      legend->AddEntry(graph,Form("%s%s",sample?"MC":"data",low?" (N<100)":""),"p");
    }
    legend->Draw();
    text.SetTextAlign(13); text.SetTextSize(.029);
    text.DrawLatex(.14,.975,caption.c_str());
    if (q==1) {
      text.SetTextSize(.025);
      text.DrawLatex(.14,.925,"Missing Gaussian point: N<100 or invalid fit");
    }
    pad->RedrawAxis();
  }
  canvas.Modified(); canvas.Update();
}
