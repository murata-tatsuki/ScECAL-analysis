#pragma once
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace fs=std::filesystem;
inline void check(bool ok,const std::string& s) {if(!ok)throw std::runtime_error(s);}
struct Page {
  TCanvas canvas;
  TPad header,body;
  TLegend legend;
  Page(const std::string& name,const std::string& title,const std::string& subtitle,int nx,int ny):
    canvas(name.c_str(),"",nx>2?2400:1800,ny>2?1800:1100),
    header((name+"_header").c_str(),"",0,.84,1,1),body((name+"_body").c_str(),"",0,0,1,.84),
    legend(.27,.05,.73,.32) {
    canvas.cd();header.Draw();body.Draw();header.cd();
    const double fit=.90*(nx>2?2400.:1800.)/(.16*(ny>2?1800.:1100.)*.55);
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(std::min(.21,fit/std::max<size_t>(1,title.size())));text.DrawLatex(.03,.72,title.c_str());
    text.SetTextSize(std::min(.15,fit/std::max<size_t>(1,subtitle.size())));text.DrawLatex(.03,.43,subtitle.c_str());
    legend.SetNColumns(2);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.17);
    body.Divide(nx,ny,.008,.008);
  }
  void save(TFile& out,const fs::path& dir) {
    header.cd();legend.Draw();canvas.cd();canvas.SaveAs((dir/(std::string(canvas.GetName())+".png")).c_str());
    out.cd();canvas.Write();
  }
};
inline void pad() {gPad->SetLeftMargin(.15);gPad->SetRightMargin(.08);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.12);gPad->SetTicks();}
inline std::unique_ptr<TH1> clone(TFile& f,const std::string& key,bool normalize) {
  auto* h=dynamic_cast<TH1*>(f.Get(key.c_str()));check(h,"Missing histogram "+key);
  static int serial=0;
  auto v=std::unique_ptr<TH1>(dynamic_cast<TH1*>(h->Clone(("plot_"+std::to_string(serial++)).c_str())));
  v->SetDirectory(nullptr);const double sum=v->Integral(0,v->GetNbinsX()+1);
  if(normalize&&sum>0)v->Scale(1./sum);
  return v;
}
inline void pairDraw(TH1& d,TH1& m,const std::string& title,bool log) {
  pad();gPad->SetLogy(log);
  d.SetLineColor(kBlack);m.SetLineColor(kRed+1);d.SetLineWidth(2);m.SetLineWidth(2);m.SetLineStyle(2);
  d.SetTitle(title.c_str());d.GetXaxis()->SetTitleOffset(1.05);d.GetYaxis()->SetTitleOffset(1.5);const double maximum=std::max(d.GetMaximum(),m.GetMaximum());
  d.SetMaximum(maximum*(log?2:1.2)+1e-12);d.SetMinimum(log?std::max(1e-8,maximum*1e-5):0);
  d.Draw("hist");m.Draw("hist same");
}
