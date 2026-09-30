#pragma once
#include "StudyCommon.hh"
#include "../../calibration_parameters/EBUdecode.h"
#include <TBox.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <array>
#include <set>
#include <unistd.h>

// Use the same strip mapping and pitch bins as calibration_parameters maps.
// The electronics-index histograms remain the authoritative rate values.
inline void drawTailLayerMaps(const TH2D& data, const TH2D& mc,
                              const std::set<int>& common,
                              const fs::path& out, double beam) {
 TDirectory::TContext context;
 fs::create_directories(out/"figures");
 const auto temporary=out/("hg_lg_tail_layer_maps.root.tmp."+std::to_string(getpid()));
 TFile file(temporary.c_str(),"RECREATE");
 require(!file.IsZombie(),"Cannot create layer map ROOT");
 std::array<std::array<std::unique_ptr<TH2D>,30>,2> maps;
 std::array<std::unique_ptr<TH2D>,30> masks;
 const std::array<const TH2D*,2> inputs={&data,&mc};
 const std::array<std::string,2> samples={"data","mc"};
 double maximum=0;
 for(int layer=0;layer<30;++layer){
  const bool even=layer%2==0;
  const int nx=even?5:42,ny=even?42:5;
  const double xmax=even?113.5:111.3,ymax=even?111.3:113.5;
  masks[layer]=std::make_unique<TH2D>(Form("common_strip_mask_layer%02d",layer),";x [mm];y [mm];Common channel",nx,-xmax,xmax,ny,-ymax,ymax);
  masks[layer]->SetDirectory(nullptr);
  for(int s=0;s<2;++s){
   maps[s][layer]=std::make_unique<TH2D>(Form("%s_low_tail_rate_xy_layer%02d",samples[s].c_str(),layer),Form("%s, layer %d;x [mm];y [mm];Low-tail fraction",samples[s].c_str(),layer),nx,-xmax,xmax,ny,-ymax,ymax);
   maps[s][layer]->SetDirectory(nullptr);
  }
  std::set<int> occupied;
  for(int chip=0;chip<6;++chip)for(int channel=0;channel<(chip==5?30:36);++channel){
   const int id=layer*100000+chip*10000+channel;
   const auto* position=EBUdecode(layer,chip,channel);
   const int ix=masks[layer]->GetXaxis()->FindBin(position[0]);
   const int iy=masks[layer]->GetYaxis()->FindBin(position[1]);
   require(ix>=1&&ix<=nx&&iy>=1&&iy<=ny,"Strip position outside layer map");
   require(occupied.insert(masks[layer]->GetBin(ix,iy)).second,"Duplicate strip position");
   if(!common.count(id))continue;
   masks[layer]->SetBinContent(ix,iy,1);
   for(int s=0;s<2;++s){
    double rate=inputs[s]->GetBinContent(layer+1,chip*36+channel+1);
    require(std::isfinite(rate)&&rate>=0&&rate<=1,"Invalid tail fraction");
    maps[s][layer]->SetBinContent(ix,iy,rate);maximum=std::max(maximum,rate);
   }
  }
  require(occupied.size()==210,"Incomplete strip geometry");
 }
 maximum=std::max(maximum,1e-6);
 TParameter<double>("beam_GeV",beam).Write();
 TParameter<double>("color_max_fraction",maximum).Write();
 TNamed("geometry","calibration_parameters/EBUdecode.cxx; pitch 5.3 x 45.4 mm; even layers x/y swapped; gray = outside common channel set").Write();
 for(auto& mask:masks)mask->Write();
 for(int s=0;s<2;++s)for(int page=0;page<2;++page){
  TCanvas canvas(Form("tail_xy_%s_%d",samples[s].c_str(),page),"",2560,2560);
  canvas.SetCanvasSize(2560,2560);
  canvas.Divide(4,4,.001,.001);
  for(int pad=0;pad<16;++pad){
   canvas.cd(pad+1);const int layer=page*16+pad;
   if(layer>=30){
    TLatex text;text.SetNDC();text.SetTextSize(.043);
    text.DrawLatex(.12,.65,Form("%s, %g GeV",samples[s].c_str(),beam));
    text.DrawLatex(.12,.53,"Gray: outside common channels");
    text.DrawLatex(.12,.43,"White: zero low-tail fraction");
    text.DrawLatex(.12,.33,"Color scale shared by all layers / samples");
    break;
   }
   // Equal x/y extent and a square frame in pixels give identical mm scales.
   const double padWidth=canvas.GetWw()*gPad->GetAbsWNDC();
   const double padHeight=canvas.GetWh()*gPad->GetAbsHNDC();
   const double side=std::min(.66*padWidth,.66*padHeight);
   gPad->SetLeftMargin(.14);gPad->SetRightMargin(1-.14-side/padWidth);
   gPad->SetBottomMargin(.16);gPad->SetTopMargin(1-.16-side/padHeight);
   gPad->SetLogz(false);
   auto& h=*maps[s][layer];h.SetStats(false);h.SetMinimum(0);h.SetMaximum(maximum);
   h.SetTitle(Form("%s %g GeV, layer %d",samples[s].c_str(),beam,layer));
   for(auto* axis:{h.GetXaxis(),h.GetYaxis(),h.GetZaxis()}){axis->SetLabelSize(.037);axis->SetTitleSize(.042);axis->SetNdivisions(505);}
   h.GetZaxis()->SetTitleOffset(1.45);h.SetContour(100);
   // ROOT omits palette labels for an entirely zero histogram. Seed only
   // the disposable display copy; zero strips are painted white below.
   TH2D display(h);display.SetDirectory(nullptr);
   if(h.GetSumOfWeights()==0)display.SetBinContent(1,1,maximum*1e-12);
   auto* frame=gPad->DrawFrame(-120,-120,120,120,h.GetTitle());
   frame->GetXaxis()->SetTitle("x [mm]");frame->GetYaxis()->SetTitle("y [mm]");
   for(auto* axis:{frame->GetXaxis(),frame->GetYaxis()}){axis->SetLabelSize(.037);axis->SetTitleSize(.042);axis->SetNdivisions(505);}
   display.DrawCopy("COLZ SAME");
   for(int ix=1;ix<=h.GetNbinsX();++ix)for(int iy=1;iy<=h.GetNbinsY();++iy){
    TBox strip(h.GetXaxis()->GetBinLowEdge(ix),h.GetYaxis()->GetBinLowEdge(iy),h.GetXaxis()->GetBinUpEdge(ix),h.GetYaxis()->GetBinUpEdge(iy));
    strip.SetLineColor(kGray+1);strip.SetLineWidth(1);
    const bool valid=masks[layer]->GetBinContent(ix,iy)>0;
    strip.SetFillStyle(!valid||h.GetBinContent(ix,iy)==0?1001:0);strip.SetFillColor(valid?kWhite:kGray);
    strip.DrawClone();
   }
   gPad->RedrawAxis();file.cd();h.Write();
  }
  canvas.cd();canvas.Update();
  const std::string name="hg_lg_tail_layer_map_"+samples[s]+(page?"_layers16_29":"_layers00_15");
  auto png=out/"figures"/(name+".tmp."+std::to_string(getpid())+".png");
  canvas.SaveAs(png.c_str());require(fs::exists(png)&&fs::file_size(png)>0,"Layer PNG save failed");
  fs::rename(png,out/"figures"/(name+".png"));
 }
 require(!file.TestBit(TFile::kWriteError),"Layer ROOT write failed");file.Close();
 fs::rename(temporary,out/"hg_lg_tail_layer_maps.root");
}
