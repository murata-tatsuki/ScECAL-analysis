#pragma once
#include "Common.hh"
// Presentation only. Hybrid interval axes are laid out in IntervalPlots.hh.
inline void readableAxes(TCanvas& canvas) {
 canvas.SetEditable(true);
 if(canvas.GetListOfPrimitives()->FindObject("interval_axis"))return;
 const bool map=std::string(canvas.GetName()).find("strip_map")!=std::string::npos;
 TIter pads(canvas.GetListOfPrimitives());while(auto* obj=pads()) {
  auto* pad=dynamic_cast<TPad*>(obj);if(!pad)continue;
  bool color=false,frame=false;
  TIter probe(pad->GetListOfPrimitives());while(auto*p=probe())if(auto*h=dynamic_cast<TH1*>(p)){frame=true;color|=h->GetDimension()==2&&!map;}
  if(!frame)continue;
  pad->SetLeftMargin(map?.22:.18);pad->SetRightMargin(color?.24:.06);pad->SetBottomMargin(map?.20:.19);pad->SetTopMargin(.12);
  TIter primitives(pad->GetListOfPrimitives());while(auto*p=primitives())if(auto*h=dynamic_cast<TH1*>(p)) {
   for(auto*axis:{h->GetXaxis(),h->GetYaxis(),h->GetZaxis()}) {
    axis->SetLabelFont(42);axis->SetTitleFont(42);
    if(axis->GetLabelSize()!=0)axis->SetLabelSize(map?.065:.055);
    axis->SetTitleSize(map?.070:.062);axis->SetLabelOffset(.008);axis->CenterTitle(true);
   }
   h->GetXaxis()->SetTitleOffset(1.25);h->GetYaxis()->SetTitleOffset(1.30);
   if(map){h->GetXaxis()->SetNdivisions(505);h->GetYaxis()->SetNdivisions(505);}
   if(color){h->GetZaxis()->SetLabelSize(.045);h->GetZaxis()->SetTitleSize(.050);h->GetZaxis()->SetTitleOffset(1.45);}
  }
  pad->Modified();
 }
 canvas.Modified();canvas.Update();
}
