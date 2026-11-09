#include "EventHistograms.hh"
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
int main(int argc,char**argv)try{
 require(argc==2,"PlotResponse OUTPUT_DIRECTORY");fs::path out=argv[1];gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);
 for(auto&entry:fs::directory_iterator(out)){auto dir=entry.path();if(!entry.is_directory()||!fs::exists(dir/"response.root"))continue;auto f=openFile((dir/"response.root").string());TCanvas canvas("conditional","",1300,600);canvas.Divide(2,1);
  for(int region=0;region<2;++region){canvas.cd(region+1);gPad->SetLeftMargin(.16);gPad->SetBottomMargin(.13);std::array<TProfile*,3> p{};for(int s=0;s<3;++s)p[s]=object<TProfile>(*f,eventplots::samples[s]+(region?"_central_energy_vs_nhit":"_energy_vs_nhit"));double xmin=1e9,xmax=-1,ymin=1e9,ymax=-1;std::array<TGraphErrors,3> graphs;
   for(int b=1;b<=p[0]->GetNbinsX();++b){if(p[0]->GetBinEntries(b)<30||p[1]->GetBinEntries(b)<30)continue;double x=p[0]->GetBinCenter(b);xmin=std::min(xmin,x);xmax=std::max(xmax,x);for(int s=0;s<3;++s){int n=graphs[s].GetN();double y=p[s]->GetBinContent(b),error=p[s]->GetBinError(b);graphs[s].SetPoint(n,x,y);graphs[s].SetPointError(n,0,error);ymin=std::min(ymin,y-error);ymax=std::max(ymax,y+error);}}
   if(xmax<0)continue;
   double dx=std::max(p[0]->GetBinWidth(1),.1*(xmax-xmin)),dy=std::max(.001,.55*(ymax-ymin));gPad->DrawFrame(xmin-dx,std::max(0.,ymin-.2*dy),xmax+dx,ymax+dy,(dir.filename().string()+": shower + CoG20;"+(region?"Central positive good hits":"Positive good hits")+";Mean energy [GeV]").c_str());TLegend legend(.18,.74,.56,.89);legend.SetTextSize(.027);for(int s=0;s<3;++s){graphs[s].SetMarkerStyle(20+s);graphs[s].SetMarkerSize(.65);graphs[s].SetMarkerColor(eventplots::colors[s]);graphs[s].SetLineColor(eventplots::colors[s]);auto*drawn=graphs[s].DrawClone("PZ SAME");legend.AddEntry(drawn,eventplots::samples[s].c_str(),"lp");}legend.DrawClone();
  }canvas.SaveAs((dir/"figures/energy_at_fixed_nhit.png").c_str());
  std::ifstream table(dir/"layers.tsv");std::string line;std::getline(table,line);std::array<TGraph,3> energies;std::array<TGraph,2> hits;double em=0,nm=0;while(std::getline(table,line)){std::istringstream row(line);std::string sample;int layer;double n,e,z,h,central,ch;row>>sample>>layer>>n>>e>>z>>h>>central>>ch;int s=sample=="data"?0:1;energies[s].SetPoint(layer,layer,e);hits[s].SetPoint(layer,layer,h);if(s)energies[2].SetPoint(layer,layer,z);em=std::max(em,e);nm=std::max(nm,h);}
  canvas.Clear();canvas.Divide(2,1);for(int r=0;r<2;++r){canvas.cd(r+1);gPad->SetLeftMargin(.15);gPad->DrawFrame(-.5,0,29.5,(r?nm:em)*1.3,(dir.filename().string()+": shower + CoG20;Layer;"+(r?"Mean positive good hits":"Mean energy [GeV]")).c_str());TLegend leg(.6,.7,.89,.89);leg.SetTextSize(.025);for(int s=0;s<(r?2:3);++s){auto&g=r?hits[s]:energies[s];g.SetLineColor(eventplots::colors[s]);g.SetMarkerColor(eventplots::colors[s]);g.SetMarkerStyle(20+s);auto*drawn=g.DrawClone("PL SAME");leg.AddEntry(drawn,eventplots::samples[s].c_str(),"lp");}leg.DrawClone();}canvas.SaveAs((dir/"figures/layer_response.png").c_str());
  if(fs::exists(dir/"topology.root")){auto topology=openFile((dir/"topology.root").string());std::array<TH2D*,2> xy={object<TH2D>(*topology,"data_cog"),object<TH2D>(*topology,"mc_cog")};double max=0;for(auto*h:xy){if(h->GetEntries())h->Scale(1./h->GetEntries());max=std::max(max,h->GetMaximum());}canvas.Clear();canvas.Divide(2,1);for(int s=0;s<2;++s){canvas.cd(s+1);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.17);gPad->SetBottomMargin(.13);xy[s]->GetYaxis()->SetTitleOffset(1.6);xy[s]->GetYaxis()->CenterTitle();xy[s]->SetTitle((dir.filename().string()+": "+(s?"MC":"data")+", shower selection").c_str());xy[s]->SetMaximum(max);xy[s]->SetMinimum(0);xy[s]->Draw("COLZ");}canvas.SaveAs((dir/"figures/shower_position.png").c_str());}
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
