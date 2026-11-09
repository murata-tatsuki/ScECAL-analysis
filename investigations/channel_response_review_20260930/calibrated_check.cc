#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
int main(int argc,char**argv){
 if(argc!=3)return 1;
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 gStyle->SetTextFont(42);gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleFont(42,"XYZ");
 gStyle->SetLabelSize(.035,"XYZ");gStyle->SetTitleSize(.04,"XYZ");
 TH2D* h[2];TProfile* p[2];long long nev[2]={},nall[2]={},n[2]={};double sum[2]={},sum2[2]={};
 std::ofstream out(std::string(argv[2])+"/channel8_cell240008_calibrated_hg_lg.tsv");
 out<<"sample\tevents\tchannel_hits\tband_min_EH_MeV\tband_max_EH_MeV\tband_hits\tmean_EL_minus_EH_MeV\tSEM_MeV\n"<<std::setprecision(10);
 for(int s=0;s<2;++s){
  std::string sample=s?"mc":"data";
  h[s]=new TH2D((sample+"_h").c_str(),(sample+" : 100 GeV, CellID 240008;E_{HG} [MeV];E_{LG} [MeV]").c_str(),100,0,2,130,-.3,2.3);
  p[s]=new TProfile((sample+"_p").c_str(),"100 GeV, CellID 240008;E_{HG} [MeV];Mean E_{LG} - E_{HG} [MeV]",40,0,2);
  std::ifstream manifest(std::string(argv[1])+"/inputs/100GeV_threshold_"+sample+".tsv");if(!manifest)return 2;
  std::string line;int files=0;
  while(std::getline(manifest,line)){
   auto sep=line.find('\t');if(sep==std::string::npos)return 3;TFile f(line.substr(sep+1).c_str());auto*t=f.Get<TTree>("Calib_Hit");if(!t)return 4;
   t->SetBranchStatus("*",0);std::vector<int>*id=nullptr;std::vector<double>*eh=nullptr,*el=nullptr;
   auto bind=[&](const char*key,auto*ptr){t->SetBranchStatus(key,1);if(t->SetBranchAddress(key,ptr)<0)throw std::runtime_error(key);};
   bind("CellID",&id);bind("Hit_HG_Energy",&eh);bind("Hit_LG_Energy",&el);
   for(Long64_t e=0;e<t->GetEntries();++e){t->GetEntry(e);++nev[s];if(id->size()!=eh->size()||id->size()!=el->size())return 5;
    for(size_t i=0;i<id->size();++i){int v=id->at(i);int cell=v/100000*100000+(v%100000/10000)*10000+v%100;if(cell!=240008)continue;
     double x=eh->at(i),y=el->at(i);if(!std::isfinite(x)||!std::isfinite(y))return 6;
     ++nall[s];h[s]->Fill(x,y);p[s]->Fill(x,y-x);
     if(x>=.3&&x<1.2){++n[s];sum[s]+=y-x;sum2[s]+=(y-x)*(y-x);}
    }
   }
   t->ResetBranchAddresses();delete id;delete eh;delete el;
   if(++files%10==0)std::cout<<sample<<" files="<<files<<std::endl;
  }
  double mean=sum[s]/n[s],sem=std::sqrt((sum2[s]-sum[s]*sum[s]/n[s])/(n[s]*(n[s]-1)));
  out<<sample<<'\t'<<nev[s]<<'\t'<<nall[s]<<"\t0.3\t1.2\t"<<n[s]<<'\t'<<mean<<'\t'<<sem<<'\n';
  std::cout<<sample<<" events="<<nev[s]<<" channel hits="<<nall[s]<<" band hits="<<n[s]<<" residual="<<mean<<" +- "<<sem<<std::endl;
 }
 TCanvas c("c","calibrated HG LG",1500,500);c.Divide(3,1);TLine diag(0,0,2,2),zero(0,0,2,0);diag.SetLineStyle(2);zero.SetLineStyle(2);
 for(int s=0;s<2;++s){c.cd(s+1);gPad->SetRightMargin(.15);gPad->SetLeftMargin(.13);gPad->SetLogz();h[s]->Draw("COLZ");diag.DrawClone();}
 c.cd(3);gPad->SetLeftMargin(.16);p[0]->SetLineColor(kBlack);p[0]->SetMarkerColor(kBlack);p[0]->SetMarkerStyle(20);p[1]->SetLineColor(kRed+1);p[1]->SetMarkerColor(kRed+1);p[1]->SetMarkerStyle(24);
 p[0]->SetMinimum(-.15);p[0]->SetMaximum(.35);p[0]->Draw("E1");p[1]->Draw("E1 SAME");zero.DrawClone();
 TLegend legend(.18,.72,.59,.87);legend.SetTextFont(42);legend.SetTextSize(.035);legend.AddEntry(p[0],"Data (mean +/- SEM)","lp");legend.AddEntry(p[1],"MC (mean +/- SEM)","lp");legend.Draw();
 c.SaveAs((std::string(argv[2])+"/channel8_cell240008_calibrated_hg_lg.png").c_str());
 TFile result((std::string(argv[2])+"/channel8_cell240008_calibrated_hg_lg.root").c_str(),"RECREATE");for(int s=0;s<2;++s){h[s]->Write();p[s]->Write();}c.Write();
}
