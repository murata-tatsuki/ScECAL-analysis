// Read-only rendering of retained joint-response histograms; no scientific tables rewritten.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <array>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace fs=std::filesystem;
int main(int argc,char**argv)try{
 if(argc!=2)throw std::runtime_error("Usage: JointResponseReplot JOINT_DIGITIZATION_DIR");fs::path out=argv[1];gROOT->SetBatch();gStyle->SetOptStat(0);TH1::AddDirectory(false);
 const std::array<int,7>colors={kBlack,kGreen+2,kBlue+1,kRed+1,kMagenta+1,kOrange+7,kCyan+2};
 const std::array<const char*,7>labels={"data, tail veto","original threshold","CT input baseline","norm only","norm + threshold","norm + pixels","norm + both"};
 for(int beam:{5,20,100}){std::vector<std::string>names={"data","original_threshold","CT_input_baseline","muon_norm"};if(beam==100)for(auto n:{"norm_threshold10","norm_pixels19600","norm_pixels19600_threshold10"})names.push_back(n);fs::path dest=out/(std::to_string(beam)+"GeV");TFile f((dest/"joint_response.root").c_str(),"READ");if(f.IsZombie())throw std::runtime_error("Missing retained ROOT");TCanvas canvas("joint_response","",2250,1500);canvas.SetCanvasSize(2250,1500);canvas.Divide(3,2);std::array<std::vector<std::unique_ptr<TH1D>>,6>hist;
  for(int p=0;p<6;++p){canvas.cd(p+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.03);gPad->SetBottomMargin(.13);double ymax=0;for(auto&name:names){std::string id="panel"+std::to_string(p)+"_"+name;TH1D*stored=nullptr;f.GetObject(id.c_str(),stored);if(!stored||stored->GetEntries()<=0)throw std::runtime_error("Missing stored histogram "+id);hist[p].emplace_back(static_cast<TH1D*>(stored->Clone()));auto&h=hist[p].back();h->SetDirectory(nullptr);h->Scale(1./h->GetEntries());ymax=std::max(ymax,h->GetMaximum());}TLegend legend(.42,.64,.92,.91);legend.SetTextSize(.027);for(size_t s=0;s<names.size();++s){auto&h=hist[p][s];h->SetLineColor(colors[s]);h->SetLineWidth(2);h->SetMaximum(1.8*ymax);h->Draw(s?"HIST SAME":"HIST");legend.AddEntry(h.get(),labels[s],"l");}legend.DrawClone();}
  canvas.SaveAs((dest/"joint_response.png").c_str());canvas.Clear();f.Close();
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
