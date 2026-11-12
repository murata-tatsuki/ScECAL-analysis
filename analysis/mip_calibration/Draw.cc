#include "Common.hh"
#include "../channel_response/adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLatex.h>
namespace {
using namespace mip;
void title(TLegend&l,const std::string&s){l.AddEntry((TObject*)nullptr,s.c_str(),"");}
void style(TH1*h,int color){h->SetLineColor(color);h->SetLineWidth(3);h->GetXaxis()->SetTitle("Pedestal-subtracted HG ADC (2 mm reference)");h->GetYaxis()->SetTitle("Acquired crossings / bin");h->GetXaxis()->SetTitleSize(.042);h->GetYaxis()->SetTitleSize(.042);}
void drawFit(TFile&fit,TDirectory&dir,const Result&r,const std::string&prefix,const char*label,TH1D*fallback=nullptr){
 auto*source=dynamic_cast<TH1D*>(dir.Get((prefix+"_used").c_str()));auto*g=dynamic_cast<TGraph*>(dir.Get((prefix+"_expected").c_str()));
 if(!source&&!fallback){TLatex text;text.SetNDC();text.SetTextSize(.035);text.DrawLatex(.12,.8,(std::string(label)+": no fit / insufficient statistics").c_str());return;}
 auto*h=static_cast<TH1D*>((source?source:fallback)->Clone());h->SetDirectory(nullptr);h->SetBit(kCanDelete);
 style(h,kBlue+2);h->SetTitle(label);h->SetMaximum(std::max(1.,h->GetMaximum()*1.42));h->Draw("E HIST");if(g){g->SetLineColor(kRed);g->SetLineWidth(3);g->SetLineStyle(2);g->Draw("L SAME");}
 auto*legend=new TLegend(.43,.62,.94,.91);legend->SetBit(kCanDelete);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.029);legend->AddEntry(h,"Data muon HG ADC","l");if(g)legend->AddEntry(g,prefix=="nominal"?"Landau-Gaussian #times efficiency":"Landau-Gaussian (plateau)","l");if(source){title(*legend,Form("MPV %.2f #pm %.2f; w_{L} %.2f; #sigma_{G} %.2f",r.LandauMPV,r.LandauMPV_error,r.LandauWidth,r.GauSigma));title(*legend,Form("#chi^{2}/NDF %.2f / %d = %.2f; convolved peak %.2f",r.ChiSquare,r.NDF,r.NDF>0?r.ChiSquare/r.NDF:mip::nan,r.convolution_peak));}else title(*legend,"No fitted curve");title(*legend,r.reason);legend->Draw();(void)fit;
}
// Reconstruct the saved fit input for channels without a nominal_used histogram.
// This is only a sum/rebin of stored run/factor histograms; no fitting is done.
std::unique_ptr<TH1D> unfittedHistogram(TFile&hf,int id,const Result&r,
                                     const std::vector<Condition>&conditions,int mode){
 const double width=std::isfinite(r.bin_width)&&r.bin_width>0?r.bin_width:(mmr::group(id)?4:8);
 const double low=std::isfinite(r.fit_low)?r.fit_low:0;
 const double high=std::isfinite(r.fit_high)&&r.fit_high>low?r.fit_high:(mmr::group(id)?1000:2000);
 auto h=hist(channelName(id)+"_unfitted",std::max(1,int(std::lround((high-low)/width))),low,high);
 for(const auto&c:conditions)if(c.mode==mode){
  auto input=getHist(hf,channelName(id)+"/conditions/"+c.name());
  for(int b=1;b<=input->GetNbinsX();++b){double x=input->GetBinCenter(b);if(x>=low&&x<high)h->AddBinContent(h->FindFixBin(x),input->GetBinContent(b));}
 }
 for(int b=1;b<=h->GetNbinsX();++b)h->SetBinError(b,std::sqrt(h->GetBinContent(b)));
 h->SetEntries(h->Integral());return h;
}
void drawChips(TFile&hf,TFile&ff,const std::map<int,Result>&results,
               const std::set<int>&selected,const std::string&out){
 auto conditions=readConditions(hf);mmr::Config config(out+"/inputs/configuration.tsv");
 const int mode=config.text("temperature-mode","current")=="current"?1:0;
 std::set<std::pair<int,int>>chips;for(int id:selected)chips.emplace(id/100000,id%100000/10000);
 // Match the actual pixel size of each pad in the existing 2250x1500, 2x2 canvas.
 // ROOT subtracts window decorations from that constructor size; use its drawable size.
 TCanvas reference("chip_size_reference","",2250,1500);
 const int width=3*reference.GetWw(),height=3*reference.GetWh();reference.Close();
 for(const auto&[layer,chip]:chips){
  const std::string path=out+"/figures/layer"+std::to_string(layer);
  fs::create_directories(path);TCanvas canvas("chip","Nominal fits by channel",width,height);
  canvas.SetCanvasSize(width,height);canvas.Divide(6,6,.01/3,.01/3);
  for(int channel=0;channel<36;++channel){
   canvas.cd(channel+1);const int id=layer*100000+chip*10000+channel;
   if(chip==5&&channel>=30){TLatex text;text.SetNDC();text.SetTextSize(.035);text.DrawLatex(.12,.8,Form("Channel %d: not connected",channel));continue;}
   const auto&r=results.at(id);auto*d=ff.GetDirectory(channelName(id).c_str());
   const std::string label=Form("Cell %d / channel %d: nominal efficiency fit",id,channel);
   std::unique_ptr<TH1D>fallback;
   if(!d||!d->GetKey("nominal_used"))fallback=unfittedHistogram(hf,id,r,conditions[id],mode);
   drawFit(ff,d?*d:static_cast<TDirectory&>(ff),r,"nominal",label.c_str(),fallback.get());
  }
  canvas.Update();const std::string file=path+"/layer"+std::to_string(layer)+"_chip"+std::to_string(chip)+".png";
  fastPng(canvas,file.c_str());require(fs::is_regular_file(file)&&fs::file_size(file)>0,"Cannot write "+file);
  std::cout<<"Chip plot: "<<file<<" ("<<width<<"x"<<height<<", "<<(chip==5?30:36)<<" channels)"<<std::endl;
 }
}
}
int main(int argc,char**argv)try{
 using namespace mip;require(argc==5||(argc==6&&std::string(argv[5])=="--chip-plots-only"),"Usage: Draw histograms.root fits.root output-directory all|CellID,... [--chip-plots-only]");gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);gStyle->SetPadLeftMargin(.12);gStyle->SetPadRightMargin(.08);
 auto hf=mmr::openFile(argv[1]),ff=mmr::openFile(argv[2]);auto meta=readMeta(*hf);auto results=readResults(*ff);auto selected=selection(argv[4]);std::string out=argv[3];fs::create_directories(out+"/figures");
 drawChips(*hf,*ff,results,selected,out);if(argc==6)return 0;
 std::map<std::pair<int,int>,Result>variants;auto*vt=mmr::tree(*ff,"FitVariants");Result vr;vr.branches(*vt,true);std::string*why=nullptr;mmr::bind(vt,"reason",&why);for(Long64_t i=0;i<vt->GetEntries();++i){mmr::read(vt,i);vr.reason=*why;variants[{vr.CellID,vr.variant}]=vr;}vt->ResetBranchAddresses();delete why;
 TH2D mpv("mpv","Accepted Landau MPV;layer;chip #times 36 + channel",30,0,30,216,0,216),error("error","Accepted MPV relative error;layer;chip #times 36 + channel",30,0,30,216,0,216),status("status","Decision: 0 undecided, 1 accepted, 2 rejected;layer;chip #times 36 + channel",30,0,30,216,0,216);
 TH1D mpv10("mpv10","Accepted MPV;ADC;Channels",100,0,1000),mpv15("mpv15","Accepted MPV;ADC;Channels",100,0,2000);
 for(auto&[id,r]:results){int l=id/100000,c=id%100000/10000,ch=id%100;status.SetBinContent(l+1,c*36+ch+1,r.decision);if(r.decision==1){mpv.SetBinContent(l+1,c*36+ch+1,r.LandauMPV);error.SetBinContent(l+1,c*36+ch+1,r.relative_error);(mmr::group(id)?mpv10:mpv15).Fill(r.LandauMPV);}
 if(!selected.count(id)||(std::string(argv[4])=="all"&&meta.at(id).usable==0))continue;
 std::string path=out+"/figures/layer"+std::to_string(l)+"/chip"+std::to_string(c);fs::create_directories(path);
 TCanvas canvas("channel","Data muon MIP",2250,1500);canvas.Divide(2,2);canvas.cd(1);
 auto before=getHist(*hf,channelName(id)+"/path_adc"),after=getHist(*hf,channelName(id)+"/current_adc");before=rebinForPlot(*before,mmr::group(id)?2:4);after=rebinForPlot(*after,mmr::group(id)?2:4);style(before.get(),kBlue+2);style(after.get(),kAzure+7);after->SetLineStyle(2);before->SetTitle(Form("Cell %d: no temperature correction / current formula",id));before->GetXaxis()->SetRangeUser(0,mmr::group(id)?1000:2000);before->SetMaximum(std::max(before->GetMaximum(),after->GetMaximum())*1.35);before->Draw("HIST");after->Draw("HIST SAME");TLegend legend(.43,.72,.93,.9);legend.SetTextSize(.031);legend.SetBorderSize(0);legend.AddEntry(before.get(),"2 mm; no temperature correction","l");legend.AddEntry(after.get(),"2 mm; current temperature formula","l");legend.Draw();
 auto*d=ff->GetDirectory(channelName(id).c_str());canvas.cd(2);if(d)drawFit(*ff,*d,r,"nominal","Nominal efficiency fit");else {TLatex txt;txt.SetNDC();txt.DrawLatex(.1,.8,r.reason.c_str());}
 canvas.cd(3);if(d&&variants.count({id,2}))drawFit(*ff,*d,variants.at({id,2}),"variant2","Diagnostic: above threshold + 3 sigma");
 canvas.cd(4);TLatex text;text.SetNDC();text.SetTextSize(.034);auto&m=meta.at(id);std::vector<std::string>lines={Form("Cell %d; layer %d / chip %d / channel %d; %d #mum",id,l,c,ch,mmr::group(id)?10:15),"Decision: "+r.reason,Form("Predicted %lld; acquired %lld; unique %lld; usable %lld",m.predicted,m.observed,m.unique,m.usable),Form("Conditional acquisition efficiency: %.4f",m.predicted?double(m.observed)/m.predicted:mip::nan),Form("Threshold %.2f; turn-on width %.2f ADC (ped. sub.)",m.threshold,m.threshold_sigma),Form("Relative MPV error %.4f; max fit variation %.4f",r.relative_error,r.stability),Form("Old %s MIP: %.2f ADC",m.old_fallback?"mean fallback":"individual",m.old_mip),"MPV = Landau parameter; convolved maximum is separate.","Current temperature formula is a reader convention.","Crosstalk retained; no electron / MC scale constraint."};double y=.9;for(auto&s:lines){text.DrawLatex(.07,y,s.c_str());y-=.073;}
 canvas.Update();fastPng(canvas,(path+"/cell"+std::to_string(id)+"_fit.png").c_str());
 }
 // Maps always represent all saved channels, even in targeted redraw mode.
 TCanvas maps("maps","MIP channel maps",2250,1500);maps.Divide(3,1);maps.cd(1);mpv.Draw("COLZ");maps.cd(2);error.Draw("COLZ");maps.cd(3);status.SetMinimum(0);status.SetMaximum(2);status.Draw("COLZ");maps.Update();fastPng(maps,(out+"/figures/channel_maps.png").c_str());
 TCanvas distributions("distributions","Accepted MPVs",2250,1500);distributions.Divide(2,1);distributions.cd(1);mpv10.SetLineColor(kBlue+2);mpv10.Draw("HIST");distributions.cd(2);mpv15.SetLineColor(kAzure+7);mpv15.Draw("HIST");distributions.Update();fastPng(distributions,(out+"/figures/MPV_distributions.png").c_str());return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
