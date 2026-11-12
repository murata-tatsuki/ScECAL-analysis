#include "FitModel.hh"
#include "../channel_response/adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
using namespace mip;
namespace {
void drawSaved(TDirectory&d,const std::string&name,const Result&r){
 auto*h=dynamic_cast<TH1D*>(d.Get((name+"_used").c_str()));auto*g=dynamic_cast<TGraph*>(d.Get((name+"_expected").c_str()));
 if(!h)return;
 h->SetDirectory(nullptr);h->SetBit(kCanDelete);h->SetTitle(name.c_str());h->SetLineColor(kBlue+2);h->SetMaximum(h->GetMaximum()*1.45);h->GetXaxis()->SetRangeUser(0,mmr::group(r.CellID)?400:1200);h->Draw("E HIST");
 if(g){g->SetLineColor(kRed);g->SetLineWidth(2);g->Draw("L SAME");}
 auto*l=new TLegend(.39,.64,.93,.90);l->SetBit(kCanDelete);l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(.029);
 for(std::string t:{std::string(Form("MPV %.2f #pm %.2f",r.LandauMPV,r.LandauMPV_error)),std::string(Form("w %.2f, #sigma %.2f, dev %.1f",r.LandauWidth,r.GauSigma,r.objective)),r.reason})l->AddEntry((TObject*)nullptr,t.c_str(),"");
 l->Draw();
}
}
int main(int argc,char**argv)try{
 require(argc==5,"Usage: Diagnose histograms.root configuration.tsv output-directory CellID,...");
 gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPadLeftMargin(.13);
 auto source=mmr::openFile(argv[1]);mmr::Config cfg(argv[2]);std::string out=argv[3];auto selected=selection(argv[4]);mmr::Config runCfg(out+"/inputs/run_configuration.tsv");
 fs::create_directories(out+"/figures");fresh(out+"/diagnostics.root");TFile output((out+"/diagnostics.root").c_str(),"CREATE");
 auto meta=readMeta(*source);auto conditions=readConditions(*source);int mode=cfg.text("temperature-mode","current")=="current"?1:0;
 std::ofstream noise(out+"/noise.tsv");noise<<"CellID\tHG_pedestal_sigma\n";
 {auto pf=mmr::openFile(cfg.text("pedestal"));auto*t=mmr::tree(*pf,"ChnLevel");std::vector<int>*ids=nullptr;std::vector<double>*sd=nullptr;std::vector<float>*sf=nullptr;mmr::bind(t,"CellID",&ids);bool isFloat=std::string(t->GetBranch("PedHighSig")->GetClassName()).find("vector<float")!=std::string::npos;if(isFloat)mmr::bind(t,"PedHighSig",&sf);else mmr::bind(t,"PedHighSig",&sd);for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);for(size_t j=0;j<ids->size();++j)if(selected.count(ids->at(j)))noise<<ids->at(j)<<'\t'<<(isFloat?sf->at(j):sd->at(j))<<'\n';}t->ResetBranchAddresses();delete ids;delete sd;delete sf;}
 std::ofstream models(out+"/models.tsv"),profiles(out+"/profiles.tsv"),runs(out+"/runs.tsv");models<<"hypothesis\t";resultHeader(models);
 profiles<<"CellID\tmodel\tMPV_fixed\tdeviance\tdelta_deviance\tLandauWidth\tGauSigma\tstatus\n";
 runs<<"CellID\trun\tentries\tbelow20\t20to40\t40to60\t60to80\t80to100\tthreshold_shift_sigma\tMPV\tMPV_error\tdeviance\tstatus\n";
 std::vector<std::pair<std::string,FitSpec>> specs;
 FitSpec s;specs.push_back({"nominal",s});
 for(double low:{20.,40.,60.}){s=FitSpec{};s.rawLow=low;specs.push_back({"cut"+std::to_string(int(low)),s});}
 for(double sig:{2.,3.,4.}){s=FitSpec{};s.plateauSigmas=sig;specs.push_back({"plateau"+std::to_string(int(sig)),s});}
 for(double shift:{-1.,1.}){s=FitSpec{};s.rawLow=40;s.thresholdShift=shift;specs.push_back({shift<0?"cut40_threshold_minus":"cut40_threshold_plus",s});}
 s=FitSpec{};s.rawLow=40;s.sigmaScale=1.5;specs.push_back({"cut40_threshold_wide",s});
 for(double low:{0.,40.}){s=FitSpec{};s.rawLow=low;s.noGaussian=true;specs.push_back({low==0?"landau_nominal":"landau_cut40",s});}
 s=FitSpec{};s.mpvLower=-800;specs.push_back({"extended_MPV",s});
 for(int id:selected){auto&m=meta.at(id);std::vector<InputCondition>input;
  for(auto&c:conditions[id])if(c.mode==mode)input.push_back({c,getHist(*source,channelName(id)+"/conditions/"+c.name())});
  auto*d=directory(output,channelName(id));std::map<std::string,Result>results;
  for(auto&[name,spec]:specs){auto f=fitBest(input,m,cfg,spec);results[name]=f.r;models<<name<<'\t';resultRow(models,f.r);models.flush();saveFit(*d,f,name);std::cout<<id<<' '<<name<<" MPV="<<f.r.LandauMPV<<" dev="<<f.r.objective<<" status="<<f.r.status<<std::endl;}
  // Profiles compare only identical data intervals and models. Different cuts have different likelihood normalizations.
  for(const std::string name:{"nominal","cut40","plateau3"}){
   FitSpec spec;for(auto&v:specs)if(v.first==name)spec=v.second;spec.mpvLower=-800;spec.hesse=false;
   const double step=mmr::group(id)?10:35,lo=-10*step,hi=20*step;
   std::vector<double>x,y;std::vector<Result>points;double minimum=results.at(name).objective;
   std::set<double>grid;for(double mpv=lo;mpv<=hi;mpv+=step)grid.insert(mpv);
   const auto&center=results.at(name);
   if(std::isfinite(center.LandauMPV)){grid.insert(center.LandauMPV);if(std::isfinite(center.LandauMPV_error))for(double offset:{-2.,-1.,1.,2.})grid.insert(center.LandauMPV+offset*center.LandauMPV_error);}
   auto extended=fitBest(input,m,cfg,spec);
   if(extended.r.status==0&&std::isfinite(extended.r.objective)){minimum=std::min(minimum,extended.r.objective);grid.insert(extended.r.LandauMPV);}
   for(double mpv:grid){auto f=fitBest(input,m,cfg,spec,mpv);points.push_back(f.r);if(f.r.status==0&&std::isfinite(f.r.objective))minimum=std::min(minimum,f.r.objective);}
   for(auto&r:points){double delta=r.objective-minimum;profiles<<id<<'\t'<<name<<'\t'<<r.LandauMPV<<'\t'<<std::setprecision(12)<<r.objective<<'\t'<<delta<<'\t'<<r.LandauWidth<<'\t'<<r.GauSigma<<'\t'<<r.status<<'\n';if(r.status==0&&std::isfinite(delta)){x.push_back(r.LandauMPV);y.push_back(delta);}}
   profiles.flush();d->cd();TGraph g(x.size(),x.data(),y.data());g.SetName((name+"_profile").c_str());g.Write();
  }
  std::set<int>runNumbers;for(auto&i:input)runNumbers.insert(i.c.run);
  for(int run:runNumbers){std::vector<InputCondition>one;double bins[6]={};
   for(auto&i:input)if(i.c.run==run){auto h=std::unique_ptr<TH1D>(static_cast<TH1D*>(i.h->Clone()));h->SetDirectory(nullptr);for(int b=1;b<=h->GetNbinsX();++b){double raw=h->GetBinCenter(b)/i.c.factor,n=h->GetBinContent(b);bins[0]+=n;for(int k=1;k<=5;++k)if(raw<(20*k)&&(k==1||raw>=20*(k-1)))bins[k]+=n;}one.push_back({i.c,std::move(h)});}
   for(double shift:{-2.,-1.,0.,1.,2.}){FitSpec spec;spec.rawLow=40;spec.thresholdShift=shift;auto f=fitBest(one,m,runCfg,spec);runs<<id<<'\t'<<run;for(double b:bins)runs<<'\t'<<b;runs<<'\t'<<shift<<'\t'<<f.r.LandauMPV<<'\t'<<f.r.LandauMPV_error<<'\t'<<f.r.objective<<'\t'<<f.r.status<<'\n';}
  }runs.flush();
  TCanvas c("diagnosis","MPV/threshold diagnosis",2200,1600);c.Divide(3,3);
  int pad=0;for(const std::string name:{"nominal","cut20","cut40","plateau3","landau_cut40","extended_MPV"}){c.cd(++pad);drawSaved(*d,name,results.at(name));}
  for(const std::string name:{"nominal","cut40","plateau3"}){c.cd(++pad);auto*g=dynamic_cast<TGraph*>(d->Get((name+"_profile").c_str()));if(!g||g->GetN()==0)continue;g->SetTitle((name+";Fixed Landau MPV [ADC];Profile #Delta deviance").c_str());g->SetLineColor(kRed);g->SetMarkerColor(kRed);g->SetMarkerStyle(20);g->SetMinimum(0);g->SetMaximum(25);g->Draw("ALP");auto*l=new TLine(g->GetX()[0],3.84,g->GetX()[g->GetN()-1],3.84);l->SetBit(kCanDelete);l->SetLineStyle(2);l->Draw();}
  c.Update();fastPng(c,(out+"/figures/cell"+std::to_string(id)+"_diagnosis.png").c_str());
 }
 output.Write();output.Close();require(bool(models)&&bool(profiles)&&bool(runs),"Diagnostic output write failed");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
