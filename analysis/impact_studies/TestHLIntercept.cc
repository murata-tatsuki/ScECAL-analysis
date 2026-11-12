// Numerical fixtures: exact gain boundary, signed offsets, fallback constants,
// non-reference temperature, raw/calibrated hit reordering and empty events.
#include "StudyCommon.hh"
#include "EventHistograms.hh"
#include <TProfile.h>
#include <TParameter.h>
#include <iostream>
#include <array>
#include <algorithm>
void near(double a,double b){require(std::abs(a-b)<1e-10,"Numerical check failed: "+std::to_string(a)+" vs "+std::to_string(b));}
void fixtures(const fs::path&dir,bool corrupt){
 fs::create_directories(dir);
 {TFile f((dir/"hl.root").c_str(),"RECREATE");TTree t("InterCalib","");int id;double slope=.025,b; t.Branch("CellID",&id);t.Branch("Slope",&slope);t.Branch("Intercept",&b);id=420008;b=-2;t.Fill();id=240008;b=1;t.Fill();t.Write();}
 {TFile f((dir/"mip.root").c_str(),"RECREATE");TTree t("MIP_Fit","");int id,ndf=10;double mpv,chi=10;t.Branch("CellID",&id);t.Branch("LandauMPV",&mpv);t.Branch("ChiSquare",&chi);t.Branch("NDF",&ndf);id=420008;mpv=100;t.Fill();id=240008;mpv=200;t.Fill();t.Write();}
 {TFile f((dir/"ped.root").c_str(),"RECREATE");TTree t("ChnLevel","");std::vector<int>id={420008,420009,240008};std::vector<double>ph(3,400),pl(3,300);t.Branch("CellID",&id);t.Branch("PedHighMean",&ph);t.Branch("PedLowMean",&pl);t.Fill();t.Write();}
 for(bool calibrated:{false,true}){TFile f((dir/(calibrated?"cal.root":"raw.root")).c_str(),"RECREATE");TTree t(calibrated?"Calib_Hit":"Raw_Hit","");int run=7,time=100,trigger;t.Branch("Run_Num",&run);t.Branch("Event_Time",&time);t.Branch(calibrated?"Event_Num":"TriggerID",&trigger);
  std::vector<int>ids,tags;std::vector<double>h,l,e,temp;t.Branch("CellID",&ids);t.Branch(calibrated?"Hit_HG_Energy":"HG_Charge",&h);t.Branch(calibrated?"Hit_LG_Energy":"LG_Charge",&l);if(calibrated){t.Branch("Hit_Energy",&e);t.Branch("NewTemperature",&temp);}else t.Branch("HitTag",&tags);
  for(int event=0;event<4;++event){trigger=event;ids.clear();h.clear();l.clear();e.clear();tags.clear();temp.clear();std::vector<std::pair<int,double>> hits;
   if(event==0)hits={{420008,2599}};
   if(event==1)hits={{420108,2600},{240008,2700}};
   if(event==2)hits={{420009,4000},{420008,2600}};
   if(calibrated)std::reverse(hits.begin(),hits.end());
   for(auto[id,adc]:hits){ids.push_back(id);tags.push_back(1);double temperature=25;temp.push_back(temperature);double C=1-5*(group(id)?1.6/135:3.5/230),mpv=group(id)?100:200,b=group(id)?-2:1;if(id==420009)b=0;double eh=.305/mpv*C*(adc-400),el=.305/mpv*(C*20/.025-b/.025);
    if(corrupt&&calibrated&&event==1)el+=1;
    h.push_back(calibrated?eh:adc);l.push_back(calibrated?el:320);e.push_back(adc<2600?eh:el);
   }t.Fill();
  }t.Write();
 }
}
int main(int argc,char**argv)try{
 require(argc==3,"TestHLIntercept fixtures|bad|check|legacy|binning DIRECTORY");fs::path dir=argv[2];std::string mode=argv[1];
 if(mode=="legacy"){TFile f((dir/"study.root").c_str(),"UPDATE");require(!f.IsZombie(),"Cannot edit fixture");f.Delete("data_events;*");return 0;}
 if(mode=="binning"){
  eventplots::Values values;
  for(int i=0;i<10000;++i){values[0].push_back(3.9+.4*std::sin(i));values[1].push_back(4.13+.22*std::sin(i));values[2].push_back(values[1].back()-.087);}
  values[0][0]=-1000;values[0][1]=1000;eventplots::Histograms hist(values);
  require(hist.totals[1]->GetBinWidth(1)<.05,"MC peak collapsed into coarse bins");
  require(hist.totals[0]->GetBinContent(0)==1&&hist.totals[0]->GetBinContent(hist.totals[0]->GetNbinsX()+1)==1,"Tail flow counts lost");
  for(int s=0;s<3;++s){near(hist.fullTotals[s]->Integral(),values[s].size());near(hist.totals[s]->GetXaxis()->GetXmin(),hist.totals[0]->GetXaxis()->GetXmin());}
  eventplots::Values emptyEvents;for(auto&v:emptyEvents)v.assign(10,0);eventplots::Histograms zero(emptyEvents);
  near(zero.change->Integral(),10);near(zero.relative->GetEntries(),0);near(zero.totals[0]->Integral(),10);
  std::cout<<"Adaptive histogram checks passed\n";return 0;
 }
 if(mode=="fixtures"||mode=="bad"){fixtures(dir,mode=="bad");return 0;}
 require(mode=="check","Unknown mode");auto f=openFile((dir/"study.root").string());auto*t=object<TTree>(*f,"paired_events");require(t->GetEntries()==4,"Empty events lost");double delta;int nlg;bind(t,"delta_GeV",&delta);bind(t,"n_lg",&nlg);std::array<double,4> expected={0,-.000183,-.000244,0};std::array<int,4>counts={0,2,2,0};for(int i=0;i<4;++i){t->GetEntry(i);near(delta,expected[i]);require(nlg==counts[i],"HG boundary or LG count wrong");}
 near(object<TParameter<Long64_t>>(*f,"validated_mc_hits")->GetVal(),5);
 auto*n=object<TProfile>(*f,"channel_420008/mc_nominal_lg_minus_hg_vs_raw_hg");auto*z=object<TProfile>(*f,"channel_420008/mc_intercept_zero_lg_minus_hg_vs_raw_hg");int bin=n->FindBin(2600);near(n->GetBinContent(bin)-z->GetBinContent(bin),.244);
 require(fs::file_size(dir/"figures/event_energy.png")>100,"Missing plot");std::cout<<"Numerical checks passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
