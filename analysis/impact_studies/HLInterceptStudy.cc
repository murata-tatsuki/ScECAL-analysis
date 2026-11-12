// Paired MC variation: remove only the additive HL intercept correction.
// No input ROOT file is modified. Saved calibrated hits define the selection.
#include "StudyCommon.hh"
#include "EventHistograms.hh"
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <array>
const std::array<std::string,3> names={"data","mc_nominal","mc_intercept_zero"};
const std::array<int,3> colors={kBlack,kBlue+1,kRed+1};
struct Channel {
 std::array<std::unique_ptr<TProfile>,3> residual,energy;
 std::array<Stats,3> overlap;
 Stats eventCorrection; // actually per LG hit; explicitly labelled in output
 long long mcHits=0,lgHits=0;
 Channel(int id){for(int s=0;s<3;++s){std::string n=names[s]+"_"+std::to_string(id);
  residual[s]=std::make_unique<TProfile>((n+"_residual").c_str(),";Raw HG ADC;Mean E_{LG} - E_{HG} [MeV]",84,0,4200);
  energy[s]=std::make_unique<TProfile>((n+"_selected").c_str(),";Raw HG ADC;Mean selected hit energy [MeV]",84,0,4200);
 }}
};
void drawProfiles(std::array<std::unique_ptr<TProfile>,3>&p,const std::string&title,double xmin,double xmax,bool zero){
 double lo=0,hi=0;for(auto&h:p)for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=xmin&&h->GetBinCenter(b)<=xmax&&h->GetBinEntries(b)){lo=std::min(lo,h->GetBinContent(b)-h->GetBinError(b));hi=std::max(hi,h->GetBinContent(b)+h->GetBinError(b));}
 double d=std::max(.01,hi-lo);
 gPad->DrawFrame(xmin,lo-.15*d,xmax,hi+.25*d,(title+";Raw HG ADC;"+p[0]->GetYaxis()->GetTitle()).c_str());
 TLegend leg(.15,zero?.71:.15,.57,zero?.9:.34);leg.SetTextFont(42);leg.SetTextSize(.028);
 for(int s=0;s<3;++s){TGraphErrors graph;auto*h=p[s].get();for(int b=1;b<=h->GetNbinsX();++b){double x=h->GetBinCenter(b);if(x<xmin||x>xmax||h->GetBinEntries(b)==0)continue;int n=graph.GetN();graph.SetPoint(n,x,h->GetBinContent(b));graph.SetPointError(n,0,h->GetBinError(b));}graph.SetLineColor(colors[s]);graph.SetMarkerColor(colors[s]);graph.SetMarkerStyle(20+s);graph.SetMarkerSize(.65);auto*drawn=graph.DrawClone("PZ SAME");leg.AddEntry(drawn,names[s].c_str(),"lp");}
 leg.DrawClone();
 if(zero){TLine line(xmin,0,xmax,0);line.SetLineStyle(2);line.DrawClone();}TLine cut(2600,lo-.15*d,2600,hi+.25*d);cut.SetLineStyle(3);cut.DrawClone();
}
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"mc-manifest",""},{"data-manifest",""},{"hl",""},{"mip",""},{"pedestal",""},{"out",""},{"energy","100"},{"max-events","0"},{"plot-channels","all"},{"plots","1"}});
 const double beam=std::stod(o["energy"]);const long long maxEvents=std::stoll(o["max-events"]);require(beam>0&&std::isfinite(beam)&&maxEvents>=0,"Invalid energy/max-events");
 fs::path out=o["out"];require(!out.empty(),"Output directory required");fs::create_directories(out/"figures");require(!fs::exists(out/"study.root"),"Output exists; choose a new --label");
 auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);
 std::set<int> plotIds;if(o["plot-channels"]!="all"){std::istringstream ss(o["plot-channels"]);std::string v;while(std::getline(ss,v,','))plotIds.insert(std::stoi(v));}
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetTextFont(42);gStyle->SetTitleFont(42,"XYZ");gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleSize(.04,"XYZ");gStyle->SetLabelSize(.035,"XYZ");
 TFile result((out/"study.root.tmp").c_str(),"RECREATE");require(!result.IsZombie(),"Cannot create output");
 eventplots::Values eventValues;std::array<Stats,3> stats;Stats delta,fraction;long long lgHitTotal=0;double maxBaselineError=0;long long validated=0;
 std::map<int,Channel> channels;
 TTree events("paired_events","Same MC events, all retained connected hits in layers 0-29");
 Long64_t eventEntry=0;int fileIndex=0,run=0,eventTime=0,trigger=0,nhit=0,nlg=0;double nominal=0,zero=0,correction=0;
 events.Branch("file_index",&fileIndex);events.Branch("entry",&eventEntry);events.Branch("Run_Num",&run);events.Branch("Event_Time",&eventTime);events.Branch("TriggerID",&trigger);events.Branch("nhit",&nhit);events.Branch("n_lg",&nlg);events.Branch("nominal_GeV",&nominal);events.Branch("intercept_zero_GeV",&zero);events.Branch("delta_GeV",&correction);
 for(int sample=0;sample<2;++sample){bool mc=sample==1;std::ifstream manifest(o[mc?"mc-manifest":"data-manifest"]);require(bool(manifest),"Cannot read manifest");std::string line;long long processed=0;fileIndex=-1;
  while(std::getline(manifest,line)){
   if(maxEvents&&processed>=maxEvents)break;
   ++fileIndex;size_t sep=line.find('\t');require(sep!=std::string::npos,"Manifest needs raw TAB calib");
   auto raw=openFile(line.substr(0,sep)),cal=openFile(line.substr(sep+1));auto*r=object<TTree>(*raw,"Raw_Hit"),*c=object<TTree>(*cal,"Calib_Hit");require(r->GetEntries()==c->GetEntries(),"Raw/calib entry counts differ");r->SetBranchStatus("*",0);c->SetBranchStatus("*",0);
   std::vector<int>*rid=nullptr,*tag=nullptr,*cid=nullptr;std::vector<double>*rh=nullptr,*rl=nullptr,*eh=nullptr,*el=nullptr,*en=nullptr,*temp=nullptr;
   int rr,rt,re,cr,ct,ce;bind(r,"Run_Num",&rr);bind(r,"Event_Time",&rt);bind(r,"TriggerID",&re);bind(c,"Run_Num",&cr);bind(c,"Event_Time",&ct);bind(c,"Event_Num",&ce);
   bind(r,"CellID",&rid);bind(r,"HitTag",&tag);bind(r,"HG_Charge",&rh);if(mc)bind(r,"LG_Charge",&rl);
   bind(c,"CellID",&cid);bind(c,"Hit_HG_Energy",&eh);bind(c,"Hit_LG_Energy",&el);bind(c,"Hit_Energy",&en);if(mc)bind(c,"NewTemperature",&temp);
   for(eventEntry=0;eventEntry<c->GetEntries()&&(!maxEvents||processed<maxEvents);++eventEntry){
    r->GetEntry(eventEntry);c->GetEntry(eventEntry);require(rr==cr&&rt==ct&&re==ce,"Raw/calib event order mismatch");require(rid->size()==tag->size()&&rid->size()==rh->size()&&(!mc||rid->size()==rl->size()),"Raw vector size mismatch");require(cid->size()==eh->size()&&cid->size()==el->size()&&cid->size()==en->size()&&(!mc||cid->size()==temp->size()),"Calib vector size mismatch");
    std::map<int,size_t> index;for(size_t j=0;j<rid->size();++j)if(tag->at(j)==1)require(index.emplace(rid->at(j),j).second,"Duplicate full CellID in raw event");
    std::set<int> used;nominal=0;zero=0;nhit=0;nlg=0;run=rr;eventTime=rt;trigger=re;
    for(size_t j=0;j<cid->size();++j){int full=cid->at(j),id=physical(full);if(!connected(id))continue;require(used.insert(full).second,"Duplicate full CellID in calib event");require(index.count(full),"Calibrated hit missing from raw event");size_t k=index.at(full);double adc=rh->at(k),h=eh->at(j),l=el->at(j),e=en->at(j);require(std::isfinite(adc)&&std::isfinite(h)&&std::isfinite(l)&&std::isfinite(e),"Nonfinite hit");bool useLG=adc>=2600;
     auto near=[&](double a,double b){double err=std::abs(a-b);maxBaselineError=std::max(maxBaselineError,err);require(std::isfinite(b)&&err<=1e-5+1e-6*std::abs(b),"Baseline calibration mismatch at CellID "+std::to_string(id)+" saved="+std::to_string(a)+" expected="+std::to_string(b));};
     near(e,useLG?l:h);auto [it,inserted]=channels.try_emplace(id,id);(void)inserted;auto&ch=it->second;int s=mc?1:0;
     ch.residual[s]->Fill(adc,l-h);ch.energy[s]->Fill(adc,e);if(adc>=800&&adc<2200)ch.overlap[s].add(l-h);
     nominal+=e/1000.;++nhit;nlg+=useLG;
     if(mc){const auto&p=pars.at(id);require(p.ped,"Missing pedestal for MC hit");double t=temp->at(j),C=1-(t-20)*(group(id)?1.6/135:3.5/230),factor=.305/p.mip;
      try { near(h,factor*C*(adc-p.ph));near(l,factor*(C*(rl->at(k)-p.pl)*p.gain+p.intercept)); }
      catch(const std::exception& error){std::cerr<<"raw HG="<<adc<<" LG="<<rl->at(k)<<" pedestal="<<p.ph<<","<<p.pl<<" mip="<<p.mip<<" gain="<<p.gain<<" intercept="<<p.intercept<<" temperature="<<t<<" saved HG="<<h<<" LG="<<l<<std::endl;throw;}
      ++validated;
      double lzero=l-p.offset(),ezero=useLG?lzero:h;zero+=ezero/1000.;ch.residual[2]->Fill(adc,lzero-h);ch.energy[2]->Fill(adc,ezero);if(adc>=800&&adc<2200)ch.overlap[2].add(lzero-h);++ch.mcHits;if(useLG){++ch.lgHits;ch.eventCorrection.add(ezero-e);}
     }
    }
    if(mc){stats[1].add(nominal);stats[2].add(zero);eventValues[1].push_back(nominal);eventValues[2].push_back(zero);correction=zero-nominal;delta.add(correction);if(nominal!=0)fraction.add(correction/nominal);lgHitTotal+=nlg;events.Fill();}
    else {stats[0].add(nominal);eventValues[0].push_back(nominal);}++processed;
   }
   r->ResetBranchAddresses();c->ResetBranchAddresses();delete rid;delete tag;delete cid;delete rh;delete rl;delete eh;delete el;delete en;delete temp;
   std::cout<<(mc?"MC":"data")<<" files="<<fileIndex+1<<" events="<<processed<<std::endl;
  }
  require(processed>0,"No events processed");
 }
 result.cd();events.Write();
 TTree dataEvents("data_events","Data event totals with the same retained hit selection");double dataEnergy=0;dataEvents.Branch("energy_GeV",&dataEnergy);for(double x:eventValues[0]){dataEnergy=x;dataEvents.Fill();}dataEvents.Write();
 eventplots::Histograms histograms(eventValues);auto&totals=histograms.totals;histograms.write();histograms.metadata(out);
 TParameter<double>("beam_energy_GeV",beam).Write();TParameter<double>("max_baseline_error_MeV",maxBaselineError).Write();TParameter<Long64_t>("validated_mc_hits",validated).Write();
 TNamed("variation","E_LG_zero=E_LG_nominal-0.305*interFactorIntercept/mipFactor; raw HG>=2600 selects LG; no additional event cuts").Write();
 std::ofstream summary(out/"summary.tsv");summary<<std::setprecision(12)<<"sample\tevents\tmean_GeV\tsem_GeV\trms_GeV\trms_over_mean\tunderflow\toverflow\n";
 for(int s=0;s<3;++s)summary<<names[s]<<'\t'<<stats[s].n<<'\t'<<stats[s].mean<<'\t'<<stats[s].sem()<<'\t'<<stats[s].rms()<<'\t'<<stats[s].rms()/stats[s].mean<<'\t'<<totals[s]->GetBinContent(0)<<'\t'<<totals[s]->GetBinContent(totals[s]->GetNbinsX()+1)<<'\n';
 std::ofstream paired(out/"paired_effect.tsv");paired<<std::setprecision(12)<<"events\tmean_delta_GeV\tsem_delta_GeV\trms_delta_GeV\tmean_event_relative_delta\tmean_relative_events\tlg_hits\tvalidated_mc_hits\tmax_baseline_error_MeV\n"<<delta.n<<'\t'<<delta.mean<<'\t'<<delta.sem()<<'\t'<<delta.rms()<<'\t'<<fraction.mean<<'\t'<<fraction.n<<'\t'<<lgHitTotal<<'\t'<<validated<<'\t'<<maxBaselineError<<'\n';
 std::ofstream table(out/"channels.tsv");table<<std::setprecision(12)<<"cell\toffset_MeV\tmip_fallback\thl_fallback\tmc_hits\tlg_hits\ttotal_delta_MeV\tdata_overlap_n\tdata_residual_MeV\tdata_sem\tmc_overlap_n\tmc_residual_MeV\tmc_sem\tzero_overlap_n\tzero_residual_MeV\tzero_sem\n";
 TH2D map("channel_mean_event_delta",";Layer;Chip * 36 + channel;Mean event energy change [MeV]",30,-.5,29.5,216,-.5,215.5);
 for(auto&[id,ch]:channels){const auto&p=pars.at(id);table<<id<<'\t'<<p.offset()<<'\t'<<p.mipFallback<<'\t'<<p.hlFallback<<'\t'<<ch.mcHits<<'\t'<<ch.lgHits<<'\t'<<ch.eventCorrection.mean*ch.eventCorrection.n;for(auto&s:ch.overlap)table<<'\t'<<s.n<<'\t'<<s.mean<<'\t'<<s.sem();table<<'\n';
  map.SetBinContent(id/100000+1,id%100000/10000*36+id%100+1,ch.eventCorrection.mean*ch.eventCorrection.n/stats[1].n);
  result.cd();auto*d=result.mkdir(("channel_"+std::to_string(id)).c_str());d->cd();for(int s=0;s<3;++s){ch.residual[s]->Write((names[s]+"_lg_minus_hg_vs_raw_hg").c_str());ch.energy[s]->Write((names[s]+"_selected_energy_vs_raw_hg").c_str());}
  if(o["plots"]=="1"&&(o["plot-channels"]=="all"||plotIds.count(id))){fs::path dir=out/"figures"/"channels"/("layer"+std::to_string(id/100000))/("chip"+std::to_string(id%100000/10000));fs::create_directories(dir);TCanvas canvas("channel","",1300,550);canvas.Divide(2,1);canvas.cd(1);gPad->SetLeftMargin(.14);drawProfiles(ch.residual,std::to_string(id)+": HG-LG, "+o["energy"]+" GeV",600,2800,true);canvas.cd(2);gPad->SetLeftMargin(.14);drawProfiles(ch.energy,std::to_string(id)+": gain switch, "+o["energy"]+" GeV",2200,3000,false);canvas.SaveAs((dir/("channel"+std::to_string(id%100)+"_cell"+std::to_string(id)+"_gain_consistency.png")).c_str());}
 }
 result.cd();map.Write();
 if(o["plots"]=="1"){
  histograms.draw(out,beam);
  TCanvas canvas("overview","",1100,650);canvas.SetRightMargin(.18);map.Draw("COLZ");canvas.SaveAs((out/"figures/channel_energy_change.png").c_str());
 }
 result.Close();fs::rename(out/"study.root.tmp",out/"study.root");std::cout<<"Finished: "<<out<<"\nMean paired change [GeV]: "<<delta.mean<<" +/- "<<delta.sem()<<std::endl;
 return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
