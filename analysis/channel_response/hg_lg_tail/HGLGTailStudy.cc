// Joint hit/event analysis. Original raw/calibrated trees are read once.
// A temporary ROOT cache permits exact channel medians and one shared tail rule.
#include "StudyCommon.hh"
#include "HGLGTailLayerMaps.hh"
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <set>
#include <unordered_map>

using Pair = std::pair<std::string,std::string>;
std::vector<Pair> manifest(const std::string& path) {
 std::ifstream in(path); require(bool(in),"Cannot read manifest: "+path);
 std::vector<Pair> result; std::string line;
 while(std::getline(in,line)) { auto p=line.find('\t'); require(p!=std::string::npos,"Expected raw TAB calib");result.emplace_back(line.substr(0,p),line.substr(p+1)); }
 require(!result.empty(),"Empty manifest");return result;
}
std::set<int> readMask(const std::string& path) {
 if(path=="none") return {};
 std::ifstream in(path); require(bool(in),"Cannot read channel mask");std::string text((std::istreambuf_iterator<char>(in)),{});
 auto a=text.find("const unordered_set<int> bad_channels = {");require(a!=std::string::npos,"Missing bad_channels definition");a=text.find('{',a);auto b=text.find("};",a);require(b!=std::string::npos,"Invalid channel mask");text=text.substr(a,b-a);
 std::set<int> result;std::regex number("[0-9]+");for(std::sregex_iterator it(text.begin(),text.end(),number),end;it!=end;++it)result.insert(std::stoi(it->str()));return result;
}
double quantile(const std::vector<double>& v,double q) {
 require(!v.empty(),"Empty quantile");double k=q*(v.size()-1);size_t i=k;return v[i]+(k-i)*(v[std::min(i+1,v.size()-1)]-v[i]);
}
struct Distribution {
 std::vector<double> values;Stats stats;double center=0,width68=0;long long low=0,high=0;
};
struct Channel {
 Parameter p;double sigma=0;bool good=false,common=false;std::array<Distribution,2> d;
};
struct CacheEvent {
 int sample=0,file=0,run=0,time=0,trigger=0,bcid=-1;
 Long64_t entry=0;double energy=0,lgEnergy=0,negativeLG=0;
 int goodHits=0,lgHits=0,negativeLGHits=0;bool cog=false;
 std::vector<int> ids,hitBCID,gainTag;
 std::vector<double> residual,hg,lg,eh,el;
 void scalars(TTree& t) {
  t.Branch("sample",&sample);t.Branch("file_index",&file);t.Branch("entry",&entry);t.Branch("Run_Num",&run);t.Branch("Event_Time",&time);t.Branch("TriggerID",&trigger);t.Branch("modal_BCID",&bcid);
  t.Branch("energy_GeV",&energy);t.Branch("lg_energy_GeV",&lgEnergy);t.Branch("negative_lg_energy_MeV",&negativeLG);
  t.Branch("good_hits",&goodHits);t.Branch("lg_hits",&lgHits);t.Branch("negative_lg_hits",&negativeLGHits);t.Branch("cog20",&cog);
 }
 void branches(TTree& t) {
  scalars(t);t.Branch("ids",&ids);t.Branch("residual",&residual);t.Branch("hg",&hg);t.Branch("lg",&lg);t.Branch("eh",&eh);t.Branch("el",&el);t.Branch("hitBCID",&hitBCID);t.Branch("gainTag",&gainTag);
 }
 void clear() {energy=lgEnergy=negativeLG=0;goodHits=lgHits=negativeLGHits=0;cog=false;ids.clear();residual.clear();hg.clear();lg.clear();eh.clear();el.clear();hitBCID.clear();gainTag.clear();}
};
struct Population {Stats energy;long long low=0,high=0,overlap=0,negativeLG=0;};
struct SelectedEvent {double energy;int lowChannels,negativeLG;};
int main(int argc,char**argv)try {
 auto o=options(argc,argv,{{"data-manifest",""},{"mc-manifest",""},{"pedestal",""},{"hl",""},{"mip",""},{"bad-channel-source","none"},{"out",""},{"scratch",""},{"energy","100"},{"max-events","0"},{"events-per-file","0"},{"max-files","0"},{"hg-min","800"},{"hg-max","2200"},{"nsigma","5"},{"min-hits","100"},{"min-tail-channels","3"},{"selection","shower-cog20"},{"reference-energy","0"}});
 const long long maxEvents=std::stoll(o["max-events"]),perFile=std::stoll(o["events-per-file"]),maxFiles=std::stoll(o["max-files"]),minHits=std::stoll(o["min-hits"]);
 const int minTail=std::stoi(o["min-tail-channels"]);const double hmin=std::stod(o["hg-min"]),hmax=std::stod(o["hg-max"]),nsigma=std::stod(o["nsigma"]),beam=std::stod(o["energy"]);
 require(maxEvents>=0&&perFile>=0&&maxFiles>=0&&minHits>=2&&minTail>=1,"Invalid count option");require(std::isfinite(hmin)&&std::isfinite(hmax)&&std::isfinite(nsigma)&&std::isfinite(beam)&&hmin>=0&&hmax>hmin&&hmax<=2600&&nsigma>0&&beam>0,"Invalid ADC/significance/beam range");
 require(o["selection"]=="all"||o["selection"]=="cog20"||o["selection"]=="shower-cog20","Unknown event selection");
 fs::path out=o["out"];require(!out.empty()&&!o["scratch"].empty(),"Output and scratch required");require(!fs::exists(out/"tail_study.root"),"Completed output exists");fs::create_directories(out/"figures");
 auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);auto excluded=readMask(o["bad-channel-source"]);std::map<int,Channel> channels;
 for(auto&[id,p]:pars){channels[id].p=p;channels[id].good=p.mip>0&&!excluded.count(id);}
 {auto f=openFile(o["pedestal"]);auto*t=object<TTree>(*f,"ChnLevel");std::vector<int>*id=nullptr;std::vector<double>*sh=nullptr,*sl=nullptr;bind(t,"CellID",&id);bind(t,"PedHighSig",&sh);bind(t,"PedLowSig",&sl);
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);require(id->size()==sh->size()&&id->size()==sl->size(),"Pedestal vector size mismatch");for(size_t j=0;j<id->size();++j){int cell=physical(id->at(j));if(!channels.count(cell))continue;auto&c=channels.at(cell);require(std::isfinite(sh->at(j))&&std::isfinite(sl->at(j))&&sh->at(j)>0&&sl->at(j)>0,"Invalid pedestal noise");c.sigma=std::sqrt(sl->at(j)*sl->at(j)+std::pow(sh->at(j)/c.p.gain,2)+1./12);}}
  t->ResetBranchAddresses();delete id;delete sh;delete sl;
 }
 for(auto&[id,c]:channels)if(c.good)require(c.p.ped&&c.sigma>0,"Missing pedestal: "+std::to_string(id));
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 TFile cacheFile((fs::path(o["scratch"])/"events.root").c_str(),"RECREATE","",1);require(!cacheFile.IsZombie(),"Cannot create cache");TTree cache("cache","Single raw/calib read, before tail classification");cache.SetAutoFlush(-16*1024*1024);CacheEvent event;event.branches(cache);
 std::vector<double> mcEnergy;std::array<long long,2> checked{},processed{};std::array<double,2> maxH{},maxL{};
 std::ofstream samples(out/"inputs/sampled_files.tsv");samples<<"sample\tfile_index\tentries\tsampled\traw\tcalib\n";
 for(int s=0;s<2;++s){auto files=manifest(o[s?"mc-manifest":"data-manifest"]);if(maxFiles&&files.size()>size_t(maxFiles))files.resize(maxFiles);
  for(size_t fi=0;fi<files.size();++fi){if(maxEvents&&processed[s]>=maxEvents)break;
   auto raw=openFile(files[fi].first),cal=openFile(files[fi].second);auto*r=object<TTree>(*raw,"Raw_Hit"),*c=object<TTree>(*cal,"Calib_Hit");require(r->GetEntries()==c->GetEntries(),"Raw/calib event counts differ");r->SetBranchStatus("*",0);c->SetBranchStatus("*",0);
   int rr,rt,re,cr,ct,ce;std::vector<int>*ri=nullptr,*tags=nullptr,*gain=nullptr,*bcid=nullptr,*ci=nullptr;std::vector<double>*rh=nullptr,*rl=nullptr,*eh=nullptr,*el=nullptr,*en=nullptr,*temp=nullptr,*x=nullptr,*y=nullptr;
   bind(r,"Run_Num",&rr);bind(r,"Event_Time",&rt);bind(r,"TriggerID",&re);bind(r,"CellID",&ri);bind(r,"HitTag",&tags);bind(r,"GainTag",&gain);bind(r,"BCID",&bcid);bind(r,"HG_Charge",&rh);bind(r,"LG_Charge",&rl);
   bind(c,"Run_Num",&cr);bind(c,"Event_Time",&ct);bind(c,"Event_Num",&ce);bind(c,"CellID",&ci);bind(c,"Hit_HG_Energy",&eh);bind(c,"Hit_LG_Energy",&el);bind(c,"Hit_Energy",&en);bind(c,"NewTemperature",&temp);bind(c,"Hit_X",&x);bind(c,"Hit_Y",&y);
   r->SetCacheSize(16*1024*1024);c->SetCacheSize(16*1024*1024);
   for(const char*k:{"Run_Num","Event_Time","TriggerID","CellID","HitTag","GainTag","BCID","HG_Charge","LG_Charge"})r->AddBranchToCache(k,true);
   for(const char*k:{"Run_Num","Event_Time","Event_Num","CellID","Hit_HG_Energy","Hit_LG_Energy","Hit_Energy","NewTemperature","Hit_X","Hit_Y"})c->AddBranchToCache(k,true);
   r->StopCacheLearningPhase();c->StopCacheLearningPhase();
   Long64_t count=c->GetEntries(),take=perFile?std::min<Long64_t>(perFile,count):count;if(maxEvents)take=std::min<Long64_t>(take,maxEvents-processed[s]);
   int blocks=perFile?std::min<Long64_t>(10,take):(take?1:0);
   for(int b=0;b<blocks;++b){Long64_t start=perFile?count*b/blocks:0,end=perFile?count*(b+1)/blocks:count,n=take/blocks+(b<take%blocks);if(perFile)start+=(end-start-n)/2;
    for(Long64_t entry=start;entry<start+n;++entry){require(r->GetEntry(entry)>0&&c->GetEntry(entry)>0,"Failed event read");require(rr==cr&&rt==ct&&re==ce,"Raw/calib event key mismatch");
     require(ri->size()==rh->size()&&ri->size()==rl->size()&&ri->size()==tags->size()&&ri->size()==gain->size()&&ri->size()==bcid->size(),"Raw vector sizes differ");require(ci->size()==eh->size()&&ci->size()==el->size()&&ci->size()==en->size()&&ci->size()==temp->size()&&ci->size()==x->size()&&ci->size()==y->size(),"Calibrated vector sizes differ");
     event.clear();event.sample=s;event.file=fi;event.entry=entry;event.run=rr;event.time=rt;event.trigger=re;std::unordered_map<int,size_t> index;index.reserve(ri->size());std::map<int,int> times;
     for(size_t k=0;k<ri->size();++k)if(tags->at(k)==1){require(index.emplace(ri->at(k),k).second,"Duplicate tagged full CellID");++times[bcid->at(k)];}
     event.bcid=-1;int maxBC=0;for(auto&[t,nBC]:times)if(nBC>maxBC){maxBC=nBC;event.bcid=t;}
     std::set<int> used;std::array<double,30> le{},lx{},ly{};
     for(size_t j=0;j<ci->size();++j){int full=ci->at(j),id=physical(full);if(!connected(id))continue;require(used.insert(full).second,"Duplicate calibrated full CellID");require(index.count(full),"Calibrated hit missing in raw");size_t k=index.at(full);auto&ch=channels.at(id);const auto&p=ch.p;
      double h=rh->at(k),l=rl->at(k),H=h-p.ph,L=l-p.pl,C=1-(temp->at(j)-20)*(group(id)?1.6/135:3.5/230),f=.305/p.mip;
      double expectedH=f*C*H,expectedL=f*(C*L*p.gain+p.intercept),errH=std::abs(eh->at(j)-expectedH),errL=std::abs(el->at(j)-expectedL);
      require(std::isfinite(expectedH)&&std::isfinite(expectedL)&&std::isfinite(en->at(j))&&std::isfinite(eh->at(j))&&std::isfinite(el->at(j))&&std::isfinite(x->at(j))&&std::isfinite(y->at(j)),"Nonfinite hit");require(errH<=1e-5+1e-6*std::abs(expectedH)&&errL<=1e-5+1e-6*std::abs(expectedL),"Calibration closure failed: "+std::to_string(id));require(std::abs(en->at(j)-(h>=2600?el->at(j):eh->at(j)))<1e-9,"Unexpected gain selection");++checked[s];maxH[s]=std::max(maxH[s],errH);maxL[s]=std::max(maxL[s],errL);
      if(!ch.good)continue;
      if(h>=2600){++event.lgHits;if(en->at(j)<0){++event.negativeLGHits;event.negativeLG+=en->at(j);}else event.lgEnergy+=en->at(j)/1000;}
      if(en->at(j)>=0){event.energy+=en->at(j)/1000;++event.goodHits;int layer=id/100000;le[layer]+=en->at(j);lx[layer]+=x->at(j)*en->at(j);ly[layer]+=y->at(j)*en->at(j);}
      if(h<hmin||h>=hmax)continue;
      // -b = B/A. Median centering removes any constant data/MC offset,
      // including the missing digitization intercept in the old MC.
      double residual=L-H/p.gain+p.intercept/p.gain;
      ch.d[s].values.push_back(residual);ch.d[s].stats.add(residual);
      event.ids.push_back(full);event.residual.push_back(residual);event.hg.push_back(h);event.lg.push_back(l);event.eh.push_back(eh->at(j));event.el.push_back(el->at(j));event.hitBCID.push_back(bcid->at(k));event.gainTag.push_back(gain->at(k));
     }
     event.cog=true;for(int l:{9,10})event.cog=event.cog&&le[l]>0&&std::abs(lx[l]/le[l])<20&&std::abs(ly[l]/le[l])<20;
     cache.Fill();++processed[s];if(s)mcEnergy.push_back(event.energy);
    }
   }
   samples<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<count<<'\t'<<take<<'\t'<<files[fi].first<<'\t'<<files[fi].second<<'\n';
   r->ResetBranchAddresses();c->ResetBranchAddresses();delete ri;delete tags;delete gain;delete bcid;delete ci;delete rh;delete rl;delete eh;delete el;delete en;delete temp;delete x;delete y;
   std::cout<<(s?"mc":"data")<<" files="<<fi+1<<" processed events="<<processed[s]<<std::endl;
  }
  require(processed[s]>0,"No sampled events");
  for(auto&[id,ch]:channels){auto&d=ch.d[s];if(!d.values.empty()){std::sort(d.values.begin(),d.values.end());d.center=quantile(d.values,.5);d.width68=(quantile(d.values,.84)-quantile(d.values,.16))/2;}std::vector<double>().swap(d.values);}
 }
 cacheFile.cd();cache.Write();cacheFile.Flush();
 int common=0;for(auto&[id,ch]:channels){ch.common=ch.good&&ch.d[0].stats.n>=minHits&&ch.d[1].stats.n>=minHits;common+=ch.common;}require(common>0,"No common channels; increase statistics or reduce --min-hits");
 std::sort(mcEnergy.begin(),mcEnergy.end());double reference=std::stod(o["reference-energy"]);require(std::isfinite(reference)&&reference>=0,"Invalid reference energy");if(reference==0)reference=quantile(mcEnergy,.5);require(reference>0,"MC energy reference is zero");std::vector<double>().swap(mcEnergy);
 TFile result((out/"tail_study.root.tmp").c_str(),"RECREATE");require(!result.IsZombie(),"Cannot create result");
 bool selected=false,tailRich=false;int overlapAll=0,overlap=0,nlow=0,nhigh=0,nlowChannels=0,maxMultiplicity=0;Long64_t cacheIndex=0;
 std::vector<int> lowCells,highCells;
 TTree events("events","Same tail definition and common channel set as channel statistics");event.scalars(events);events.Branch("cache_index",&cacheIndex);events.Branch("selected",&selected);events.Branch("tail_rich",&tailRich);events.Branch("overlap_all_good",&overlapAll);events.Branch("overlap_common",&overlap);events.Branch("low_tail_hits",&nlow);events.Branch("high_tail_hits",&nhigh);events.Branch("low_tail_channels",&nlowChannels);events.Branch("low_tail_cellid",&lowCells);events.Branch("high_tail_cellid",&highCells);events.SetAutoFlush(-16*1024*1024);
 int tailCell=0,tailBCID=0,tailGain=0;double tailResidual=0,z=0,h=0,l=0,ehValue=0,elValue=0;bool low=false;
 TTree tailHits("tail_hits","All low and high tails in common channels, with event identity");tailHits.Branch("sample",&event.sample);tailHits.Branch("file_index",&event.file);tailHits.Branch("entry",&event.entry);tailHits.Branch("Run_Num",&event.run);tailHits.Branch("TriggerID",&event.trigger);tailHits.Branch("event_index",&cacheIndex);tailHits.Branch("CellID",&tailCell);tailHits.Branch("BCID",&tailBCID);tailHits.Branch("GainTag",&tailGain);tailHits.Branch("residual_ADC",&tailResidual);tailHits.Branch("centered_z",&z);tailHits.Branch("HG_raw",&h);tailHits.Branch("LG_raw",&l);tailHits.Branch("EH_MeV",&ehValue);tailHits.Branch("EL_MeV",&elValue);tailHits.Branch("low",&low);tailHits.SetAutoFlush(-16*1024*1024);
 std::ofstream table(out/"events.tsv");table<<std::setprecision(12)<<"sample\tfile_index\tentry\trun\tevent_time\ttrigger\tmodal_BCID\tenergy_GeV\tselected\toverlap_all_good\toverlap_common\tlow_tail_hits\thigh_tail_hits\tlow_tail_channels\ttail_rich\tLG_hits\tnegative_LG_hits\tnegative_LG_energy_MeV\n";
 std::array<std::array<Population,6>,2> pop;std::array<std::vector<double>,3> energies;std::array<long long,2> totalsLow{},totalsHigh{},totalsOverlap{};
 std::array<std::vector<SelectedEvent>,2> selectedEvents;
 std::array<std::unique_ptr<TH1D>,2> residualPlots,multi;std::array<std::unique_ptr<TH2D>,2> maps;
 for(int s=0;s<2;++s){std::string name=s?"mc":"data";residualPlots[s]=std::make_unique<TH1D>((name+"_residual_z").c_str(),";Median-centered HG/LG residual / expected noise;Weighted probability / bin",200,-25,25);residualPlots[s]->Sumw2();multi[s]=std::make_unique<TH1D>((name+"_tail_channels").c_str(),";Distinct low-tail channels / selected event;Probability / bin",6301,-.5,6300.5);maps[s]=std::make_unique<TH2D>((name+"_low_tail_rate").c_str(),(name+";Layer;Chip * 36 + channel").c_str(),30,-.5,29.5,216,-.5,215.5);}
 for(cacheIndex=0;cacheIndex<cache.GetEntries();++cacheIndex){require(cache.GetEntry(cacheIndex)>0,"Cannot read cached event");overlapAll=event.ids.size();overlap=nlow=nhigh=0;lowCells.clear();highCells.clear();std::set<int> distinct;
  int s=event.sample;
  for(size_t j=0;j<event.ids.size();++j){tailCell=event.ids[j];auto&ch=channels.at(physical(tailCell));if(!ch.common)continue;++overlap;tailResidual=event.residual[j];z=(tailResidual-ch.d[s].center)/ch.sigma;
   double weight=double(std::min(ch.d[0].stats.n,ch.d[1].stats.n))/ch.d[s].stats.n;residualPlots[s]->Fill(z,weight);
   low=z < -nsigma;bool high=z>nsigma;
   if(low){++nlow;++ch.d[s].low;lowCells.push_back(tailCell);distinct.insert(physical(tailCell));}
   if(high){++nhigh;++ch.d[s].high;highCells.push_back(tailCell);}
   if(low||high){h=event.hg[j];l=event.lg[j];ehValue=event.eh[j];elValue=event.el[j];tailBCID=event.hitBCID[j];tailGain=event.gainTag[j];tailHits.Fill();}
  }
  nlowChannels=distinct.size();tailRich=nlowChannels>=minTail;selected=o["selection"]=="all"||(event.cog&&(o["selection"]=="cog20"||(event.energy>.5*reference&&event.energy<1.5*reference)));
  totalsLow[s]+=nlow;totalsHigh[s]+=nhigh;totalsOverlap[s]+=overlap;
  for(int k=0;k<6;++k)if((k<3||selected)&&(k%3==0||(k%3==1&&tailRich)||(k%3==2&&!tailRich))){auto&p=pop[s][k];p.energy.add(event.energy);p.low+=nlow;p.high+=nhigh;p.overlap+=overlap;p.negativeLG+=event.negativeLGHits;}
  if(selected){maxMultiplicity=std::max(maxMultiplicity,nlowChannels);multi[s]->Fill(nlowChannels);energies[s].push_back(event.energy);selectedEvents[s].push_back({event.energy,nlowChannels,event.negativeLGHits});if(!s&&!tailRich)energies[2].push_back(event.energy);}
  events.Fill();table<<(s?"mc":"data")<<'\t'<<event.file<<'\t'<<event.entry<<'\t'<<event.run<<'\t'<<event.time<<'\t'<<event.trigger<<'\t'<<event.bcid<<'\t'<<event.energy<<'\t'<<selected<<'\t'<<overlapAll<<'\t'<<overlap<<'\t'<<nlow<<'\t'<<nhigh<<'\t'<<nlowChannels<<'\t'<<tailRich<<'\t'<<event.lgHits<<'\t'<<event.negativeLGHits<<'\t'<<event.negativeLG<<'\n';
 }
 std::ofstream channelTable(out/"channels.tsv");channelTable<<std::setprecision(12)<<"cell\tgood\tcommon\tmip\tHL_A\tHL_B\tnoise_sigma_LG_ADC\tdata_n\tdata_median\tdata_width68\tdata_RMS\tdata_low\tdata_high\tmc_n\tmc_median\tmc_width68\tmc_RMS\tmc_low\tmc_high\n";
 std::array<long long,2> channelLow{},channelHigh{},channelOverlap{};std::array<double,2> weightedLow{},weightedHigh{};double totalWeight=0;
 for(auto&[id,ch]:channels){channelTable<<id<<'\t'<<ch.good<<'\t'<<ch.common<<'\t'<<ch.p.mip<<'\t'<<ch.p.gain<<'\t'<<ch.p.intercept<<'\t'<<ch.sigma;for(int s=0;s<2;++s){auto&d=ch.d[s];channelTable<<'\t'<<d.stats.n<<'\t'<<d.center<<'\t'<<d.width68<<'\t'<<d.stats.rms()<<'\t'<<d.low<<'\t'<<d.high;if(ch.common){channelLow[s]+=d.low;channelHigh[s]+=d.high;channelOverlap[s]+=d.stats.n;double w=std::min(ch.d[0].stats.n,ch.d[1].stats.n);weightedLow[s]+=w*d.low/d.stats.n;weightedHigh[s]+=w*d.high/d.stats.n;maps[s]->SetBinContent(id/100000+1,id%100000/10000*36+id%100+1,double(d.low)/d.stats.n);}}channelTable<<'\n';if(ch.common)totalWeight+=std::min(ch.d[0].stats.n,ch.d[1].stats.n);}
 std::ofstream validation(out/"validation.tsv");validation<<std::setprecision(12)<<"sample\tevents\tcalibration_checked_hits\tmax_HG_error_MeV\tmax_LG_error_MeV\tevent_low_sum\tchannel_low_sum\tevent_high_sum\tchannel_high_sum\tevent_overlap_sum\tchannel_overlap_sum\n";
 for(int s=0;s<2;++s){require(totalsLow[s]==channelLow[s]&&totalsHigh[s]==channelHigh[s]&&totalsOverlap[s]==channelOverlap[s],"Hit/event tail bookkeeping mismatch");validation<<(s?"mc":"data")<<'\t'<<processed[s]<<'\t'<<checked[s]<<'\t'<<maxH[s]<<'\t'<<maxL[s]<<'\t'<<totalsLow[s]<<'\t'<<channelLow[s]<<'\t'<<totalsHigh[s]<<'\t'<<channelHigh[s]<<'\t'<<totalsOverlap[s]<<'\t'<<channelOverlap[s]<<'\n';}
 std::ofstream summary(out/"summary.tsv");summary<<std::setprecision(12)<<"sample\tselection\tevents\tmean_GeV\tRMS_GeV\tRMS_over_mean\tSEM_GeV\toverlap_hits\tlow_tail_hits\thigh_tail_hits\tnegative_LG_hits\n";const std::array<std::string,6> names={"all","all_tail_rich","all_tail_poor","selected","selected_tail_rich","selected_tail_poor"};
 for(int s=0;s<2;++s)for(int k=0;k<6;++k){auto&p=pop[s][k];summary<<(s?"mc":"data")<<'\t'<<names[k]<<'\t'<<p.energy.n<<'\t'<<p.energy.mean<<'\t'<<p.energy.rms()<<'\t'<<(p.energy.mean?p.energy.rms()/p.energy.mean:0)<<'\t'<<p.energy.sem()<<'\t'<<p.overlap<<'\t'<<p.low<<'\t'<<p.high<<'\t'<<p.negativeLG<<'\n';}
 std::ofstream rate(out/"tail_rates.tsv");rate<<std::setprecision(12)<<"sample\tcommon_channels\toverlap_hits\tlow_tail_hits\thigh_tail_hits\tmatched_weight_low_fraction\tmatched_weight_high_fraction\n";for(int s=0;s<2;++s)rate<<(s?"mc":"data")<<'\t'<<common<<'\t'<<totalsOverlap[s]<<'\t'<<totalsLow[s]<<'\t'<<totalsHigh[s]<<'\t'<<weightedLow[s]/totalWeight<<'\t'<<weightedHigh[s]/totalWeight<<'\n';
 result.cd();events.Write();tailHits.Write();TParameter<double>("beam_GeV",beam).Write();TParameter<double>("reference_energy_GeV",reference).Write();TParameter<double>("nsigma",nsigma).Write();TParameter<int>("min_tail_channels",minTail).Write();TParameter<int>("common_channels",common).Write();TNamed("event_selection",o["selection"].c_str()).Write();TNamed("tail_definition","Same for hit/event: common channels with >=min-hits per sample; z=(LG-pedLG-(HG-pedHG)/A+B/A-channel_sample_median)/sqrt(sigLG^2+(sigHG/A)^2+1/12); low iff z<-nsigma; event flag counts distinct physical channels.").Write();
 TCanvas canvas("overview","",1300,550);canvas.Divide(2,1);
 for(int panel=0;panel<2;++panel){canvas.cd(panel+1);gPad->SetLogy();double max=0;for(int s=0;s<2;++s){auto*h=panel?multi[s].get():residualPlots[s].get();h->Write();double norm=h->Integral(0,h->GetNbinsX()+1);if(norm>0)h->Scale(1./norm);max=std::max(max,h->GetMaximum());}TLegend leg(.6,.73,.89,.89);for(int s=0;s<2;++s){auto*h=panel?multi[s].get():residualPlots[s].get();h->SetLineColor(s?kBlue+1:kBlack);h->SetLineWidth(2);h->SetMinimum(1e-7);h->SetMaximum(std::max(1e-6,max*2));if(panel)h->GetXaxis()->SetRangeUser(-.5,std::max(10.,maxMultiplicity*1.1));h->Draw(s?"HIST SAME":"HIST");leg.AddEntry(h,s?"MC":"data","l");}leg.DrawClone();}
 canvas.SaveAs((out/"figures/hg_lg_tail_overview.png").c_str());
 canvas.Clear();canvas.Divide(2,1);double mapMax=std::max(maps[0]->GetMaximum(),maps[1]->GetMaximum());for(int s=0;s<2;++s){canvas.cd(s+1);gPad->SetRightMargin(.16);maps[s]->SetMaximum(std::max(mapMax,1e-6));maps[s]->Write();maps[s]->Draw("COLZ");}canvas.SaveAs((out/"figures/hg_lg_tail_channel_map.png").c_str());
 {std::set<int> commonCells;for(const auto& [id,ch]:channels)if(ch.common)commonCells.insert(id);drawTailLayerMaps(*maps[0],*maps[1],commonCells,out,beam);}
 std::ofstream ranges(out/"histogram_ranges.tsv");ranges<<"view\tsample\tbins\tmin_GeV\tmax_GeV\tunderflow\toverflow\n";
 bool any=false;double lo=1e99,hi=-1,fullLo=1e99,fullHi=-1,width=1e99;for(auto&v:energies)if(!v.empty()){any=true;std::sort(v.begin(),v.end());lo=std::min(lo,quantile(v,.001));hi=std::max(hi,quantile(v,.999));fullLo=std::min(fullLo,v.front());fullHi=std::max(fullHi,v.back());double w=(quantile(v,.84)-quantile(v,.16))/20;if(w>0)width=std::min(width,w);}
 if(any){double pad=std::max(.0001,(hi-lo)*.05);lo=std::max(0.,lo-pad);hi+=pad;pad=std::max(.0001,(fullHi-fullLo)*.05);fullLo=std::max(0.,fullLo-pad);fullHi+=pad;
  for(int view=0;view<2;++view){double xmin=view?fullLo:lo,xmax=view?fullHi:hi;int bins=view?400:(width<1e90?std::clamp(int(std::ceil((xmax-xmin)/width)),80,800):160);canvas.Clear();canvas.SetLogy(view);std::array<std::unique_ptr<TH1D>,3> hist;double max=0;const std::array<std::string,3> labels={"data selected","MC selected","data selected, tail-poor"};
   for(int i=0;i<3;++i){hist[i]=std::make_unique<TH1D>(("event_energy_"+std::to_string(view)+"_"+std::to_string(i)).c_str(),(std::to_string(beam)+" GeV;Positive good event energy [GeV];Probability / bin").c_str(),bins,xmin,xmax);for(double e:energies[i])hist[i]->Fill(e);hist[i]->Write();ranges<<(view?"full":"main")<<'\t'<<i<<'\t'<<bins<<'\t'<<xmin<<'\t'<<xmax<<'\t'<<hist[i]->GetBinContent(0)<<'\t'<<hist[i]->GetBinContent(bins+1)<<'\n';if(!energies[i].empty())hist[i]->Scale(1./energies[i].size());max=std::max(max,hist[i]->GetMaximum());}
   TLegend leg(.58,.72,.89,.89);for(int i=0;i<3;++i){hist[i]->SetLineColor(i==0?kBlack:i==1?kBlue+1:kRed+1);hist[i]->SetLineWidth(2);hist[i]->SetMaximum(std::max(1e-6,max*(view?3:1.3)));hist[i]->SetMinimum(view?1e-7:0);hist[i]->Draw(i?"HIST SAME":"HIST");leg.AddEntry(hist[i].get(),labels[i].c_str(),"l");}leg.DrawClone();canvas.SaveAs((out/"figures"/(view?"event_energy_full_range.png":"event_energy.png")).c_str());
  }
 }
 // Plot the event association directly, using common axes for data and MC.
 if(any){
  int maxNegative=0;for(const auto& sample:selectedEvents)for(const auto&e:sample)maxNegative=std::max(maxNegative,e.negativeLG);
  const int nx=std::min(200,std::max(5,maxMultiplicity)+1),ny=std::min(200,std::max(5,maxNegative)+1);
  const double xmax=std::max(5,maxMultiplicity)+.5,ymax=std::max(5,maxNegative)+.5;
  TCanvas correlations("correlations","",1300,1000);correlations.Divide(2,2);
  std::array<std::unique_ptr<TH2D>,4> plots;
  for(int s=0;s<2;++s){std::string name=s?"mc":"data";
   plots[s]=std::make_unique<TH2D>((name+"_energy_vs_tail_channels").c_str(),(name+" selected;Distinct low-tail channels;Positive good event energy [GeV]").c_str(),nx,-.5,xmax,100,fullLo,fullHi);
   plots[s+2]=std::make_unique<TH2D>((name+"_negative_lg_vs_tail_channels").c_str(),(name+" selected;Distinct low-tail channels;Negative-energy LG hits").c_str(),nx,-.5,xmax,ny,-.5,ymax);
   for(const auto&e:selectedEvents[s]){plots[s]->Fill(e.lowChannels,e.energy);plots[s+2]->Fill(e.lowChannels,e.negativeLG);}
  }
  for(int row=0;row<2;++row){double maximum=std::max(plots[row*2]->GetMaximum(),plots[row*2+1]->GetMaximum());for(int s=0;s<2;++s){auto& h=plots[row*2+s];correlations.cd(row*2+s+1);gPad->SetRightMargin(.16);gPad->SetLogz();h->SetMinimum(.5);h->SetMaximum(std::max(1.,maximum));h->GetZaxis()->SetTitle("Events");h->Write();h->Draw("COLZ");}}
  correlations.SaveAs((out/"figures/event_tail_correlations.png").c_str());
 }
 events.SetDirectory(nullptr);tailHits.SetDirectory(nullptr);result.Close();cache.SetDirectory(nullptr);cacheFile.Close();fs::rename(out/"tail_study.root.tmp",out/"tail_study.root");std::cout<<"Finished "<<out<<"; common channels="<<common<<"; hit/event counts agree"<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
