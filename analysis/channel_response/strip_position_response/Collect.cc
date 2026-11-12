#include "Common.hh"
#include <EBUdecode.h>
namespace {
using namespace spr;
using ID=std::tuple<int,int,int>;
struct Raw {
 TTree*t;int run=0,time=0,cycle=0,event=0;
 std::vector<int>*id=nullptr,*hit=nullptr,*gain=nullptr;std::vector<double>*hg=nullptr,*lg=nullptr;std::vector<std::vector<double>>*temp=nullptr;
 Raw(TTree*p):t(p){mmr::bind(t,"Run_Num",&run);mmr::bind(t,"Event_Time",&time);mmr::bind(t,"CycleID",&cycle);mmr::bind(t,"TriggerID",&event);}
 ID key()const{return {run,time,event};}
 void vectors(){mmr::bind(t,"CellID",&id);mmr::bind(t,"HitTag",&hit);mmr::bind(t,"GainTag",&gain);mmr::bind(t,"HG_Charge",&hg);mmr::bind(t,"LG_Charge",&lg);mmr::bind(t,"Temperature",&temp);t->SetCacheSize(32*1024*1024);t->AddBranchToCache("*",true);}
 void read(Long64_t i){mmr::read(t,i);require(id&&hit&&gain&&hg&&lg&&temp&&id->size()==hit->size()&&id->size()==gain->size()&&id->size()==hg->size()&&id->size()==lg->size(),"Raw vector size mismatch");}
 ~Raw(){t->ResetBranchAddresses();delete id;delete hit;delete gain;delete hg;delete lg;delete temp;}
};
struct Calib {
 TTree*t;int run=0,time=0,event=0;std::vector<int>*id=nullptr;std::vector<double>*energy=nullptr,*temp=nullptr;
 Calib(TTree*p):t(p){mmr::bind(t,"Run_Num",&run);mmr::bind(t,"Event_Time",&time);mmr::bind(t,"Event_Num",&event);}
 ID key()const{return {run,time,event};}
 void vectors(){mmr::bind(t,"CellID",&id);mmr::bind(t,"Hit_Energy",&energy);mmr::bind(t,"NewTemperature",&temp);}
 void read(Long64_t i){mmr::read(t,i);require(id&&energy&&temp&&id->size()==energy->size()&&id->size()==temp->size(),"Calibration vector size mismatch");}
 ~Calib(){t->ResetBranchAddresses();delete id;delete energy;delete temp;}
};
struct Event {
 int sample=0,run=0,event=0,time=0,cycle=0,paired=0,recorded=0,hit_records=0,physical_strips=0,duplicate_strips=0,calib_records=-1,calib_strips=-1,negative_energy=-1,zero_energy=-1,ssa_virtual_hits=-1,tail_cut_applied=0,energy_cut_applied=0;
 double energy_sum=qnan,positive_energy_sum=qnan;
 void branches(TTree&t){
#define B(x) t.Branch(#x,&x)
 B(sample);B(run);B(event);B(time);B(cycle);B(paired);B(recorded);B(hit_records);B(physical_strips);B(duplicate_strips);B(calib_records);B(calib_strips);B(negative_energy);B(zero_energy);B(ssa_virtual_hits);B(tail_cut_applied);B(energy_cut_applied);B(energy_sum);B(positive_energy_sum);
#undef B
 }
};
}
int main(int argc,char**argv)try{
 using namespace spr;using namespace mmr_track;
 require(argc==4,"Collect manifest.tsv configuration.tsv output-directory");mmr::Config cfg(argv[2]);Settings settings(cfg);auto selected=selection(cfg.text("channels","all"));std::set<int>mask;if(cfg.text("mask","none")!="none")mask=selection(cfg.text("mask"));
 std::array<std::map<int,mmr::Parameter>,2>pars;
 for(int s=0;s<2;++s){auto p=sample(s);pars[s]=mmr::loadParameters(cfg.text(p+"-pedestal"),cfg.text(p+"-mip"),cfg.text(p+"-hl"),cfg.text(p+"-threshold"));}
 std::string out=argv[3];require(!fs::exists(out+"/crossings.root"),"Collection already exists");
 auto params=table(out+"/parameters.tsv");params<<"sample\tcellid\tpedestal\tpedestal_present\tthreshold\tsigma\tthreshold_fallback\tcurrent_MIP\tMIP_fallback\tHL_gain\tHL_intercept\tx_mm\ty_mm\tz_mm\tmasked\n";
 for(int s=0;s<2;++s)for(auto&[id,p]:pars[s]){auto g=mmr::geometry(id);params<<sample(s)<<'\t'<<id<<'\t'<<p.ph<<'\t'<<p.pedPresent<<'\t'<<p.threshold<<'\t'<<p.sigma<<'\t'<<p.thresholdFallback<<'\t'<<p.mip<<'\t'<<p.mipFallback<<'\t'<<p.gain<<'\t'<<p.offset<<'\t'<<g[0]<<'\t'<<g[1]<<'\t'<<g[2]<<'\t'<<mask.count(id)<<'\n';}
 std::array<std::array<std::array<int,5>,42>,30>grid;for(auto&l:grid)for(auto&r:l)r.fill(-1);std::array<double,30>z{};std::map<int,StripRectangle>strips;
 for(int id:mmr::channels()){StripRectangle s{id,id/100000,mmr::geometry(id)};double fine=s.layer%2?s.p[0]:s.p[1],along=s.layer%2?s.p[1]:s.p[0];int row=std::lround(20.5-fine/5.3),col=std::lround(2-along/45.4);require(row>=0&&row<42&&col>=0&&col<5&&grid[s.layer][row][col]<0,"Invalid geometry grid");grid[s.layer][row][col]=id;strips.emplace(id,s);z[s.layer]=s.p[2];}
 std::ifstream manifest(argv[1]);require(bool(manifest),"Missing manifest");std::vector<std::vector<std::string>>files;std::set<std::pair<int,int>>seen;
 for(std::string line;std::getline(manifest,line);){if(line.empty()||line[0]=='#')continue;auto v=mmr::split(line);require(v.size()==6,"Expected 6 manifest columns");int s=std::stoi(v[0]),run=std::stoi(v[1]);require((s==0||s==1)&&seen.insert({s,run}).second,"Duplicate sample/run or invalid sample");files.push_back(v);}
 auto audit=table(out+"/event_audit.tsv");audit<<"sample\trun\traw_entries\tcalib_entries\tsampled\tunique\tpaired_by_ID\treordered\tmissing_or_ambiguous_calib_ID\tduplicate_raw_ID\tpredicted\tobserved\tusable\n";
 auto provenance=table(out+"/production_audit.tsv");provenance<<"sample\trun\ttruth_entries\tprimary_PDG\tprimary_energy_stored\ttruth_join\tproduction_link\tCT_state\n";
 TFile flat((out+"/crossings.root").c_str(),"CREATE");require(!flat.IsZombie(),"Cannot create crossings");Crossing x;Event e;TTree ct("Crossings","Leave entire target layer out; temperature from raw sensors, ADC independent of calibration survival"),et("Events","Recorded hits and physical strips; SSA not read (-1)");x.branches(ct);e.branches(et);
 for(auto&v:files){int s=std::stoi(v[0]),run=std::stoi(v[1]);auto rf=mmr::openFile(v[2]),cf=mmr::openFile(v[3]);Raw r(mmr::tree(*rf,"Raw_Hit"));Calib c(mmr::tree(*cf,"Calib_Hit"));
  std::map<ID,int> rawCounts;std::map<ID,std::vector<Long64_t>>calIndex;
  for(Long64_t i=0;i<r.t->GetEntries();++i){mmr::read(r.t,i);require(s||r.run==run,"Data run differs from manifest");++rawCounts[r.key()];}
  for(Long64_t i=0;i<c.t->GetEntries();++i){mmr::read(c.t,i);require(s||c.run==run,"Data calibration run differs from manifest");calIndex[c.key()].push_back(i);}
  if(s){require(v[5]!="-","Muon truth input is required for primary audit");auto tf=mmr::openFile(v[5]);auto*t=mmr::tree(*tf,"MC_Truth");int pdg=0;double en=0;mmr::bind(t,"PrimaryPDG",&pdg);mmr::bind(t,"PrimaryEnergy",&en);require(t->GetEntries()>0,"Empty truth");for(Long64_t j=0;j<std::min<Long64_t>(100,t->GetEntries());++j){mmr::read(t,(2*j+1)*t->GetEntries()/(2*std::min<Long64_t>(100,t->GetEntries())));require(std::abs(pdg)==13,"Electron/non-muon truth rejected");require(std::isfinite(en)&&en>0,"Invalid primary energy");}provenance<<"mc\t"<<run<<'\t'<<t->GetEntries()<<'\t'<<pdg<<'\t'<<en<<"\tnot_joined_no_saved_truth_event_ID\tmanifest_pair_not_generation_proof\tunverified_for_historical_input\n";t->ResetBranchAddresses();}
  else provenance<<"data\t"<<run<<"\t-1\t-1\tnan\tnot_applicable\tID_checked_decode_calib\tnot_applied_by_this_analysis\n";
  r.vectors();c.vectors();Long64_t n=r.t->GetEntries(),take=settings.maxEvents?std::min(n,settings.maxEvents):n,unique=0,paired=0,reordered=0,unpaired=0,dup=0,pred=0,observed=0,usable=0;
  for(Long64_t j=0;j<take;++j){Long64_t entry=(2*j+1)*n/(2*take);r.read(entry);if(rawCounts.at(r.key())!=1){++dup;continue;}++unique;
   auto ci=calIndex.find(r.key());bool has=ci!=calIndex.end()&&ci->second.size()==1;std::map<int,std::vector<size_t>>calCells;
   e=Event{};e.sample=s;e.run=run;e.event=r.event;e.time=r.time;e.cycle=r.cycle;e.recorded=r.id->size();e.paired=has;
   if(has){c.read(ci->second[0]);require(c.key()==r.key(),"Event ID join mismatch");++paired;reordered+=ci->second[0]!=entry;e.calib_records=c.id->size();e.negative_energy=e.zero_energy=0;e.energy_sum=e.positive_energy_sum=0;std::set<int>physical;
    for(size_t k=0;k<c.id->size();++k){int id=mmr::physical(c.id->at(k));if(!mmr::connected(id)||mask.count(id))continue;physical.insert(id);calCells[c.id->at(k)].push_back(k);double val=c.energy->at(k);e.negative_energy+=val<0;e.zero_energy+=val==0;e.energy_sum+=val;e.positive_energy_sum+=std::isfinite(val)?std::max(0.,val):qnan;}e.calib_strips=physical.size();
   }else ++unpaired;
   std::map<int,std::vector<size_t>>hits;std::array<std::vector<double>,30>fine;
   for(size_t k=0;k<r.id->size();++k)if(r.hit->at(k)==1){int id=mmr::physical(r.id->at(k));if(!mmr::connected(id)||mask.count(id))continue;hits[id].push_back(k);auto&p=strips.at(id);fine[p.layer].push_back(p.layer%2?p.p[0]:p.p[1]);++e.hit_records;}
   e.physical_strips=hits.size();for(auto&[id,h]:hits)e.duplicate_strips+=h.size()>1;et.Fill();
   for(auto&q:fine){std::sort(q.begin(),q.end());q.erase(std::unique(q.begin(),q.end()),q.end());}
   auto fullX=track(fine,z,1,-1,settings),fullY=track(fine,z,0,-1,settings);
   for(int l=0;l<30;++l){auto first=selected.lower_bound(l*100000);if(first==selected.end()||*first>=(l+1)*100000)continue;
    auto fx=l%2?track(fine,z,1,l,settings):fullX,fy=l%2?fullY:track(fine,z,0,l,settings);if(!good(fx,settings)||!good(fy,settings))continue;
    double px=fx.a+fx.b*z[l],py=fy.a+fy.b*z[l];int row=std::lround(20.5-(l%2?px:py)/5.3),col=std::lround(2-(l%2?py:px)/45.4);if(row<0||row>=42||col<0||col>=5)continue;int id=grid[l][row][col];if(!selected.count(id)||mask.count(id))continue;
    double len,u,w;auto&strip=strips.at(id);if(!chord(strip,px,py,fx.b,fy.b,settings,len,u,w))continue;
    x=Crossing{};x.sample=s;x.run=run;x.source_run=r.run;x.event=r.event;x.event_time=r.time;x.cycle=r.cycle;x.entry=entry;x.cellid=id;x.calib_paired=has;x.sx=fx.b;x.sy=fy.b;x.length=len;x.longitudinal=u;x.transverse=w;x.track_rms=std::sqrt((fx.rms*fx.rms*fx.n+fy.rms*fy.rms*fy.n)/(fx.n+fy.n));
    x.temperature_state=0;x.temperature=qnan;
    if(size_t(l)>=r.temp->size()||r.temp->at(l).empty())x.temperature_state=1;
    else if(r.temp->at(l).size()!=16||!std::all_of(r.temp->at(l).begin(),r.temp->at(l).end(),[](double v){return std::isfinite(v)&&v>-20&&v<60;}))x.temperature_state=2;
    else{auto pos=strip.p;x.temperature=tempReconstruction(l,pos.data(),r.temp->at(l));}
    double factor=mmr::tempFactor(id,x.temperature);if(x.temperature_state==0&&!(std::isfinite(factor)&&factor>0))x.temperature_state=2;
    auto h=hits.find(id);x.observed=h!=hits.end();x.multiplicity=x.observed?h->second.size():0;
    if(x.multiplicity==1){auto k=h->second[0];x.full_cellid=r.id->at(k);x.raw_hg=r.hg->at(k);x.raw_lg=r.lg->at(k);x.gain=r.gain->at(k);auto&p=pars[s].at(id);
     x.raw_adc=x.raw_hg-p.ph;x.path_adc=x.raw_adc*2/len;x.current_adc=x.path_adc*factor;
     x.usable=p.pedPresent&&std::isfinite(x.raw_hg)&&std::isfinite(x.raw_lg)&&x.raw_hg>=0&&x.raw_hg<2600&&x.temperature_state==0;
     auto cc=calCells.find(x.full_cellid);if(cc!=calCells.end()&&cc->second.size()==1){x.saved_temperature=c.temp->at(cc->second[0]);x.saved_energy=c.energy->at(cc->second[0]);}
    }
    ct.Fill();++pred;observed+=x.observed;usable+=x.usable;
   }
  }
  audit<<sample(s)<<'\t'<<run<<'\t'<<n<<'\t'<<c.t->GetEntries()<<'\t'<<take<<'\t'<<unique<<'\t'<<paired<<'\t'<<reordered<<'\t'<<unpaired<<'\t'<<dup<<'\t'<<pred<<'\t'<<observed<<'\t'<<usable<<'\n';std::cout<<sample(s)<<" run="<<run<<" sampled="<<take<<" paired="<<paired<<" crossings="<<pred<<" usable="<<usable<<std::endl;
 }
 flat.cd();ct.Write();et.Write();TParameter<int>("strip_position_schema",schema).Write();TNamed("selection","Muon only; no energy/tail/CoG selection; ADC response requires one HG<2600 hit and valid raw temperature; all three ADC stages use identical hits. Target layer excluded before track construction. No SSA.").Write();flat.Close();require(bool(audit)&&bool(params)&&bool(provenance),"Output table write failed");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
