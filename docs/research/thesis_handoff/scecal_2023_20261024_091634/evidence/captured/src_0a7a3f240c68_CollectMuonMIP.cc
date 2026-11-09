#include "MuonMIPCommon.hh"
#include "TrackGeometry.hh"
#include <TNamed.h>
#include <TParameter.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
namespace collect {
using mmr::require;
using Identity=std::pair<int,int>;
using namespace mmr_track;
struct Raw {
 TTree*t;int run=0,time=0,event=0;std::vector<int>*id=nullptr,*tag=nullptr;std::vector<double>*h=nullptr,*l=nullptr;
 Raw(TTree*p):t(p){mmr::bind(t,"Run_Num",&run);mmr::bind(t,"Event_Time",&time);mmr::bind(t,"TriggerID",&event);mmr::bind(t,"CellID",&id);mmr::bind(t,"HitTag",&tag);mmr::bind(t,"HG_Charge",&h);mmr::bind(t,"LG_Charge",&l);}
 void read(Long64_t i){mmr::read(t,i);require(id&&tag&&h&&l&&id->size()==tag->size()&&id->size()==h->size()&&id->size()==l->size(),"Raw vector lengths differ");}
 ~Raw(){t->ResetBranchAddresses();delete id;delete tag;delete h;delete l;}
};
struct Calib {
 TTree*t;int run=0,time=0,event=0;std::vector<int>*id=nullptr;std::vector<double>*e=nullptr,*h=nullptr,*temp=nullptr;
 Calib(TTree*p):t(p){mmr::bind(t,"Run_Num",&run);mmr::bind(t,"Event_Time",&time);mmr::bind(t,"Event_Num",&event);mmr::bind(t,"CellID",&id);mmr::bind(t,"Hit_Energy",&e);mmr::bind(t,"Hit_HG_Energy",&h);mmr::bind(t,"NewTemperature",&temp);}
 void read(Long64_t i){mmr::read(t,i);require(id&&e&&h&&temp&&id->size()==e->size()&&id->size()==h->size()&&id->size()==temp->size(),"Calib vector lengths differ");}
 ~Calib(){t->ResetBranchAddresses();delete id;delete e;delete h;delete temp;}
};
struct Truth {
 TTree*t;int pdg=0;double primary=NaN,x=NaN,y=NaN,z=NaN;std::vector<int>*id=nullptr;std::vector<double>*e=nullptr,*dir=nullptr;
 Truth(TTree*p):t(p){mmr::bind(t,"CellID",&id);mmr::bind(t,"Hit_Energy",&e);mmr::bind(t,"PrimaryPDG",&pdg);mmr::bind(t,"PrimaryEnergy",&primary);mmr::bind(t,"PrimaryPosX",&x);mmr::bind(t,"PrimaryPosY",&y);mmr::bind(t,"PrimaryPosZ",&z);mmr::bind(t,"PrimaryDirection",&dir);}
 void read(Long64_t i){mmr::read(t,i);require(id&&e&&id->size()==e->size(),"Truth vector lengths differ");require(std::abs(pdg)==13,"MC primary is not a muon");require(std::isfinite(primary)&&primary>0&&std::isfinite(x)&&std::isfinite(y)&&std::isfinite(z)&&dir&&dir->size()==3&&std::all_of(dir->begin(),dir->end(),[](double v){return std::isfinite(v);}),"Invalid muon primary");}
 ~Truth(){t->ResetBranchAddresses();delete id;delete e;delete dir;}
};
std::vector<Long64_t>sampled(Long64_t n,Long64_t maximum){Long64_t take=maximum>0?std::min(n,maximum):n;std::vector<Long64_t>v;v.reserve(take);for(Long64_t i=0;i<take;++i)v.push_back((2*i+1)*n/(2*take));return v;}
struct Counts{Long64_t raw=0,calib=0,truth=0,sampled=0,paired=0,missing=0,duplicateRaw=0,duplicateCalib=0,reordered=0,ambiguousRaw=0,ambiguousCalib=0,cross=0,observed=0,response=0,truthHits=0,invalidChord=0,truthMissingRaw=0,recovered=0,nonfinite=0;};
struct Crossing {
 int source_run=0,file_index=0,sample=0,run=0,event=0,event_time=0,cellid=0,observed=0,ambiguous=0,calib_ok=0,response_ok=0;
 Long64_t event_index=0;
 double raw_hg=NaN,raw_lg=NaN,saved_hg=NaN,saved_energy=NaN,temperature=NaN,sx=NaN,sy=NaN,length=NaN,longitudinal=NaN,transverse=NaN,track_rms=NaN,truth_edep=NaN;
 void branches(TTree&t){t.Branch("source_run",&source_run);t.Branch("file_index",&file_index);t.Branch("sample",&sample);t.Branch("run",&run);t.Branch("event",&event);t.Branch("event_time",&event_time);t.Branch("cellid",&cellid);t.Branch("observed",&observed);t.Branch("ambiguous",&ambiguous);t.Branch("calib_ok",&calib_ok);t.Branch("response_ok",&response_ok);t.Branch("event_index",&event_index);t.Branch("raw_hg",&raw_hg);t.Branch("raw_lg",&raw_lg);t.Branch("saved_hg",&saved_hg);t.Branch("saved_energy",&saved_energy);t.Branch("temperature",&temperature);t.Branch("sx",&sx);t.Branch("sy",&sy);t.Branch("length",&length);t.Branch("longitudinal",&longitudinal);t.Branch("transverse",&transverse);t.Branch("track_rms",&track_rms);t.Branch("truth_edep",&truth_edep);}
};
struct TruthHit {
 int source_run=0,file_index=0,sample=1,run=0,event=0,event_time=0,cellid=0,truth_detected=0;Long64_t event_index=0;double truth_edep=NaN,truth_length=NaN;
 void branches(TTree&t){t.Branch("source_run",&source_run);t.Branch("file_index",&file_index);t.Branch("sample",&sample);t.Branch("run",&run);t.Branch("event",&event);t.Branch("event_time",&event_time);t.Branch("cellid",&cellid);t.Branch("event_index",&event_index);t.Branch("truth_edep",&truth_edep);t.Branch("truth_length",&truth_length);t.Branch("truth_detected",&truth_detected);}
};
struct EventSummary {
 int source_run=0,file_index=0,sample=0,run=0,event=0,event_time=0,calib_paired=0,event_hits=0,selected_crossings=0,observed_crossings=0,response_crossings=0;Long64_t event_index=0;double event_energy=NaN,event_positive_energy=NaN;
 void branches(TTree&t){t.Branch("source_run",&source_run);t.Branch("file_index",&file_index);t.Branch("sample",&sample);t.Branch("run",&run);t.Branch("event",&event);t.Branch("event_time",&event_time);t.Branch("event_index",&event_index);t.Branch("calib_paired",&calib_paired);t.Branch("event_hits",&event_hits);t.Branch("selected_crossings",&selected_crossings);t.Branch("observed_crossings",&observed_crossings);t.Branch("response_crossings",&response_crossings);t.Branch("event_energy",&event_energy);t.Branch("event_positive_energy",&event_positive_energy);}
};
// Build identities once using scalars only, then restore the vector branches.
std::map<Identity,std::vector<Long64_t>>rawIndex(Raw&r){auto*t=r.t;t->SetBranchStatus("*",false);for(auto*n:{"Run_Num","Event_Time","TriggerID"})t->SetBranchStatus(n,true);std::map<Identity,std::vector<Long64_t>>idx;for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);idx[{r.time,r.event}].push_back(i);}for(auto*n:{"CellID","HitTag","HG_Charge","LG_Charge"})t->SetBranchStatus(n,true);return idx;}
std::map<Identity,std::vector<Long64_t>>calibIndex(Calib&r){auto*t=r.t;t->SetBranchStatus("*",false);for(auto*n:{"Run_Num","Event_Time","Event_Num"})t->SetBranchStatus(n,true);std::map<Identity,std::vector<Long64_t>>idx;for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);idx[{r.time,r.event}].push_back(i);}for(auto*n:{"CellID","Hit_Energy","Hit_HG_Energy","NewTemperature"})t->SetBranchStatus(n,true);return idx;}
std::map<int,double>deposits(Truth&t){std::map<int,double>v;for(size_t i=0;i<t.id->size();++i){int id=mmr::physical(t.id->at(i));if(!mmr::connected(id))continue;double e=t.e->at(i);require(std::isfinite(e)&&e>=0,"Invalid truth deposit");v[id]+=e;}return v;}
void validation(std::ostream&o,int file,int sample,int run,const Counts&c){o<<file<<'\t'<<sample<<'\t'<<run;for(auto v:{c.raw,c.calib,c.truth,c.sampled,c.paired,c.missing,c.duplicateRaw,c.duplicateCalib,c.reordered,c.ambiguousRaw,c.ambiguousCalib,c.cross,c.observed,c.response,c.truthHits,c.invalidChord,c.truthMissingRaw,c.recovered,c.nonfinite})o<<'\t'<<v;o<<'\n';}
} // namespace collect
int main(int argc,char**argv)try{
 using namespace collect;
 require(argc==5,"Usage: CollectMuonMIP manifest.tsv config.tsv flat.root validation.tsv");
 mmr::Config config(argv[2]);Settings settings(config);
 require(!mmr::fs::exists(argv[3])&&!mmr::fs::exists(argv[4]),"Output ROOT or validation table already exists");
 std::ifstream manifest(argv[1]);require(bool(manifest),"Cannot read input manifest");
 std::ofstream audit(argv[4]);require(bool(audit),"Cannot create validation table");
 audit<<"file\tsample\trun\traw_entries\tcalib_entries\ttruth_entries\tsampled_raw\tpaired\tmissing_calib\tduplicate_raw_events\tduplicate_calib_events\treordered_pairs\tambiguous_raw_cells\tambiguous_calib_cells\tcrossings\tobserved\tresponse\ttruth_hits\ttruth_invalid_primary_chord\ttruth_without_decode_identity\trecovered_files\tnonfinite_saved\n";
 TFile output(argv[3],"CREATE");require(!output.IsZombie(),"Cannot create output ROOT");
 Crossing crossing;TruthHit truthHit;EventSummary summary;
 TTree crossings("Crossings","Target-layer-independent raw crossings; conditional efficiency"),truthHits("TruthHits","Generated deposits independent of tracking; primary-line chord proxy"),events("Events","Actual saved connected hit sums, one count per raw event");
 crossing.branches(crossings);truthHit.branches(truthHits);summary.branches(events);
 std::array<std::array<std::array<int,5>,42>,30>stripGrid;for(auto&layer:stripGrid)for(auto&row:layer)row.fill(-1);std::array<double,30>layerZ{};std::map<int,StripRectangle>stripById;
 for(int id:mmr::channels()){StripRectangle s{id,id/100000,mmr::geometry(id)};double fine=s.layer%2?s.p[0]:s.p[1],longAt=s.layer%2?s.p[1]:s.p[0];int row=std::lround(20.5-fine/5.3),column=std::lround(2-longAt/45.4);require(row>=0&&row<42&&column>=0&&column<5&&stripGrid[s.layer][row][column]<0,"Invalid geometry grid");require(std::abs(fine-(20.5-row)*5.3)<1e-6&&std::abs(longAt-(2-column)*45.4)<1e-6,"Geometry does not match regular strip grid");stripGrid[s.layer][row][column]=id;stripById.emplace(id,s);layerZ[s.layer]=s.p[2];}
 int fileIndex=-1;
 for(std::string line;std::getline(manifest,line);){
  if(line.empty()||line[0]=='#')continue;
  auto fields=mmr::split(line);if(fields.size()==6&&fields[0]=="sample")continue;
  require(fields.size()==6,"Manifest must contain sample/run/raw/calib/track/truth");
  ++fileIndex;int sample=std::stoi(fields[0]),manifestRun=std::stoi(fields[1]);require(sample==0||sample==1,"Sample must be 0(data) or 1(MC)");
  auto rawFile=mmr::openFile(fields[2],settings.allow),calibFile=mmr::openFile(fields[3],settings.allow);
  Raw raw(mmr::tree(*rawFile,"Raw_Hit"));Calib calib(mmr::tree(*calibFile,"Calib_Hit"));
  mmr::File truthFile;std::unique_ptr<Truth>truth;
  if(sample==1){require(fields[5]!="-","MC generated truth is required");truthFile=mmr::openFile(fields[5],settings.allow);truth=std::make_unique<Truth>(mmr::tree(*truthFile,"MC_Truth"));}
  Counts count;count.raw=raw.t->GetEntries();count.calib=calib.t->GetEntries();count.truth=truth?truth->t->GetEntries():0;
  count.recovered=rawFile->TestBit(TFile::kRecovered)+calibFile->TestBit(TFile::kRecovered)+(truthFile?truthFile->TestBit(TFile::kRecovered):0);
  auto rawIds=rawIndex(raw),calibIds=calibIndex(calib);
  for(auto&v:rawIds)if(v.second.size()>1)count.duplicateRaw+=v.second.size();
  for(auto&v:calibIds)if(v.second.size()>1)count.duplicateCalib+=v.second.size();
  for(Long64_t entry:sampled(count.raw,settings.maxEvents)){
   raw.read(entry);++count.sampled;require(sample==1||raw.run==manifestRun,"Data raw run differs from manifest");
   Identity key{raw.time,raw.event};bool rawUnique=rawIds.at(key).size()==1;auto paired=calibIds.find(key);bool calibPaired=rawUnique&&paired!=calibIds.end()&&paired->second.size()==1;
   if(calibPaired){Long64_t ci=paired->second.front();calib.read(ci);require(calib.run==raw.run&&calib.time==raw.time&&calib.event==raw.event,"Raw/calib identity mismatch");++count.paired;count.reordered+=ci!=entry;}else ++count.missing;
   // Duplicate event identities cannot establish a response denominator.
   if(!rawUnique)continue;
   std::map<int,std::vector<size_t>>rawCells,calibFull;std::array<std::vector<double>,30>fine;
   for(size_t j=0;j<raw.id->size();++j)if(raw.tag->at(j)==1){int id=mmr::physical(raw.id->at(j));if(!mmr::connected(id))continue;require(std::isfinite(raw.h->at(j))&&std::isfinite(raw.l->at(j)),"Invalid raw charge");rawCells[id].push_back(j);auto p=stripById.at(id).p;fine[id/100000].push_back(id/100000%2?p[0]:p[1]);}
   for(auto&q:fine){std::sort(q.begin(),q.end());q.erase(std::unique(q.begin(),q.end()),q.end());}
   for(auto&v:rawCells)count.ambiguousRaw+=v.second.size()>1;
   summary=EventSummary{};summary.file_index=fileIndex;summary.sample=sample;summary.run=manifestRun;summary.source_run=raw.run;summary.event=raw.event;summary.event_time=raw.time;summary.event_index=entry;summary.calib_paired=calibPaired;
   if(calibPaired){summary.event_energy=summary.event_positive_energy=0;for(size_t j=0;j<calib.id->size();++j){int id=mmr::physical(calib.id->at(j));if(!mmr::connected(id))continue;calibFull[calib.id->at(j)].push_back(j);++summary.event_hits;double e=calib.e->at(j);if(!std::isfinite(e)){++count.nonfinite;summary.event_energy=summary.event_positive_energy=NaN;continue;}summary.event_energy+=e;summary.event_positive_energy+=std::max(0.,e);}for(auto&v:calibFull)count.ambiguousCalib+=v.second.size()>1;}
   std::map<int,double>truthDeposits;
   if(truth){require(raw.time==raw.event&&raw.event>=0&&raw.event<count.truth,"MC truth/decode trigger mapping is not established");truth->read(raw.event);truthDeposits=deposits(*truth);}
   const Fit fullX=track(fine,layerZ,1,-1,settings),fullY=track(fine,layerZ,0,-1,settings);
   for(int target=0;target<30;++target){
    // The other projection contains no target-layer coordinates by construction.
    const Fit fx=target%2?track(fine,layerZ,1,target,settings):fullX;
    const Fit fy=target%2?fullY:track(fine,layerZ,0,target,settings);
    if(!good(fx,settings)||!good(fy,settings))continue;
    double px=fx.a+fx.b*layerZ[target],py=fy.a+fy.b*layerZ[target];
    int row=std::lround(20.5-(target%2?px:py)/5.3),column=std::lround(2-(target%2?py:px)/45.4);
    if(row<0||row>=42||column<0||column>=5)continue;
    int targetId=stripGrid[target][row][column];if(targetId<0)continue;
    const auto&strip=stripById.at(targetId);
    {
     double length,longAt,transAt;if(!chord(strip,px,py,fx.b,fy.b,settings,length,longAt,transAt))continue;
     crossing=Crossing{};crossing.file_index=fileIndex;crossing.sample=sample;crossing.run=manifestRun;crossing.source_run=raw.run;crossing.event=raw.event;crossing.event_time=raw.time;crossing.event_index=entry;crossing.cellid=strip.id;crossing.sx=fx.b;crossing.sy=fy.b;crossing.length=length;crossing.longitudinal=longAt;crossing.transverse=transAt;crossing.track_rms=std::sqrt((fx.rms*fx.rms*fx.n+fy.rms*fy.rms*fy.n)/(fx.n+fy.n));
     crossing.truth_edep=truth?(truthDeposits.count(strip.id)?truthDeposits.at(strip.id):0.):NaN;
     auto hit=rawCells.find(strip.id);crossing.observed=hit!=rawCells.end();crossing.ambiguous=crossing.observed&&hit->second.size()!=1;
     if(crossing.observed&&!crossing.ambiguous){
      size_t j=hit->second.front();crossing.raw_hg=raw.h->at(j);crossing.raw_lg=raw.l->at(j);auto saved=calibFull.find(raw.id->at(j));
      if(calibPaired&&saved!=calibFull.end()&&saved->second.size()==1){size_t k=saved->second.front();crossing.saved_hg=calib.h->at(k);crossing.saved_energy=calib.e->at(k);crossing.temperature=calib.temp->at(k);crossing.calib_ok=std::isfinite(crossing.saved_hg)&&std::isfinite(crossing.saved_energy);double C=mmr::tempFactor(strip.id,crossing.temperature);crossing.response_ok=crossing.calib_ok&&std::isfinite(crossing.temperature)&&std::isfinite(C)&&C>0&&crossing.raw_hg<2600&&crossing.raw_hg!=4000;}
      else if(saved!=calibFull.end()&&saved->second.size()>1)crossing.ambiguous=1;
     }
     crossings.Fill();++count.cross;++summary.selected_crossings;count.observed+=crossing.observed;summary.observed_crossings+=crossing.observed;count.response+=crossing.response_ok;summary.response_crossings+=crossing.response_ok;
    }
   }
   events.Fill();
  }
  // Independent generated-event sample: no reconstructed track/target required.
  if(truth)for(Long64_t entry:sampled(count.truth,settings.maxEvents)){
   truth->read(entry);auto dep=deposits(*truth);auto decoded=rawIds.find({int(entry),int(entry)});std::set<int>detected;
   if(decoded!=rawIds.end()&&decoded->second.size()==1){raw.read(decoded->second.front());require(raw.time==entry&&raw.event==entry,"Truth/decode identity mismatch");for(size_t j=0;j<raw.id->size();++j)if(raw.tag->at(j)==1){int id=mmr::physical(raw.id->at(j));if(mmr::connected(id))detected.insert(id);}}else ++count.truthMissingRaw;
   double norm=std::sqrt(std::inner_product(truth->dir->begin(),truth->dir->end(),truth->dir->begin(),0.));bool direction=std::isfinite(norm)&&norm>0&&std::abs(truth->dir->at(2)/norm)>1e-9;
   double sx=direction?truth->dir->at(0)/truth->dir->at(2):NaN,sy=direction?truth->dir->at(1)/truth->dir->at(2):NaN;
   for(auto&v:dep){truthHit=TruthHit{};truthHit.file_index=fileIndex;truthHit.run=manifestRun;truthHit.event=entry;truthHit.event_time=entry;truthHit.event_index=entry;truthHit.cellid=v.first;truthHit.truth_edep=v.second;truthHit.truth_detected=detected.count(v.first);auto&strip=stripById.at(v.first);double length=NaN,longAt=NaN,transAt=NaN;double px=direction?truth->x+(strip.p[2]-truth->z)*sx:NaN,py=direction?truth->y+(strip.p[2]-truth->z)*sy:NaN;if(direction&&chord(strip,px,py,sx,sy,settings,length,longAt,transAt))truthHit.truth_length=length;else ++count.invalidChord;truthHits.Fill();++count.truthHits;}
  }
  validation(audit,fileIndex,sample,manifestRun,count);
  std::cout<<"file "<<fileIndex<<" sample="<<sample<<" run="<<manifestRun<<" sampled="<<count.sampled<<" crossings="<<count.cross<<" response="<<count.response<<" truth_hits="<<count.truthHits<<'\n';
 }
 output.cd();crossings.Write();truthHits.Write();events.Write();
 TNamed selection("selection_definition","Raw HitTag=1; unweighted fine-coordinate OLS/consensus; entire target layer excluded BEFORE track building; min-layers per projection; full active45x5x2mm chord with margins; denominator conditional on other-layer tracking, not unconditional hardware efficiency");selection.Write();
 TNamed definition("truth_path_definition","Primary straight-line chord with fiducial margins; NaN for deposits not crossed by that primary line. This is not a Geant4 accumulated step path.");definition.Write();
 TParameter<int>("allow_recovered",settings.allow).Write();output.Write();output.Close();require(bool(audit),"Validation write failed");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
