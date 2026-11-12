#include "Common.hh"
#include "../../resolution/HLTailSelection.hh"
#include <TKey.h>
double tempReconstruction(int,double*,std::vector<double>);
using namespace ts;
using Identity=std::tuple<int,int,int>;
struct Reader {
 TTree*t;int run=0,time=0,event=0;std::vector<int>*id=nullptr,*tag=nullptr;
 std::vector<double>*hg=nullptr,*lg=nullptr,*e=nullptr,*temp=nullptr;
 std::vector<std::vector<double>>*sensors=nullptr;bool cal;
 Reader(TTree*p,bool c):t(p),cal(c){bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,cal?"Event_Num":"TriggerID",&event);}
 Identity key(){return {run,time,event};}
 void vectors(){bind(t,"CellID",&id);if(cal){bind(t,"Hit_Energy",&e);bind(t,"NewTemperature",&temp);}else{bind(t,"HitTag",&tag);bind(t,"HG_Charge",&hg);bind(t,"LG_Charge",&lg);if(t->GetBranch("Temperature"))bind(t,"Temperature",&sensors);}t->SetCacheSize(16*1024*1024);t->AddBranchToCache("*",true);}
 void get(Long64_t i){mmr::read(t,i);if(id){if(cal)require(e&&temp&&e->size()==id->size()&&temp->size()==id->size(),"Calib vector mismatch");else require(tag&&hg&&lg&&tag->size()==id->size()&&hg->size()==id->size()&&lg->size()==id->size(),"Raw vector mismatch");}}
 double sensor(int cell){int l=cell/100000;if(!sensors||sensors->size()<=size_t(l)||sensors->at(l).size()!=16)return NaN;auto &v=sensors->at(l);for(double t:v)if(!std::isfinite(t)||t<=0||t>=60)return NaN;auto p=geometry(cell);return tempReconstruction(l,p.data(),v);}
 ~Reader(){t->ResetBranchAddresses();delete id;delete tag;delete hg;delete lg;delete e;delete temp;delete sensors;}
};
int main(int argc,char**argv)try{
 require(argc==4,"Collect manifest config output_dir");Config c(argv[2]);auto in=inputs(argv[1]);fs::path out=argv[3];require(!fs::exists(out/"flat.root"),"Existing flat ROOT");auto excluded=mask(c);mmr_track::Settings setting(c);
 std::map<int,Parameter> pars[2];for(int s=0;s<2;++s)pars[s]=loadParameters(c.text(s?"mc-pedestal":"data-pedestal"),c.text(s?"mc-mip":"data-mip"),c.text(s?"mc-hl":"data-hl"),c.text("threshold"));
 std::ofstream parameters(out/"calibration.tsv");parameters<<"sample\tcell\tpedHG\tpedLG\tMIP\tHL_gain\tHL_offset\tthreshold\tsigma\tmip_fallback\thl_fallback\tped_present\tthreshold_fallback\n";
 for(int s=0;s<2;++s)for(auto &[id,p]:pars[s])parameters<<s<<'\t'<<id<<'\t'<<p.ph<<'\t'<<p.pl<<'\t'<<p.mip<<'\t'<<p.gain<<'\t'<<p.offset<<'\t'<<p.threshold<<'\t'<<p.sigma<<'\t'<<p.mipFallback<<'\t'<<p.hlFallback<<'\t'<<p.pedPresent<<'\t'<<p.thresholdFallback<<'\n';
 std::unique_ptr<hl_tail::Selection> tails;for(auto&i:in)if(i.particle==1&&c.integer(i.sample?"mc-tail":"data-tail")) {tails=std::make_unique<hl_tail::Selection>(c.text("tail-results"));break;}
 TFile output((out/"flat.root").c_str(),"CREATE","",1);require(!output.IsZombie(),"Cannot create flat.root");Event ev;Hit h,ped;
 TTree events("Events","One sampled raw event including rejected events"),hits("Crossings","Target layer excluded before tracking; conditional denominator"),peds("Pedestal","HitTag=0 conditional population, memory cells separate");ev.io(events,false);h.io(hits,false);ped.io(peds,false);
 std::ofstream boundaries(out/"run_boundaries.tsv");boundaries<<"sample\tparticle\tprevious_run\trun\tprevious_file\tfile\tprevious_last_tick\tfirst_tick\tdifference_ticks\tstatus\n";std::map<std::pair<int,int>,std::tuple<int,int,int>>lastFile;
 std::ofstream audit(out/"input_audit.tsv"),counts(out/"counters.tsv"),provenance(out/"production_audit.tsv"),metadata(out/"root_metadata.tsv");
 audit<<"file\tsample\tparticle\trun\traw_entries\tcalib_entries\tfirst_tick\tlast_tick\tminimum_positive_delta\tdelta_gcd\trepeated_ticks\tbackward_steps\tsegments\tduplicate_raw_ids\tduplicate_calib_ids\tclock_status\tseconds_per_tick\tclock_evidence\n";
 counts<<"file\tsampled\tpaired\tmissing_or_ambiguous_id\treordered\ttail_checked\ttail_rejected\tcog_rejected\tselected\tcrossings\tobserved\tresponse\tpedestal_rows\n";
 provenance<<"file\ttruth_entries\tchecked_primaries\tprimary_pdg\tprimary_energy_MeV\ttruth_mapping\tproduction_status\tevidence\tCT_status\n";metadata<<"file\tinput\tkey\ttitle\n";
 std::array<double,30>z{};std::array<std::array<std::array<int,5>,42>,30>grid;for(auto&l:grid)for(auto&r:l)r.fill(-1);
 for(int id:channels()){auto p=geometry(id);int l=id/100000;grid[l][std::lround(20.5-(l%2?p[0]:p[1])/5.3)][std::lround(2-(l%2?p[1]:p[0])/45.4)]=id;z[l]=p[2];}
 int fi=-1;for(auto&i:in){++fi;std::cout<<"Collect file="<<fi<<" "<<particleName(i.particle)<<" "<<sampleName(i.sample)<<" run="<<i.run<<std::endl;
 auto rf=openFile(i.raw);Reader raw(tree(*rf,"Raw_Hit"),false);File cf;std::unique_ptr<Reader>cal;
 if(i.calib!="-"){cf=openFile(i.calib);cal=std::make_unique<Reader>(tree(*cf,"Calib_Hit"),true);}
 for(auto f:{rf.get(),cf.get()})if(f){TIter next(f->GetListOfKeys());while(auto*k=dynamic_cast<TKey*>(next())){if(std::string(k->GetClassName())=="TNamed"||std::string(k->GetClassName()).find("TParameter")==0){std::unique_ptr<TObject>obj(k->ReadObj());metadata<<fi<<'\t'<<f->GetName()<<'\t'<<k->GetName()<<'\t'<<obj->GetTitle()<<'\n';}}}
 std::map<Identity,Long64_t> ri,ci;std::vector<int>segments;std::vector<int>starts;
 int prev=0,first=0,last=0,seg=0;long long repeats=0,back=0,minDelta=LLONG_MAX,gcdDelta=0,dr=0,dc=0;
 for(Long64_t e=0;e<raw.t->GetEntries();++e){raw.get(e);require(i.sample||raw.run==i.run,"Data run differs from manifest");if(e==0){first=raw.time;starts.push_back(raw.time);}else{long long d=static_cast<long long>(raw.time)-prev;if(d==0)++repeats;else if(d<0){++back;++seg;starts.push_back(raw.time);}else{minDelta=std::min(minDelta,d);gcdDelta=std::gcd(gcdDelta,d);}}prev=last=raw.time;segments.push_back(seg);auto it=ri.emplace(raw.key(),e);if(!it.second){it.first->second=-1;++dr;}}
 if(cal)for(Long64_t e=0;e<cal->t->GetEntries();++e){cal->get(e);auto it=ci.emplace(cal->key(),e);if(!it.second){it.first->second=-1;++dc;}}
 double tick=i.sample?0:c.number("seconds-per-tick",0);require(tick==0||c.text("clock-evidence","-")!="-","Physical time requires --clock-evidence");
 audit<<fi<<'\t'<<i.sample<<'\t'<<particleName(i.particle)<<'\t'<<i.run<<'\t'<<raw.t->GetEntries()<<'\t'<<(cal?cal->t->GetEntries():0)<<'\t'<<first<<'\t'<<last<<'\t'<<(minDelta==LLONG_MAX?-1:minDelta)<<'\t'<<gcdDelta<<'\t'<<repeats<<'\t'<<back<<'\t'<<seg+1<<'\t'<<dr<<'\t'<<dc<<'\t'<<(i.sample?"generated_event_index":tick?"user_documented_unit_epoch_unconfirmed":"integer_ticks_unit_unconfirmed")<<'\t'<<tick<<'\t'<<c.text("clock-evidence","-")<<'\n';
 auto previous=lastFile.find({i.sample,i.particle});if(previous!=lastFile.end()){auto[pr,pf,pt]=previous->second;boundaries<<i.sample<<'\t'<<particleName(i.particle)<<'\t'<<pr<<'\t'<<i.run<<'\t'<<pf<<'\t'<<fi<<'\t'<<pt<<'\t'<<first<<'\t'<<static_cast<long long>(first)-pt<<'\t'<<(i.sample?"generated_index_restarts_per_file":"run_boundaries_kept_separate_epoch_unconfirmed")<<'\n';}lastFile[{i.sample,i.particle}]={i.run,fi,last};
 std::vector<unsigned char>tailMask;if(i.particle==1&&c.integer(i.sample?"mc-tail":"data-tail"))tailMask=tails->mask(i.calib);
 File truth;TTree*tt=nullptr;int pdg=0;double primary=NaN;Long64_t truthN=0,truthChecked=0;
 if(i.sample){truth=openFile(i.truth);tt=tree(*truth,"MC_Truth");bind(tt,"PrimaryPDG",&pdg);bind(tt,"PrimaryEnergy",&primary);truthN=tt->GetEntries();}
 raw.vectors();if(cal)cal->vectors();
 Long64_t sampled=0,paired=0,missing=0,reorder=0,tailChecked=0,tailRejected=0,cogRejected=0,selected=0,crossings=0,observed=0,response=0,pedRows=0;
 Long64_t start=c.integer("start-entry",0),end=raw.t->GetEntries();if(setting.maxEvents)end=std::min(end,start+setting.maxEvents);
 for(Long64_t e=start;e<end;++e){raw.get(e);++sampled;ev=Event{};ev.sample=i.sample;ev.particle=i.particle;ev.file=fi;ev.run=i.run;ev.source_run=raw.run;ev.event=raw.event;ev.entry=e;ev.time=raw.time;ev.segment=segments[e];ev.seconds=tick?(static_cast<double>(raw.time)-starts[ev.segment])*tick:NaN;
 if(ri.at(raw.key())<0){++missing;events.Fill();continue;}
 auto found=ci.find(raw.key());Long64_t ce=-1;if(cal&&found!=ci.end()&&found->second>=0){ce=found->second;cal->get(ce);require(raw.key()==cal->key(),"Saved event identity mismatch");ev.paired=1;++paired;reorder+=ce!=e;}else if(cal)++missing;
 if(tt){require(raw.time==raw.event&&raw.event>=0&&raw.event<truthN,"MC generated index mapping not verified");mmr::read(tt,raw.event);require(std::abs(pdg)==(i.particle==0?13:11),"Wrong MC species");double target=c.number(i.particle==0?"muon-energy":"energy")*1000;require(std::abs(primary-target)<std::max(1.,target*.005),"Wrong primary energy");++truthChecked;}
 std::map<int,std::vector<size_t>> rawCells,calFull;std::array<std::vector<double>,30>fine;
 for(size_t j=0;j<raw.id->size();++j){int full=raw.id->at(j),id=physical(full);if(!connected(id)||excluded.count(id))continue;
 if(raw.tag->at(j)==1){++ev.raw_hits;rawCells[id].push_back(j);auto p=geometry(id);fine[id/100000].push_back(id/100000%2?p[0]:p[1]);}
 if(raw.tag->at(j)==0){ped=Hit{};static_cast<Base&>(ped)=ev;ped.cell=id;ped.memory=full%10000/100;ped.temperature=i.sample?NaN:raw.sensor(id);ped.sensor=ped.temperature;ped.temp_measured=std::isfinite(ped.temperature);auto&p=pars[i.sample].at(id);ped.before=p.pedPresent?raw.hg->at(j)-p.ph:NaN;ped.after=p.pedPresent?raw.lg->at(j)-p.pl:NaN;ped.valid=p.pedPresent&&std::isfinite(ped.before)&&std::isfinite(ped.after);peds.Fill();++pedRows;}
 }
 ev.raw_strips=rawCells.size();for(auto&q:fine){std::sort(q.begin(),q.end());q.erase(std::unique(q.begin(),q.end()),q.end());}
 if(ev.paired){std::set<int>phys;Moment temperatures;double cx[30]={},cy[30]={};ev.energy=ev.signed_energy=0;bool finite=true;
 for(size_t j=0;j<cal->id->size();++j){int full=cal->id->at(j),id=physical(full);if(!connected(id)||excluded.count(id))continue;calFull[full].push_back(j);++ev.recorded;phys.insert(id);double v=cal->e->at(j);if(!std::isfinite(v)){finite=false;continue;}ev.signed_energy+=v;ev.negative+=v<0;ev.zero+=v==0;if(v<0)continue;auto p=geometry(id);ev.energy+=v;ev.layer[id/100000]+=v;cx[id/100000]+=v*p[0];cy[id/100000]+=v*p[1];if(!i.sample&&std::isfinite(raw.sensor(id)))temperatures.add(cal->temp->at(j));}
 ev.strips=phys.size();ev.duplicate=ev.recorded-ev.strips;ev.temperature=temperatures.mean();if(!finite)ev.energy=ev.signed_energy=NaN;
 // e_ssa is a production tag: Calib_Hit holds strips, no SSA virtual hits are inferred.
 ev.cog_ok=1;double range=c.number("cog",200);if(range<200)for(int l:{9,10}){if(ev.layer[l]<=0||std::abs(cx[l]/ev.layer[l])>range||std::abs(cy[l]/ev.layer[l])>range)ev.cog_ok=0;}
 if(ev.layer[10]>0){ev.cogx=cx[10]/ev.layer[10];ev.cogy=cy[10]/ev.layer[10];}
 if(!tailMask.empty()){++tailChecked;ev.tail=tailMask.at(ce);tailRejected+=ev.tail;}
 if(i.particle==1&&!ev.tail&&!ev.cog_ok)++cogRejected;
 }
 ev.selected=i.particle==2||(i.particle==0?true:(ev.paired&&!ev.tail&&ev.cog_ok&&std::isfinite(ev.energy)));selected+=ev.selected;
 if(i.particle==0){using namespace mmr_track;const auto fx0=track(fine,z,1,-1,setting),fy0=track(fine,z,0,-1,setting);
 for(int l=0;l<30;++l){auto fx=l%2?track(fine,z,1,l,setting):fx0,fy=l%2?fy0:track(fine,z,0,l,setting);if(!good(fx,setting)||!good(fy,setting))continue;
 double x=fx.a+fx.b*z[l],y=fy.a+fy.b*z[l];int row=std::lround(20.5-(l%2?x:y)/5.3),col=std::lround(2-(l%2?y:x)/45.4);if(row<0||row>=42||col<0||col>=5)continue;int id=grid[l][row][col];if(id<0||excluded.count(id))continue;
 h=Hit{};static_cast<Base&>(h)=ev;h.cell=id;h.x=x;h.y=y;h.sx=fx.b;h.sy=fy.b;StripRectangle strip{id,l,geometry(id)};if(!chord(strip,x,y,fx.b,fy.b,setting,h.length,h.longitudinal,h.transverse))continue;
 h.sensor=raw.sensor(id);h.temp_measured=!i.sample&&std::isfinite(h.sensor);h.temperature=h.temp_measured?h.sensor:ts::NaN;
 auto r=rawCells.find(id);h.observed=r!=rawCells.end();h.unique=h.observed&&r->second.size()==1;
 if(h.unique){auto j=r->second.front();h.memory=raw.id->at(j)%10000/100;auto k=calFull.find(raw.id->at(j));auto&p=pars[i.sample].at(id);
 if(ev.paired&&k!=calFull.end()&&k->second.size()==1){double savedT=cal->temp->at(k->second.front());h.factor=tempFactor(id,savedT);if(h.temp_measured)h.temperature=savedT;
 h.before=(raw.hg->at(j)-p.ph)*2/h.length;h.after=h.before*h.factor;
 double path=2/h.length;h.low=std::max({10.,(p.threshold+3*p.sigma)*path,(p.threshold+3*p.sigma)*path*h.factor});h.high=std::min((c.number("hg-max",2200)-p.ph)*path,(c.number("hg-max",2200)-p.ph)*path*h.factor);
 if(!std::isfinite(p.threshold)||!std::isfinite(p.sigma)||p.sigma<=0)h.low=h.high=ts::NaN;
 h.valid=p.pedPresent&&std::isfinite(h.before)&&std::isfinite(h.after)&&h.factor>0&&raw.hg->at(j)<c.number("hg-max",2200)&&raw.hg->at(j)!=4000&&(i.sample||h.temp_measured);
 }}hits.Fill();++crossings;observed+=h.observed;response+=h.valid;
 }}events.Fill();
 }
 counts<<fi<<'\t'<<sampled<<'\t'<<paired<<'\t'<<missing<<'\t'<<reorder<<'\t'<<tailChecked<<'\t'<<tailRejected<<'\t'<<cogRejected<<'\t'<<selected<<'\t'<<crossings<<'\t'<<observed<<'\t'<<response<<'\t'<<pedRows<<'\n';
 provenance<<fi<<'\t'<<truthN<<'\t'<<truthChecked<<'\t'<<pdg<<'\t'<<primary<<'\t'<<(i.sample?"trigger_equals_time_generated_index_checked_against_current_MCDigi":"not_applicable")<<'\t'<<(i.sample?"species_energy_and_saved_ids_compatible_historical_chain_unproven":"saved_raw_calib_ids_verified")<<'\t'<<i.evidence<<"\tunconfirmed_not_inferred_from_path_or_mtime\n";
 std::cout<<"sampled="<<sampled<<" paired="<<paired<<" selected="<<selected<<" crossings="<<crossings<<" paired_response="<<response<<std::endl;
 }
 output.cd();events.Write();hits.Write();peds.Write();TNamed("schema",schema).Write();TNamed("tail_definition",tails?tails->definition().c_str():"none").Write();TNamed("temperature_definition","Current data EBUdecode inverse distance, 16 finite 0<T<60 C required for measured; missing not measured; Tref20 alpha1.6/135 or3.5/230; saved NewTemperature used for paired correction; MC nominal only").Write();events.SetDirectory(nullptr);hits.SetDirectory(nullptr);peds.SetDirectory(nullptr);output.Close();return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
