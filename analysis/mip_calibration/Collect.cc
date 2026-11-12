#include "Common.hh"
#include <EBUdecode.h>
#include <TProfile.h>
namespace {
using namespace mip;
struct Raw {
 TTree*t;int run=0,time=0,cycle=0,event=0;
 std::vector<int>*id=nullptr,*hit=nullptr,*gain=nullptr;std::vector<double>*hg=nullptr,*lg=nullptr;
 std::vector<std::vector<double>>*temp=nullptr;
 explicit Raw(TTree*p):t(p){for(auto n:{"Run_Num","Event_Time","CycleID","TriggerID"})require(t->GetBranch(n),std::string("Missing identity ")+n);mmr::bind(t,"Run_Num",&run);mmr::bind(t,"Event_Time",&time);mmr::bind(t,"CycleID",&cycle);mmr::bind(t,"TriggerID",&event);}
 void vectors(){mmr::bind(t,"CellID",&id);mmr::bind(t,"HitTag",&hit);mmr::bind(t,"GainTag",&gain);mmr::bind(t,"HG_Charge",&hg);mmr::bind(t,"LG_Charge",&lg);mmr::bind(t,"Temperature",&temp);t->SetCacheSize(32*1024*1024);t->AddBranchToCache("*",true);}
 auto key()const{return std::make_tuple(run,time,cycle,event);}
 void read(Long64_t i){mmr::read(t,i);require(id&&hit&&gain&&hg&&lg&&temp&&id->size()==hit->size()&&id->size()==gain->size()&&id->size()==hg->size()&&id->size()==lg->size(),"Decode vector size mismatch");}
 ~Raw(){t->ResetBranchAddresses();delete id;delete hit;delete gain;delete hg;delete lg;delete temp;}
};
struct Spectrum{Condition c;std::unique_ptr<TH1D>h;};
struct Dependency{Long64_t predicted=0,observed=0,n=0;long double raw=0,corrected=0,temp=0,length=0;};
struct Channel{
 Meta m;std::map<std::tuple<int,int,int>,Spectrum> conditions;
 std::unique_ptr<TH1D>raw,ped,path,current;
 explicit Channel(Meta v):m(v){raw=hist("raw_hg",1024,0,4096);ped=hist("pedestal_adc",1250,-100,2400);path=hist("path_adc",1250,-100,2400);current=hist("current_adc",1250,-100,2400);}
};
struct Crossing {
 int run=0,event_time=0,cycle=0,trigger=0,CellID=0,full_CellID=-1,observed=0,multiplicity=0,usable=0,temperature_state=0;
 Long64_t entry=0;double hg=mip::nan,lg=mip::nan,temperature=mip::nan,length=0,longitudinal=0,transverse=0,sx=0,sy=0,track_rms=0,raw_adc=mip::nan,path_adc=mip::nan,current_adc=mip::nan;
 void branches(TTree&t){
#define B(v) t.Branch(#v,&v)
 B(run);B(event_time);B(cycle);B(trigger);B(CellID);B(full_CellID);B(observed);B(multiplicity);B(usable);B(temperature_state);B(entry);B(hg);B(lg);B(temperature);B(length);B(longitudinal);B(transverse);B(sx);B(sy);B(track_rms);B(raw_adc);B(path_adc);B(current_adc);
#undef B
 }
};
}
int main(int argc,char**argv)try{
 using namespace mip;using namespace mmr_track;
 require(argc==4,"Usage: Collect manifest.tsv configuration.tsv output-directory");
 mmr::Config cfg(argv[2]);Settings settings(cfg);auto selected=selection(cfg.text("channels","all"));
 double factorStep=cfg.number("factor-step",.005);require(factorStep>0&&factorStep<=.02,"factor-step out of range");
 std::string out=argv[3];fresh(out+"/histograms.root");fresh(out+"/crossings.root");
 auto params=mmr::loadParameters(cfg.text("pedestal"),cfg.text("current-mip"),cfg.text("hl"),cfg.text("threshold"));
 std::map<int,Channel>channels;
 for(int id:mmr::channels()){auto&p=params.at(id);Meta m;m.CellID=id;m.selected=selected.count(id);m.pedestal=p.ph;m.ped_present=p.pedPresent;m.threshold=p.threshold;m.threshold_sigma=p.sigma;m.threshold_fallback=p.thresholdFallback;m.old_mip=p.mip;m.old_raw_mip=p.rawMip;m.old_fallback=p.mipFallback;channels.emplace(id,Channel(m));}
 std::array<std::array<std::array<int,5>,42>,30>grid;for(auto&l:grid)for(auto&r:l)r.fill(-1);
 std::array<double,30>z{};std::map<int,StripRectangle>strips;
 for(int id:mmr::channels()){StripRectangle s{id,id/100000,mmr::geometry(id)};double fine=s.layer%2?s.p[0]:s.p[1],along=s.layer%2?s.p[1]:s.p[0];int row=std::lround(20.5-fine/5.3),col=std::lround(2-along/45.4);require(row>=0&&row<42&&col>=0&&col<5&&grid[s.layer][row][col]<0,"Invalid strip grid");grid[s.layer][row][col]=id;strips.emplace(id,s);z[s.layer]=s.p[2];}
 std::ifstream manifest(argv[1]);require(bool(manifest),"Cannot read manifest");std::vector<std::pair<int,std::string>>files;std::set<std::string>seen;
 for(std::string line;std::getline(manifest,line);){if(line.empty()||line[0]=='#')continue;auto v=mmr::split(line);if(v[0]=="run")continue;require(v.size()==2,"Manifest requires run TAB absolute-decode-path");require(fs::path(v[1]).is_absolute()&&seen.insert(fs::canonical(v[1]).string()).second,"Repeated or relative input path");files.emplace_back(std::stoi(v[0]),v[1]);}
 require(!files.empty(),"No input files");
 std::map<std::tuple<int,int,int,int>,int>identities;
 std::ofstream audit(out+"/event_audit.tsv");audit<<"run\tfile\tentries\tsampled\tunique_events\tduplicate_identity_events\tinvalid_temperature_crossings\tcrossings\n";
 for(auto&[run,path]:files){auto f=mmr::openFile(path);Raw r(mmr::tree(*f,"Raw_Hit"));for(Long64_t i=0;i<r.t->GetEntries();++i){mmr::read(r.t,i);require(r.run==run,"Run_Num disagrees with manifest: "+path);++identities[r.key()];}}
 TFile flat((out+"/crossings.root").c_str(),"CREATE");require(!flat.IsZombie(),"Cannot create crossings.root");TTree crossings("Crossings","Decode-only, leave-target-layer-out predictions; including absent hits");Crossing x;x.branches(crossings);
 std::map<std::tuple<int,std::string,int>,Dependency>deps;
 for(auto&[run,path]:files){auto f=mmr::openFile(path);Raw r(mmr::tree(*f,"Raw_Hit"));r.vectors();Long64_t n=r.t->GetEntries(),take=settings.maxEvents?std::min(n,settings.maxEvents):n,unique=0,duplicates=0,ncross=0,badtemp=0;
 for(Long64_t sample=0;sample<take;++sample){Long64_t entry=settings.maxEvents?(2*sample+1)*n/(2*take):sample;r.read(entry);require(r.run==run,"Run changed during read");if(identities.at(r.key())!=1){++duplicates;continue;}++unique;
 std::map<int,std::vector<size_t>>hits;std::array<std::vector<double>,30>fine;
 for(size_t j=0;j<r.id->size();++j)if(r.hit->at(j)==1){int id=mmr::physical(r.id->at(j));if(!mmr::connected(id))continue;hits[id].push_back(j);auto&s=strips.at(id);fine[s.layer].push_back(s.layer%2?s.p[0]:s.p[1]);}
 for(auto&v:fine){std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}
 auto fullX=track(fine,z,1,-1,settings),fullY=track(fine,z,0,-1,settings);
 for(int l=0;l<30;++l){auto selectedLayer=selected.lower_bound(l*100000);if(selectedLayer==selected.end()||*selectedLayer>=(l+1)*100000)continue;auto fx=l%2?track(fine,z,1,l,settings):fullX,fy=l%2?fullY:track(fine,z,0,l,settings);if(!good(fx,settings)||!good(fy,settings))continue;
 double px=fx.a+fx.b*z[l],py=fy.a+fy.b*z[l];int row=std::lround(20.5-(l%2?px:py)/5.3),col=std::lround(2-(l%2?py:px)/45.4);if(row<0||row>=42||col<0||col>=5)continue;int id=grid[l][row][col];if(!selected.count(id))continue;
 double len,along,trans;auto&s=strips.at(id);if(!chord(s,px,py,fx.b,fy.b,settings,len,along,trans))continue;
 auto&c=channels.at(id);auto&m=c.m;++m.predicted;++ncross;x=Crossing{};x.run=run;x.event_time=r.time;x.cycle=r.cycle;x.trigger=r.event;x.CellID=id;x.entry=entry;x.length=len;x.longitudinal=along;x.transverse=trans;x.sx=fx.b;x.sy=fy.b;x.track_rms=std::sqrt((fx.rms*fx.rms*fx.n+fy.rms*fy.rms*fy.n)/(fx.n+fy.n));
 // Temperature at the SiPM/strip center, as in the current data reader.
 x.temperature_state=0;x.temperature=20;
 if(size_t(l)>=r.temp->size()||r.temp->at(l).empty()){x.temperature_state=1;++m.missing_temperature;}
 else {try{auto pos=s.p;x.temperature=tempReconstruction(l,pos.data(),r.temp->at(l));}catch(const std::exception&){x.temperature=mip::nan;x.temperature_state=2;++m.invalid_temperature;++badtemp;}}
 double factor=mmr::tempFactor(id,x.temperature);if(!(std::isfinite(factor)&&factor>0)){if(x.temperature_state!=2){++m.invalid_temperature;++badtemp;}x.temperature_state=2;}
 auto hit=hits.find(id);x.observed=hit!=hits.end();x.multiplicity=x.observed?hit->second.size():0;m.observed+=x.observed;
 if(x.multiplicity>1)++m.duplicate;
 if(x.multiplicity==1){++m.unique;auto j=hit->second[0];x.full_CellID=r.id->at(j);x.hg=r.hg->at(j);x.lg=r.lg->at(j);m.gain0+=r.gain->at(j)==0;
 bool finite=std::isfinite(x.hg)&&std::isfinite(x.lg);if(!finite)++m.invalid_adc;else{c.raw->Fill(x.hg);m.forced4000+=x.hg==4000;m.saturated+=x.hg>=4095;m.lg_saturated+=x.lg>=4095;m.switched+=x.hg>=2600;
 x.raw_adc=x.hg-m.pedestal;x.path_adc=x.raw_adc*2/len;x.current_adc=x.path_adc*factor;
 if(m.ped_present&&std::isfinite(m.pedestal)&&x.hg>=0&&x.hg<2600){x.usable=1;++m.usable;c.ped->Fill(x.raw_adc);c.path->Fill(x.path_adc);if(x.temperature_state!=2)c.current->Fill(x.current_adc);}
 }}
 for(int mode=0;mode<2;++mode){if(mode==1&&x.temperature_state!=0)continue;double fct=2/len*(mode==1?factor:1.);int bin=std::lround(fct/factorStep);double rounded=bin*factorStep;require(rounded>0,"Nonpositive binned correction factor");auto key=std::make_tuple(run,mode,bin);auto it=c.conditions.find(key);if(it==c.conditions.end()){Condition cc;cc.CellID=id;cc.run=run;cc.mode=mode;cc.bin=bin;cc.factor=rounded;cc.lo=cc.hi=fct;Spectrum a{cc,hist(cc.name(),1250,-100,2400)};it=c.conditions.emplace(key,std::move(a)).first;}auto&cond=it->second;++cond.c.predicted;cond.c.lo=std::min(cond.c.lo,fct);cond.c.hi=std::max(cond.c.hi,fct);if(x.usable)cond.h->Fill(x.raw_adc*rounded);}
 std::vector<std::pair<std::string,int>>bins={{"run",run},{"longitudinal_5mm",int(std::floor(along/5))},{"transverse_1mm",int(std::floor(trans))},{"path_0.001mm",int(std::floor(len/.001))}};if(x.temperature_state!=2)bins.emplace_back("temperature_1C",int(std::floor(x.temperature)));
 for(auto&[axis,b]:bins){auto&d=deps[{id,axis,b}];++d.predicted;d.observed+=x.observed;if(x.usable&&x.temperature_state!=2){++d.n;d.raw+=x.path_adc;d.corrected+=x.current_adc;d.temp+=x.temperature;d.length+=len;}}
 crossings.Fill();
 }
 }
 audit<<run<<'\t'<<path<<'\t'<<n<<'\t'<<take<<'\t'<<unique<<'\t'<<duplicates<<'\t'<<badtemp<<'\t'<<ncross<<'\n';std::cout<<"run="<<run<<" sampled="<<take<<" unique="<<unique<<" crossings="<<ncross<<std::endl;
 }
 flat.cd();crossings.Write();TNamed("identity","Run_Num,Event_Time,CycleID,TriggerID; duplicates across manifest excluded; full CellID retained for unique target hit; no cross-tree join").Write();flat.Close();
 TFile output((out+"/histograms.root").c_str(),"CREATE");require(!output.IsZombie(),"Cannot create histograms.root");Meta meta;Condition cond;TTree mt("Channels","Connected physical channels 0-29; counters and reference constants"),ct("Conditions","Exact quantized correction strata by run; no averaged threshold");meta.branches(mt);cond.branches(ct);
 for(auto&[id,c]:channels){meta=c.m;mt.Fill();auto*d=directory(output,channelName(id));d->cd();for(auto*h:{c.raw.get(),c.ped.get(),c.path.get(),c.current.get()})h->Write();auto*cd=directory(*d,"conditions");cd->cd();for(auto&[key,a]:c.conditions){cond=a.c;ct.Fill();a.h->Write();}}
 output.cd();mt.Write();ct.Write();TParameter<double>("factor_step",factorStep).Write();TNamed("definition","HG minus float-rounded production pedestal; full 45x5x2 mm strip chord; x=(HG-ped)*2/length*C; C=1-alpha*(T-20), alpha=1.6/135 (10um),3.5/230 (15um); C=1 diagnostic also saved. No crosstalk removal. Not a verified physical 20C extrapolation.").Write();output.Close();
 std::ofstream dep(out+"/dependencies.tsv");dep<<"CellID\taxis\tbin\tpredicted\tobserved\tusable\tmean_path_ADC\tmean_current_ADC\tmean_temperature\tmean_path_mm\n";for(auto&[key,d]:deps){auto[id,axis,bin]=key;dep<<std::setprecision(12)<<id<<'\t'<<axis<<'\t'<<bin<<'\t'<<d.predicted<<'\t'<<d.observed<<'\t'<<d.n<<'\t'<<(d.n?double(d.raw/d.n):mip::nan)<<'\t'<<(d.n?double(d.corrected/d.n):mip::nan)<<'\t'<<(d.n?double(d.temp/d.n):mip::nan)<<'\t'<<(d.n?double(d.length/d.n):mip::nan)<<'\n';}
 require(bool(audit)&&bool(dep),"Failed writing audit tables");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
