#include "Common.hh"
#include <TKey.h>
#include <TParameter.h>
namespace {
using namespace de;
using Key=std::tuple<int,int,int,int>;
using Reduced=std::tuple<int,int,int>;
struct Raw{
 TTree*t;int run=0,time=0,cycle=0,trigger=0;std::vector<int>*id=nullptr,*tag=nullptr;
 explicit Raw(TTree*p):t(p){bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,"CycleID",&cycle);bind(t,"TriggerID",&trigger);}
 void vectors(){bind(t,"CellID",&id);bind(t,"HitTag",&tag);t->SetCacheSize(16*1024*1024);for(auto n:{"CellID","HitTag"})t->AddBranchToCache(n,true);}
 Key key()const{return {run,time,cycle,trigger};} Reduced reduced()const{return {run,time,trigger};}
 ~Raw(){t->ResetBranchAddresses();delete id;delete tag;}
};
struct TrackAudit{double sigma=de::NaN;bool ambiguous=false;};
TrackAudit uncertainty(const mmr_track::Fit&f,const std::array<std::vector<double>,30>&fine,const Geometry&g,int parity,int excluded,double targetZ,const Settings&s){
 std::vector<double> zs;int ties=0;
 for(int l=parity;l<30;l+=2){if(l==excluded)continue;double best=1e30,second=1e30;for(double q:fine[l]){double d=std::abs(q-f.a-f.b*g.z[l]);if(d<best){second=best;best=d;}else second=std::min(second,d);}if(best>std::max(5.3,2*s.rms))continue;zs.push_back(g.z[l]);ties+=second-best<.25;}
 if(zs.size()<3)return {};
 double mean=std::accumulate(zs.begin(),zs.end(),0.)/zs.size(),ss=0;for(double z:zs)ss+=(z-mean)*(z-mean);
 // Quantization floor is a diagnostic approximation, not an external tracker resolution.
 double variance=std::max(5.3*5.3/12,f.rms*f.rms*zs.size()/(zs.size()-2));
 return {std::sqrt(variance*(1./zs.size()+(targetZ-mean)*(targetZ-mean)/ss)),ties>=2};
}
void schemaAudit(TFile&f,const std::string&path,std::ostream&o){TIter it(f.GetListOfKeys());while(auto*k=dynamic_cast<TKey*>(it())){o<<path<<"\tkey\t"<<k->GetName()<<'\t'<<k->GetClassName()<<'\n';auto*t=dynamic_cast<TTree*>(f.Get(k->GetName()));if(!t)continue;TIter b(t->GetListOfBranches());while(auto*p=dynamic_cast<TBranch*>(b()))o<<path<<'\t'<<t->GetName()<<'\t'<<p->GetName()<<'\t'<<p->GetTitle()<<'\n';}}
}
int main(int argc,char**argv)try{
 using namespace de;require(argc==4,"Collect manifest.tsv configuration.tsv output_directory");Config c(argv[2]);auto vs=variants(c);auto gs=groups(c);auto mask=masks(c);Geometry geo;std::string out=argv[3];fresh(out+"/crossings.root");
 auto audit=table(out+"/input_audit.tsv","file\tsample\trun\traw_entries\tcalib_entries\ttruth_entries\tsampled\tunique_raw\tduplicate_raw\tpaired_by_ID\tmissing_or_ambiguous_calib_ID\treordered_pairs\ttruth_muon_entries\ttruth_energy_MeV_min\ttruth_energy_MeV_max\tproduction_pairing\tCT_status");
 auto branches=table(out+"/schema.tsv","path\ttree_or_key\tname\ttitle_or_class");
 auto policies=table(out+"/sample_selection.tsv","sample\trun\telectron_tail_cut\ttail_rejected\telectron_energy_cut\tenergy_rejected\tCoG\ttemperature_correction\tpath_correction\tmask_channels\tnegative_zero_energy\tSSA_virtual_hits");
 auto settings=table(out+"/variants.tsv","variant\tname\tmin_layers_per_projection\tmax_angle_rad\tmax_RMS_mm\tlong_margin_mm\ttrans_margin_mm");for(size_t v=0;v<vs.size();++v){auto&s=vs[v].cuts;settings<<v<<'\t'<<vs[v].name<<'\t'<<s.minLayers<<'\t'<<s.angle<<'\t'<<s.rms<<'\t'<<s.longMargin<<'\t'<<s.transMargin<<'\n';}
 auto groupTable=table(out+"/groups.tsv","mode\texcluded_layers\tselection_bias");groupTable<<"0\ttarget_layer_only\tother_target_layers_used_in_tracking\n";for(size_t m=1;m<gs.size();++m){groupTable<<m<<'\t';for(int l:gs[m])groupTable<<l<<',';groupTable<<"\tall_listed_layers_excluded_before_track;conditioned_on_remaining_layers\n";}
 auto chipAudit=table(out+"/chip_readout_audit.tsv","sample\trun\tlayer\tchip\tindependent_readout_state\tindependent_multi_strip_predictions\tconclusion");
 TFile output((out+"/crossings.root").c_str(),"CREATE");require(!output.IsZombie(),"Cannot create crossings.root");Crossing x;TTree crossings("Crossings","All target attempts, including no track/no hit; no charge or energy selection");x.branches(crossings);
 int esample,erun,etime,ecycle,etrigger,raw_records,hit_records,physical_strips,unconnected_records,masked_records,calib_id_paired,ssa_virtual_hits=-1;Long64_t eentry;
 TTree events("Events","Decode record multiplicity before track selection; SSA unavailable");
#define B(n) events.Branch(#n,&n);
 B(esample) B(erun) B(etime) B(ecycle) B(etrigger) B(eentry) B(raw_records) B(hit_records) B(physical_strips) B(unconnected_records) B(masked_records) B(calib_id_paired) B(ssa_virtual_hits)
#undef B
 std::ifstream manifest(argv[1]);require(bool(manifest),"Missing manifest");int file=-1;std::set<std::pair<int,int>>sr;std::set<std::string> paths;
 for(std::string line;std::getline(manifest,line);){if(line.empty()||line[0]=='#')continue;auto f=split(line);require(f.size()==6,"Manifest needs six columns");int sample=std::stoi(f[0]),run=std::stoi(f[1]);require((sample==0||sample==1)&&sr.emplace(sample,run).second,"Repeated sample/run; combine split input explicitly before analysis");require(paths.insert(fs::canonical(f[2]).string()).second,"Repeated decode input");++file;
 auto rf=openFile(f[2]);schemaAudit(*rf,f[2],branches);Raw r(tree(*rf,"Raw_Hit"));std::map<Key,int> ids;std::map<Reduced,int> reduced;
 Long64_t entries=r.t->GetEntries();for(Long64_t i=0;i<entries;++i){read(r.t,i);require(sample==1||r.run==run,"Data run does not match manifest");++ids[r.key()];++reduced[r.reduced()];}
 std::map<Reduced,std::vector<Long64_t>> ci;Long64_t cn=-1,tn=-1,muons=0;double emin=NaN,emax=NaN;
 if(f[3]!="-"){auto cf=openFile(f[3]);schemaAudit(*cf,f[3],branches);auto*t=tree(*cf,"Calib_Hit");int cr=0,ct=0,ce=0;bind(t,"Run_Num",&cr);bind(t,"Event_Time",&ct);bind(t,"Event_Num",&ce);cn=t->GetEntries();for(Long64_t i=0;i<cn;++i){read(t,i);require(sample==1||cr==run,"Calibration run does not match manifest");ci[{cr,ct,ce}].push_back(i);}t->ResetBranchAddresses();}
 if(sample==1){require(f[5]!="-","MC requires explicit muon truth input for particle/energy audit");auto tf=openFile(f[5]);schemaAudit(*tf,f[5],branches);auto*t=tree(*tf,"MC_Truth");int pdg;double energy;bind(t,"PrimaryPDG",&pdg);bind(t,"PrimaryEnergy",&energy);tn=t->GetEntries();for(Long64_t i=0;i<tn;++i){read(t,i);require(std::abs(pdg)==13,"MC truth contains a non-muon primary");require(std::isfinite(energy)&&energy>0,"Invalid truth primary energy");++muons;if(i==0)emin=emax=energy;else{emin=std::min(emin,energy);emax=std::max(emax,energy);}}t->ResetBranchAddresses();require(muons>0,"Empty MC truth");require(std::abs(emin-c.number("energy",100)*1000)<1.&&std::abs(emax-c.number("energy",100)*1000)<1.,"Truth energy differs from requested GeV; specify actual muon energy");}
 int nm=0;for(auto&[key,reason]:mask)nm+=key.first==sample;
 policies<<sample<<'\t'<<run<<"\tdisabled_muon\t0\tdisabled_muon\t0\tnone\tnone_counts_only\tnone_counts_only\t"<<nm<<"\tnot_read_or_selected\tnot_measured_no_SSA_input\n";
 for(int l=0;l<30;++l)for(int chip=0;chip<6;++chip)chipAudit<<sample<<'\t'<<run<<'\t'<<l<<'\t'<<chip<<"\tunavailable_hit_derived_fields_only\t0_single_track_predictor\tconditional_chip_efficiency_only_no_chip_wide_loss_inference\n";
 r.vectors();Long64_t take=vs[0].cuts.maxEvents?std::min(entries,vs[0].cuts.maxEvents):entries,uniq=0,duplicates=0,paired=0,missing=0,reordered=0;
 for(Long64_t a=0;a<take;++a){Long64_t entry=(2*a+1)*entries/(2*take);read(r.t,entry);require(r.id&&r.tag&&r.id->size()==r.tag->size(),"Decode vector sizes differ");if(ids.at(r.key())!=1){++duplicates;continue;}++uniq;
 auto pair=ci.find(r.reduced());bool pairedID=reduced.at(r.reduced())==1&&pair!=ci.end()&&pair->second.size()==1;if(pairedID){++paired;reordered+=pair->second[0]!=entry;}else ++missing;
 std::map<int,std::vector<int>>hits;std::array<std::vector<double>,30>fine;raw_records=r.id->size();hit_records=physical_strips=unconnected_records=masked_records=0;
 for(size_t j=0;j<r.id->size();++j)if(r.tag->at(j)==1){++hit_records;int id=physical(r.id->at(j));if(!connected(id)){++unconnected_records;continue;}hits[id].push_back(r.id->at(j));if(mask.count({sample,id})){++masked_records;continue;}auto&s=geo.strips.at(id);fine[s.layer].push_back(s.layer%2?s.p[0]:s.p[1]);}
 physical_strips=hits.size();for(auto&q:fine){std::sort(q.begin(),q.end());q.erase(std::unique(q.begin(),q.end()),q.end());}
 esample=sample;erun=run;etime=r.time;ecycle=r.cycle;etrigger=r.trigger;eentry=entry;calib_id_paired=pairedID;events.Fill();
 for(size_t v=0;v<vs.size();++v){auto&s=vs[v].cuts;
 for(size_t m=0;m<gs.size();++m){auto work=fine;for(int l:gs[m])work[l].clear();auto fullX=mmr_track::track(work,geo.z,1,-1,s),fullY=mmr_track::track(work,geo.z,0,-1,s);
 std::vector<int> targets=gs[m];if(m==0){targets.resize(30);std::iota(targets.begin(),targets.end(),0);}
 for(int target:targets){x=Crossing{};x.sample=sample;x.run=run;x.source_run=r.run;x.file=file;x.time=r.time;x.cycle=r.cycle;x.trigger=r.trigger;x.entry=entry;x.variant=v;x.mode=m;x.layer=target;x.cellid=-1;
 auto fx=m?fullX:(target%2?mmr_track::track(work,geo.z,1,target,s):fullX);auto fy=m?fullY:(target%2?fullY:mmr_track::track(work,geo.z,0,target,s));
 x.track_ok=mmr_track::good(fx,s)&&mmr_track::good(fy,s);
 if(x.track_ok){x.x=fx.a+fx.b*geo.z[target];x.y=fy.a+fy.b*geo.z[target];x.sx=fx.b;x.sy=fy.b;x.nx=fx.n;x.ny=fy.n;x.rms=std::sqrt((fx.rms*fx.rms*fx.n+fy.rms*fy.rms*fy.n)/(fx.n+fy.n));
 auto ux=uncertainty(fx,work,geo,1,target,geo.z[target],s),uy=uncertainty(fy,work,geo,0,target,geo.z[target],s);x.sigma_x=ux.sigma;x.sigma_y=uy.sigma;x.track_ambiguous=ux.ambiguous||uy.ambiguous;
 int row=std::lround(20.5-(target%2?x.x:x.y)/5.3),col=std::lround(2-(target%2?x.y:x.x)/45.4);x.cellid=geo.cell(target,row,col);
 if(x.cellid>=0){auto&strip=geo.strips.at(x.cellid);x.chord=mmr_track::chord(strip,x.x,x.y,x.sx,x.sy,s,x.length,x.longitudinal,x.transverse);x.masked=mask.count({sample,x.cellid});
 double ls=target%2?x.sy:x.sx,ts=target%2?x.sx:x.sy,us=target%2?uy.sigma:ux.sigma,ut=target%2?ux.sigma:uy.sigma;
 x.ambiguous=!std::isfinite(us)||!std::isfinite(ut)||std::abs(x.longitudinal)+std::abs(ls)+c.number("ambiguity-sigma",1)*us>22.5||std::abs(x.transverse)+std::abs(ts)+c.number("ambiguity-sigma",1)*ut>2.5;
 auto hi=hits.find(x.cellid);if(hi!=hits.end()){x.multiplicity=hi->second.size();std::set<int>full(hi->second.begin(),hi->second.end());x.memory_extra=full.size()-1;x.repeated_records=x.multiplicity-full.size();}
 if(x.multiplicity==0)for(auto[dr,dc]:{std::pair<int,int>{-1,0},{1,0},{0,-1},{0,1}}){int neighbor=geo.cell(target,row+dr,col+dc);if(neighbor>=0&&!mask.count({sample,neighbor})&&hits.count(neighbor))x.neighbor_only=1;}
 x.eligible=x.chord&&!x.ambiguous&&!x.track_ambiguous&&!x.masked;
 }}
 crossings.Fill();
 }}}
 }
 audit<<file<<'\t'<<sample<<'\t'<<run<<'\t'<<entries<<'\t'<<cn<<'\t'<<tn<<'\t'<<take<<'\t'<<uniq<<'\t'<<duplicates<<'\t'<<paired<<'\t'<<missing<<'\t'<<reordered<<'\t'<<muons<<'\t'<<emin<<'\t'<<emax<<"\tdecode_calib_ID_checked_truth_has_no_ID_production_unverified\tunknown_no_historical_generation_receipt\n";
 std::cout<<"sample="<<sample<<" run="<<run<<" sampled="<<take<<" unique="<<uniq<<" ID_paired="<<paired<<std::endl;
 }
 require(file>=0,"No input files");output.cd();crossings.Write();events.Write();TNamed("schema",schema).Write();TNamed("configuration",slurp(argv[2]).c_str()).Write();TNamed("manifest",slurp(argv[1]).c_str()).Write();TNamed("definition","HitTag=1 raw acquisition; entire target layer/group excluded before reconstruction; no ADC, energy, CoG, tail, temperature or path response selection; geometry uncertainty diagnostic; not unconditional hardware efficiency; truth is provenance-only, no entry-order join").Write();output.Close();return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
