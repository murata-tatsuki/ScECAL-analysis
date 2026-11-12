#include "Common.hh"
#include <TH2D.h>
namespace {
using namespace spr;
struct Spectrum {Bin b;std::unique_ptr<TH1D>h;};
// ADC never enters these nuisance strata or deterministic selection ranks.
using Nuisance=std::tuple<int,int,int,int,int>; // run, T, sx, sy, other-coordinate
Nuisance nuisance(const Crossing&x,int axis,const mmr::Config&c){return {x.run,int(std::floor((x.temperature-20)/c.number("temperature-bin",2)+.5)),int(std::floor(x.sx/c.number("slope-bin",.02)+.5)),int(std::floor(x.sy/c.number("slope-bin",.02)+.5)),posBin(axis?x.longitudinal:x.transverse,1-axis,c)};}
int split(const Crossing&x,const mmr::Config&c){if(c.text("split-mode","event")=="run"){auto v=mmr::split(c.text("validation-runs"),',');return std::find(v.begin(),v.end(),std::to_string(x.run))!=v.end();}return identityHash(x.run,x.event_time,x.event)%2;}
using ConfKey=std::tuple<int,int,int,int,int>; // sample, condition, channel, axis, variable
std::string confName(const ConfKey&k){auto[s,c,id,a,v]=k;return "conf_s"+std::to_string(s)+"c"+std::to_string(c)+"i"+std::to_string(id)+"a"+std::to_string(a)+"v"+std::to_string(v);}
}
int main(int argc,char**argv)try{
 using namespace spr;require(argc==4,"Aggregate crossings.root configuration.tsv output-directory");mmr::Config cfg(argv[2]);std::string out=argv[3];require(!fs::exists(out+"/histograms.root"),"Histograms already exist");
 auto selected=selection(cfg.text("channels","all"));for(int axis=0;axis<2;++axis){require(bins(axis,cfg)>0&&bins(axis,cfg)%2==1,"Use odd bin counts so the central reference is a disjoint bin");require(2*half(axis)/bins(axis,cfg)>=2*cfg.number("position-resolution",.75),"Position bins narrower than twice the configured resolution floor");}
 std::array<std::map<int,mmr::Parameter>,2>pars;for(int s=0;s<2;++s){auto p=sample(s);pars[s]=mmr::loadParameters(cfg.text(p+"-pedestal"),cfg.text(p+"-mip"),cfg.text(p+"-hl"),cfg.text(p+"-threshold"));}
 auto file=mmr::openFile(argv[1]);auto*version=dynamic_cast<TParameter<int>*>(file->Get("strip_position_schema"));require(version&&version->GetVal()==schema,"Incompatible Crossings schema (legacy muon flats lack raw temperature/absence semantics)");auto*t=mmr::tree(*file,"Crossings");Crossing x;x.branches(*t,true);t->SetCacheSize(32*1024*1024);
 std::map<int,std::vector<Crossing>>channels;for(Long64_t i=0;i<t->GetEntries();++i){mmr::read(t,i);require(selected.count(x.cellid),"Crossings channel selection differs from saved config");channels[x.cellid].push_back(x);}t->ResetBranchAddresses();
 auto coverage=table(out+"/channel_coverage.tsv");coverage<<"cellid\tpredicted\tstatus\n";for(int id:selected)coverage<<id<<'\t'<<channels[id].size()<<'\t'<<(channels[id].empty()?"unmeasured":"collected")<<'\n';
 auto matching=table(out+"/matching.tsv");matching<<"cellid\taxis\tsplit\trun\tTbin\tsxbin\tsybin\tother_bin\tsample\tposition_bin\tpredicted\tposition_quota\tdata_mc_quota\n";
 TFile output((out+"/histograms.root").c_str(),"CREATE");require(!output.IsZombie(),"Cannot create histograms");Bin row;TTree bt("Bins","Absolute predictions/hits and moments; deterministic subsampling, not shape normalization");row.branches(bt);
 auto btout=table(out+"/bins.tsv");btout<<"sample\tcondition\tsplit\tcellid\taxis\tbin\tstage\tpredicted\tobserved\tunique\tusable\tmean_ADC\tmean_ADC_SE\tacquisition\tacquisition_SE\tduplicate\tforced4000\tsaturated\tswitched\tinvalid_temperature\tnegative_ADC\tzero_ADC\n";
 for(auto&[id,records]:channels){if(records.empty())continue;std::map<Key,Spectrum>spectra;std::map<ConfKey,std::unique_ptr<TH2D>>conf;
  for(int s=0;s<2;++s)for(int c=0;c<3;++c)for(int d=0;d<2;++d)for(int a=0;a<2;++a)for(int b=0;b<bins(a,cfg);++b)for(int st=0;st<3;++st){Key k{s,c,d,id,a,b,st};auto&p=pars[s].at(id);auto&sp=spectra[k];sp.b.k=k;sp.b.threshold=p.threshold;sp.b.sigma=p.sigma;sp.b.pedestal=p.ph;sp.b.threshold_fallback=p.thresholdFallback;}
  auto fill=[&](const Crossing&h,int axis,int condition){int b=posBin(axis?h.transverse:h.longitudinal,axis,cfg);if(b<0)return;int d=split(h,cfg);
   for(int st=0;st<3;++st){auto&sp=spectra.at({h.sample,condition,d,id,axis,b,st});auto&q=sp.b;++q.predicted;q.observed+=h.observed;q.unique+=h.multiplicity==1;q.duplicate+=h.multiplicity>1;q.forced+=h.raw_hg==4000;q.saturated+=h.raw_hg>=4095;q.switched+=h.raw_hg>=2600;q.invalid_temperature+=h.temperature_state!=0;
    if(h.usable){double v=st==0?h.raw_adc:st==1?h.path_adc:h.current_adc;double f=st==0?1:2/h.length*(st==2?mmr::tempFactor(id,h.temperature):1);require(std::isfinite(v)&&std::isfinite(f)&&f>0,"Invalid usable ADC");++q.usable;q.sum+=v;q.sum2+=v*v;q.negative+=v<0;q.zero+=v==0;q.fmin=std::min(q.fmin,f);q.fmax=std::max(q.fmax,f);if(!sp.h)sp.h=hist(q.k.name());sp.h->Fill(v);}
   }
   for(int var=0;var<4;++var){if(var==0&&h.temperature_state!=0)continue;ConfKey k{h.sample,condition,id,axis,var};auto&hh=conf[k];if(!hh){int n=var==0?80:var==1?400:50;double lo=var==0?0:var==1?0:-.05,hi=var==0?40:var==1?400:.05;hh=std::make_unique<TH2D>(confName(k).c_str(),"",bins(axis,cfg),-half(axis),half(axis),n,lo,hi);hh->SetDirectory(nullptr);}hh->Fill(axis?h.transverse:h.longitudinal,var==0?h.temperature:var==1?h.run:var==2?h.sx:h.sy);}
  };
  for(int axis=0;axis<2;++axis){using Lists=std::array<std::vector<std::vector<size_t>>,2>;std::map<std::pair<int,Nuisance>,Lists>groups;
   for(size_t i=0;i<records.size();++i){auto&h=records[i];fill(h,axis,0);int b=posBin(axis?h.transverse:h.longitudinal,axis,cfg);if(b<0||h.temperature_state!=0)continue;auto&g=groups[{split(h,cfg),nuisance(h,axis,cfg)}];for(auto&v:g)if(v.empty())v.resize(bins(axis,cfg));g[h.sample][b].push_back(i);}
   for(auto&[key,g]:groups){auto[d,n]=key;auto[run,tt,sx,sy,other]=n;size_t quota[2]={SIZE_MAX,SIZE_MAX};for(int s=0;s<2;++s)for(auto&v:g[s]){quota[s]=std::min(quota[s],v.size());std::sort(v.begin(),v.end(),[&](size_t a,size_t b){auto&x=records[a];auto&y=records[b];auto hx=identityHash(x.run,x.event_time,x.event),hy=identityHash(y.run,y.event_time,y.event);return hx==hy?a<b:hx<hy;});}size_t common=std::min(quota[0],quota[1]);
    for(int s=0;s<2;++s)for(int b=0;b<bins(axis,cfg);++b){auto&v=g[s][b];matching<<id<<'\t'<<axis<<'\t'<<d<<'\t'<<run<<'\t'<<tt<<'\t'<<sx<<'\t'<<sy<<'\t'<<other<<'\t'<<sample(s)<<'\t'<<b<<'\t'<<v.size()<<'\t'<<quota[s]<<'\t'<<common<<'\n';for(size_t i=0;i<quota[s];++i)fill(records[v[i]],axis,1);for(size_t i=0;i<common;++i)fill(records[v[i]],axis,2);}
   }
  }
  output.cd();for(auto&[key,sp]:spectra){row=sp.b;bt.Fill();if(sp.h)sp.h->Write();auto&q=sp.b;double eff=q.predicted?double(q.observed)/q.predicted:qnan;btout<<sample(key.sample)<<'\t'<<condition(key.condition)<<'\t'<<key.split<<'\t'<<id<<'\t'<<(key.axis?"v":"u")<<'\t'<<key.bin<<'\t'<<stage(key.stage)<<'\t'<<q.predicted<<'\t'<<q.observed<<'\t'<<q.unique<<'\t'<<q.usable<<'\t'<<mean(q)<<'\t'<<meanError(q)<<'\t'<<eff<<'\t'<<(q.predicted?std::sqrt(eff*(1-eff)/q.predicted):qnan)<<'\t'<<q.duplicate<<'\t'<<q.forced<<'\t'<<q.saturated<<'\t'<<q.switched<<'\t'<<q.invalid_temperature<<'\t'<<q.negative<<'\t'<<q.zero<<'\n';}for(auto&[k,h]:conf)h->Write();
 }
 output.cd();bt.Write();TParameter<int>("strip_position_schema",schema).Write();TNamed("reference","Central disjoint bin per axis, no interpolation. Relative response uses two independent fits except reference/reference=1 exactly. Training and validation split by saved event identity or explicit held-out runs.").Write();output.Close();
 auto summary=table(out+"/event_summary.tsv");summary<<"sample\tquantity\tevents\tmean\tstandard_error\tstatus\n";auto*ev=mmr::tree(*file,"Events");int s,rec,hit,phys,cal;double energy;mmr::bind(ev,"sample",&s);mmr::bind(ev,"recorded",&rec);mmr::bind(ev,"hit_records",&hit);mmr::bind(ev,"physical_strips",&phys);mmr::bind(ev,"calib_records",&cal);mmr::bind(ev,"energy_sum",&energy);std::array<std::array<std::array<long double,3>,5>,2>m{};
 for(Long64_t i=0;i<ev->GetEntries();++i){mmr::read(ev,i);double vals[]={double(rec),double(hit),double(phys),cal>=0?double(cal):qnan,energy};for(int j=0;j<5;++j)if(std::isfinite(vals[j])){auto&a=m[s][j];++a[0];a[1]+=vals[j];a[2]+=vals[j]*vals[j];}}
 const char*names[]={"recorded_raw","HitTag1_connected","physical_strips","recorded_calib","saved_energy_sum"};for(s=0;s<2;++s){for(int j=0;j<5;++j){auto&a=m[s][j];summary<<sample(s)<<'\t'<<names[j]<<'\t'<<a[0]<<'\t'<<(a[0]?double(a[1]/a[0]):qnan)<<'\t'<<(a[0]>1?std::sqrt(double(std::max(0.L,(a[2]-a[1]*a[1]/a[0])/(a[0]*(a[0]-1))))):qnan)<<'\t'<<(a[0]?"measured":"unmeasured")<<'\n';}summary<<sample(s)<<"\tSSA_virtual_hits\t0\tnan\tnan\tnot_read_muon_no_SSA\n";}
 std::cout<<"Aggregated channels="<<channels.size()<<"; raw inputs not rescanned per channel\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
