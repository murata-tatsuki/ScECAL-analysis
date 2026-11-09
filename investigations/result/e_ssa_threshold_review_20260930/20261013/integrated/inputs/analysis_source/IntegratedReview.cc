// Fresh saved-calibration diagnostics. Shares validated width/CV algorithms with the earlier study.
#define main historical_resolution_impact_main
#include "ResolutionImpact.cc"
#undef main
#include "../../analysis/resolution/HLTailSelection.hh"
#include "CurrentParameters.hh"
#include <TParameter.h>
#include <unordered_map>
#include <regex>
#include <set>
struct Detail {
 Event e{}; int eventTime=0,eventNumber=0; bool selected=false,tail=false,window=false;
 double temperature=0,negative=0,back5=0,edge=0,depth=0,allPositive=0;
 std::array<double,30> layers{};
 std::array<double,4> cutEnergy{},cutHits{};
};
struct Gain { Stats difference,residual; long long hg4000=0,lg3000=0,early=0,hits=0; };
struct Aggregate {
 Stats energy,central,hits,chits,temp,back,edge,depth,negative;
 std::array<Stats,30> layers;
 double cross[30][30]{};
 void add(const Detail&d){energy.add(d.e.e);central.add(d.e.central);hits.add(d.e.hits);chits.add(d.e.chits);temp.add(d.temperature);back.add(d.back5);edge.add(d.edge);depth.add(d.depth);negative.add(d.negative);for(int i=0;i<30;++i){layers[i].add(d.layers[i]);for(int j=0;j<30;++j)cross[i][j]+=d.layers[i]*d.layers[j];}}
};
void truthScale(const fs::path&out);
int main(int argc,char**argv)try {
 need(argc==3,"Usage: IntegratedReview INPUT_BASE OUTPUT_BASE, or --truth-only OUTPUT_BASE");if(std::string(argv[1])=="--truth-only"){truthScale(argv[2]);return 0;}fs::path input=argv[1],out=argv[2];fs::create_directories(out);
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 const fs::path project=fs::canonical(fs::path(__FILE__)).parent_path().parent_path().parent_path();
 const std::string params="/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/";
 auto pars=currentParameters(params+"all_hl_electron2023_hlratio_v3.root",params+"all_auto_muon_v4_trackfit.root",params+"pedestal2023_SPS.root");
 std::ifstream mf(project/"analysis/cog_channel_scan/SingleEnergyAnalysis.cc");std::string source((std::istreambuf_iterator<char>(mf)),{});auto a=source.find("const unordered_set<int> bad_channels = {");need(a!=std::string::npos,"Missing mask");a=source.find('{',a);source=source.substr(a,source.find("};",a)-a);std::regex digits("[0-9]+");std::set<int> bad;for(std::sregex_iterator i(source.begin(),source.end(),digits),end;i!=end;++i)bad.insert(std::stoi(i->str()));
 std::ofstream sum(out/"summary.tsv"),imp(out/"impact.tsv"),audit(out/"validation.tsv"),gain(out/"gain_channels.tsv"),mods(out/"models.tsv"),cov(out/"covariance.tsv"),prod(out/"production_resolution.tsv"),cuts(out/"hit_cuts.tsv");
 for(auto*f:{&sum,&imp,&audit,&gain,&mods,&cov,&prod,&cuts})*f<<std::setprecision(12);
 sum<<"energy_GeV\tsample\tselection\tcondition\tevents\tmean_GeV\tRMS_over_mean\tsigma68_over_median\tgaussian_sigma_over_mean\tfit_status\tbaseline_central_GeV\tbaseline_hits\tbaseline_central_hits\tbaseline_mean_temperature_C\tbaseline_back5_fraction\tbaseline_outer_fine_strip_fraction\tbaseline_depth_layer\tbaseline_negative_GeV\n";
 imp<<"energy_GeV\tselection\tcondition\tdata_events\tmc_events\tdata_RMS\tmc_RMS\tvariance_gap_reduction\tdata_gaussian\tmc_gaussian\tdata_mean_shift\tmc_mean_shift\n";
 audit<<"energy_GeV\tsample\tsampled\ttail_rejected\tchecked_hits\tmax_HG_error_MeV\tmax_LG_error_MeV\tmax_selected_error_MeV\n";
 gain<<"energy_GeV\tsample\tcell\tHG800_2200_hits\tEL_minus_EH_mean_MeV\tEL_minus_EH_RMS_MeV\tLG_line_residual_mean_ADC\tHG4000_hits\tLG3000_hits\tearly_switch_hits\thits\n";
 mods<<"energy_GeV\tsample_selection\tmode_0_run_1_position_2_both\theld_out_fold\tcoefficient_index\tcoefficient\n";
 cov<<"energy_GeV\tsample\tselection\tnormalized_diagonal_variance\tnormalized_offdiagonal_covariance\tRMS_over_mean\n";
 prod<<"energy_GeV\tsample\tcog_mm\tevent_stage\ttail_excluded\tevents\tmean_MeV\tRMS_over_mean\tgaussian_mean_MeV\tgaussian_sigma_over_mean\tmean_hits\n";
 cuts<<"energy_GeV\tsample\tselection\tminimum_hit_MeV\tmean_GeV\tRMS_over_mean\tmean_hits\n";
 for(int beam:{5,20,100}) {
  std::string energy=std::to_string(beam)+"GeV";fs::path dest=out/"comparison/20mm"/energy;fs::create_directories(dest);
  auto baseline=openFile((input/energy/"baseline/study.root").string());double ref=object<TParameter<double>>(*baseline,"beam_energy_GeV")->GetVal();need(ref==beam,"Beam mismatch");
  std::ifstream rf(input/energy/"baseline/reference.tsv");std::string line;std::getline(rf,line);std::getline(rf,line);std::istringstream rr(line);long long rn;double rmean;rmean=0;rr>>rn>>rmean>>ref;need(bool(rr),"Missing median");
  fs::path tailPath=project/"analysis/result/channel_response/hg_lg_tail/threshold/e_ssa"/energy;
  hl_tail::Selection tail(tailPath);std::array<std::vector<Detail>,2> all;
  std::array<std::map<int,Gain>,2> gains;
  for(int s=0;s<2;++s) {
   std::ifstream files(input/energy/"baseline/inputs"/(s?"mc.tsv":"data.tsv"));need(bool(files),"Missing input manifest");int fi=-1;long long sampled=0,rejected=0,checked=0;double maxH=0,maxL=0,maxS=0;
   while(std::getline(files,line)) {
    ++fi;auto tab=line.find('\t');need(tab!=std::string::npos,"Malformed manifest");auto rawpath=line.substr(0,tab),path=line.substr(tab+1);std::vector<unsigned char> flags;
    if(!s){need(fs::last_write_time(path)<=fs::last_write_time(tailPath/"tail_study.root"),"Stale data tail flags");flags=tail.mask(path);}
    auto f=openFile(path),rawFile=openFile(rawpath);auto*t=object<TTree>(*f,"Calib_Hit"),*raw=object<TTree>(*rawFile,"Raw_Hit");need(t->GetEntries()==raw->GetEntries(),"Raw/calib event count mismatch");t->SetBranchStatus("*",0);raw->SetBranchStatus("*",0);
    int run,time,trigger,rRun,rTime,rTrigger;std::vector<int>*ids=nullptr,*rid=nullptr,*tag=nullptr;std::vector<double>*signal=nullptr,*hg=nullptr,*lg=nullptr,*x=nullptr,*y=nullptr,*temp=nullptr,*rh=nullptr,*rl=nullptr;
    bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,"Event_Num",&trigger);bind(t,"CellID",&ids);bind(t,"Hit_Energy",&signal);bind(t,"Hit_HG_Energy",&hg);bind(t,"Hit_LG_Energy",&lg);bind(t,"Hit_X",&x);bind(t,"Hit_Y",&y);bind(t,"NewTemperature",&temp);
    bind(raw,"Run_Num",&rRun);bind(raw,"Event_Time",&rTime);bind(raw,"TriggerID",&rTrigger);bind(raw,"CellID",&rid);bind(raw,"HitTag",&tag);bind(raw,"HG_Charge",&rh);bind(raw,"LG_Charge",&rl);
    Long64_t n=t->GetEntries(),take=std::min<Long64_t>(n,s?1000:2000);int nb=std::min<Long64_t>(10,take);
    for(int b=0;b<nb;++b){Long64_t begin=n*b/nb,end=n*(b+1)/nb,len=take/nb+(b<take%nb);begin+=(end-begin-len)/2;
     for(Long64_t entry=begin;entry<begin+len;++entry){need(t->GetEntry(entry)>0&&raw->GetEntry(entry)>0,"Read failed");need(run==rRun&&time==rTime&&trigger==rTrigger,"Raw/calib identity differs");++sampled;Detail d;d.e.file=fi;d.e.entry=entry;d.e.run=run;d.eventTime=time;d.eventNumber=trigger;d.tail=!s&&flags.at(entry);rejected+=d.tail;
      need(ids->size()==signal->size()&&ids->size()==hg->size()&&ids->size()==lg->size()&&ids->size()==x->size()&&ids->size()==y->size()&&ids->size()==temp->size(),"Calib vectors differ");need(rid->size()==rh->size()&&rid->size()==rl->size()&&rid->size()==tag->size(),"Raw vectors differ");std::unordered_map<int,size_t> index;for(size_t k=0;k<rid->size();++k)if(tag->at(k)==1)need(index.emplace(rid->at(k),k).second,"Duplicate raw ID");
      std::array<double,30> lx{},ly{};double tw=0;
      for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!connected(id))continue;need(index.count(ids->at(j)),"Calib hit missing from raw");size_t k=index.at(ids->at(j));auto&p=pars.at(id);double C=1-(temp->at(j)-20)*(group(id)?1.6/135:3.5/230),qh=rh->at(k)-p.ph,ql=rl->at(k)-p.pl;double eh=.305/p.mip*C*qh,el=.305/p.mip*(C*ql*p.gain+p.intercept);maxH=std::max(maxH,std::abs(eh-hg->at(j)));maxL=std::max(maxL,std::abs(el-lg->at(j)));maxS=std::max(maxS,std::abs(signal->at(j)-(rh->at(k)<2600?hg->at(j):lg->at(j))));++checked;
       need(std::abs(eh-hg->at(j))<1e-6&&std::abs(el-lg->at(j))<1e-6&&maxS<1e-6,"Current calibration closure failed at "+path+" cell="+std::to_string(id)+" HGerr="+std::to_string(eh-hg->at(j))+" LGerr="+std::to_string(el-lg->at(j)));double v=signal->at(j)/1000;need(std::isfinite(v)&&C>0,"Invalid signal or temperature factor");if(v<0){d.negative+=v;continue;}d.allPositive+=v;if(bad.count(id))continue;d.e.e+=v;++d.e.hits;int l=id/100000;d.layers[l]+=v;lx[l]+=x->at(j)*v;ly[l]+=y->at(j)*v;tw+=v*temp->at(j);d.depth+=l*v;if(l>=25)d.back5+=v;double fine=l%2?x->at(j):y->at(j);if(std::abs(fine)>103.35)d.edge+=v;if(std::abs(x->at(j))<22.5&&std::abs(y->at(j))<22.5){d.e.central+=v;++d.e.chits;}bool isLG=rh->at(k)>=2600;double offset=isLG?p.offset()/1000:0;d.e.temp+=(v-offset)/C+offset;d.e.lg+=isLG?v:0;
       for(int c=0;c<4;++c)if(v*1000>=std::array<double,4>{0,.3,1,3}[c]){d.cutEnergy[c]+=v;++d.cutHits[c];}
      }
      d.selected=d.e.e>0;for(int l:{9,10})d.selected=d.selected&&d.layers[l]>0&&std::abs(lx[l]/d.layers[l])<20&&std::abs(ly[l]/d.layers[l])<20;d.window=d.allPositive>.5*ref&&d.allPositive<1.5*ref;
      d.e.x=d.layers[9]>0?lx[9]/d.layers[9]:999;d.e.y=d.layers[10]>0?ly[10]/d.layers[10]:999;if(d.e.e>0){d.temperature=tw/d.e.e;d.back5/=d.e.e;d.edge/=d.e.e;d.depth/=d.e.e;}
      if(d.selected&&!d.tail)for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!connected(id)||bad.count(id)||signal->at(j)<0)continue;size_t k=index.at(ids->at(j));auto&p=pars.at(id);auto&g=gains[s][id];++g.hits;g.hg4000+=rh->at(k)==4000;g.lg3000+=rl->at(k)==3000;g.early+=rh->at(k)==4000&&(rl->at(k)-p.pl)*p.gain+p.intercept+p.ph<2600;if(rh->at(k)>=800&&rh->at(k)<2200){g.difference.add(lg->at(j)-hg->at(j));g.residual.add(rl->at(k)-p.pl-(rh->at(k)-p.ph-p.intercept)/p.gain);}}
      all[s].push_back(d);
     }
    }
    t->ResetBranchAddresses();raw->ResetBranchAddresses();delete ids;delete signal;delete hg;delete lg;delete x;delete y;delete temp;delete rid;delete tag;delete rh;delete rl;
   }
   audit<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<sampled<<'\t'<<rejected<<'\t'<<checked<<'\t'<<maxH<<'\t'<<maxL<<'\t'<<maxS<<'\n';
   for(auto&[id,g]:gains[s])gain<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<id<<'\t'<<g.difference.n<<'\t'<<g.difference.mean<<'\t'<<g.difference.rms()<<'\t'<<g.residual.mean<<'\t'<<g.hg4000<<'\t'<<g.lg3000<<'\t'<<g.early<<'\t'<<g.hits<<'\n';
  }
  std::ofstream detail(dest/"events.tsv");
  detail<<std::setprecision(12)<<"sample\tfile\tentry\trun\tcog20\ttail_rich\tenergy_window\tx_mm\ty_mm\tenergy_GeV\tcentral_GeV\thits\tcentral_hits\tno_temperature_GeV\tback5_fraction\tedge_fraction\tdepth_layer\tmean_temperature_C\tlg_energy_GeV\tevent_time\tevent_number";
  for(int layer=0;layer<30;++layer)detail<<"\tlayer"<<layer<<"_GeV";
  detail<<'\n';
  for(int s=0;s<2;++s)for(auto&d:all[s]){
   detail<<(s?"mc":"data")<<'\t'<<d.e.file<<'\t'<<d.e.entry<<'\t'<<d.e.run<<'\t'<<d.selected<<'\t'<<d.tail<<'\t'<<d.window<<'\t'<<d.e.x<<'\t'<<d.e.y<<'\t'<<d.e.e<<'\t'<<d.e.central<<'\t'<<d.e.hits<<'\t'<<d.e.chits<<'\t'<<d.e.temp<<'\t'<<d.back5<<'\t'<<d.edge<<'\t'<<d.depth<<'\t'<<d.temperature<<'\t'<<d.e.lg<<'\t'<<d.eventTime<<'\t'<<d.eventNumber;
   for(double energy:d.layers)detail<<'\t'<<energy;
   detail<<'\n';
  }
  for(const std::string selection:{"cog20","cog20_energy_window"}){
   std::array<std::map<std::string,Values>,2> values;std::array<std::map<std::string,Width>,2> widths;std::array<Aggregate,2> aggregate;
   for(int s=0;s<2;++s){Events ev;for(auto&d:all[s])if(d.selected&&!d.tail&&(selection=="cog20"||d.window)){ev.push_back(d.e);aggregate[s].add(d);values[s]["baseline"].push_back(d.e.e);values[s]["no_temperature"].push_back(d.e.temp);if(d.back5<=.05)values[s]["back5_le_5percent"].push_back(d.e.e);if(d.edge<=.05)values[s]["outer_fine_strip_le_5percent"].push_back(d.e.e);for(int c=0;c<4;++c)values[s]["hitcut"+std::to_string(c)].push_back(d.cutEnergy[c]);}
    need(ev.size()>50,"Too few events");for(int mode=0;mode<3;++mode)values[s][std::array<std::string,3>{"run_cv","position_cv","position_run_cv"}[mode]]=corrected(ev,mode,mods,(s?"mc:":"data:")+selection,beam);
    for(auto&[condition,v]:values[s]){if(condition.rfind("hitcut",0)==0)continue;double mu=mean(v);TH1D h("diagnostic_energy","",240,.3*mu,1.7*mu);for(double x:v)h.Fill(x);auto w=width(v,&h);if(w.fit!=0)w.gaus=NAN;widths[s][condition]=w;auto&A=aggregate[s];sum<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<selection<<'\t'<<condition<<'\t'<<v.size()<<'\t'<<w.mean<<'\t'<<w.rms<<'\t'<<w.core<<'\t'<<w.gaus<<'\t'<<w.fit<<'\t'<<A.central.mean<<'\t'<<A.hits.mean<<'\t'<<A.chits.mean<<'\t'<<A.temp.mean<<'\t'<<A.back.mean<<'\t'<<A.edge.mean<<'\t'<<A.depth.mean<<'\t'<<A.negative.mean<<'\n';}
    for(int c=0;c<4;++c){auto&v=values[s]["hitcut"+std::to_string(c)];Stats count;for(auto&d:all[s])if(d.selected&&!d.tail&&(selection=="cog20"||d.window))count.add(d.cutHits[c]);cuts<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<selection<<'\t'<<std::array<double,4>{0,.3,1,3}[c]<<'\t'<<mean(v)<<'\t'<<rms(v)<<'\t'<<count.mean<<'\n';}
    auto&A=aggregate[s];double diag=0,total=0;for(int i=0;i<30;++i)for(int j=0;j<30;++j){double cv=A.cross[i][j]/A.energy.n-A.layers[i].mean*A.layers[j].mean;total+=cv;if(i==j)diag+=cv;}need(std::abs(total-A.energy.m2/A.energy.n)<1e-8,"Covariance closure failed");cov<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<selection<<'\t'<<diag/std::pow(A.energy.mean,2)<<'\t'<<(total-diag)/std::pow(A.energy.mean,2)<<'\t'<<A.energy.rms()/A.energy.mean<<'\n';
   }
   double gap=std::pow(widths[0]["baseline"].rms,2)-std::pow(widths[1]["baseline"].rms,2);for(auto&[condition,d]:widths[0]){auto&m=widths[1].at(condition);imp<<beam<<'\t'<<selection<<'\t'<<condition<<'\t'<<values[0][condition].size()<<'\t'<<values[1][condition].size()<<'\t'<<d.rms<<'\t'<<m.rms<<'\t'<<1-(d.rms*d.rms-m.rms*m.rms)/gap<<'\t'<<d.gaus<<'\t'<<m.gaus<<'\t'<<d.mean/widths[0]["baseline"].mean-1<<'\t'<<m.mean/widths[1]["baseline"].mean-1<<'\n';}
   if(selection=="cog20_energy_window"){TCanvas c("integrated_review","",2250,1000);c.SetCanvasSize(2250,1000);c.Divide(2,1);for(int p=0;p<2;++p){c.cd(p+1);std::string condition=p?"position_run_cv":"baseline";std::array<TH1D,2> h={TH1D("d","",100,.5,1.5),TH1D("m","",100,.5,1.5)};double ymax=0;for(int s=0;s<2;++s){auto&v=values[s][condition];double mu=mean(v);for(double x:v)h[s].Fill(x/mu);h[s].Scale(1./v.size());h[s].SetLineColor(s?kRed+1:kBlack);h[s].SetLineWidth(2);ymax=std::max(ymax,h[s].GetMaximum());}h[0].SetTitle((energy+", "+condition+";Energy / sample mean;Probability / bin").c_str());h[0].SetMaximum(1.15*ymax);auto*d=h[0].DrawCopy("HIST");auto*m=h[1].DrawCopy("HIST SAME");TLegend legend(.59,.76,.88,.88);legend.AddEntry(d,"data, LG tail excluded","l");legend.AddEntry(m,"MC e_ssa threshold","l");legend.DrawClone();}c.SaveAs((dest/"resolution_position_impact.png").c_str());}
  }
  for(int s=0;s<2;++s)for(int cog:{200,20,10,5}){fs::path p=project/"analysis/result/resolution"/(s?"simulation/threshold":"data")/(std::to_string(cog)+"mm")/(energy+"_"+std::to_string(cog)+"mm.root");auto f=openFile(p.string());auto*flag=object<TParameter<int>>(*f,"hl_tail_excluded");for(const std::string stage:{"beforeEventCut","afterEventCut"}){auto*h=object<TH1>(*f,stage+"/edep_"+energy);auto*fit=object<TF1>(*f,stage+"/fit_gaus");auto*nh=object<TH1>(*f,stage+"/nhits_"+energy);prod<<beam<<'\t'<<(s?"mc":"data")<<'\t'<<cog<<'\t'<<stage<<'\t'<<flag->GetVal()<<'\t'<<h->GetEntries()<<'\t'<<h->GetMean()<<'\t'<<h->GetRMS()/h->GetMean()<<'\t'<<fit->GetParameter(1)<<'\t'<<fit->GetParameter(2)/fit->GetParameter(1)<<'\t'<<nh->GetMean()<<'\n';}}
  std::cout<<energy<<" fresh review complete"<<std::endl;
 }
 truthScale(out);
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}

void truthScale(const fs::path&out) {
 fs::create_directories(out/"inputs");
 const std::string base="/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/";
 const std::string p=base+"Analysis_edit/share/";
 auto pars=currentParameters(p+"all_hl_electron2023_hlratio_v3.root",p+"all_auto_muon_v4_trackfit.root",p+"pedestal2023_SPS.root");
 std::ofstream result(out/"truth_scale.tsv"),metadata(out/"inputs/truth_files.tsv");
 result<<std::setprecision(12)<<"energy_GeV\tSiPM_um\ttruth_band_MeV\tHG_hits\tHG_sum_over_linear_prediction\tLG_hits\tLG_sum_over_linear_prediction\n";
 metadata<<"path\tbytes\tfilesystem_mtime_ticks\n";
 for(const std::string energy:{"5","20","100"}) {
  std::string rawpath=base+"Result_MC/decode/e-/e_ssa/threshold/"+energy+"GeV/"+energy+"GeV_0.root",truthpath=base+"Result_MC/generate/e-/e_ssa/"+energy+"GeV/"+energy+"GeV_0.root";
  auto fingerprint=[](const std::string&path){return std::make_pair(fs::file_size(path),fs::last_write_time(path));};auto beforeRaw=fingerprint(rawpath),beforeTruth=fingerprint(truthpath);
  auto rf=openFile(rawpath),tf=openFile(truthpath);auto*r=object<TTree>(*rf,"Raw_Hit"),*t=object<TTree>(*tf,"MC_Truth");need(r->GetEntries()==t->GetEntries(),"Truth/raw entry count differs");r->SetBranchStatus("*",0);t->SetBranchStatus("*",0);
  std::vector<int>*rid=nullptr,*tid=nullptr,*tag=nullptr;std::vector<double>*h=nullptr,*l=nullptr,*e=nullptr;int event;
  bind(r,"CellID",&rid);bind(r,"HitTag",&tag);bind(r,"HG_Charge",&h);bind(r,"LG_Charge",&l);bind(r,"Event_Time",&event);bind(t,"CellID",&tid);bind(t,"Hit_Energy",&e);
  double sumH[2][3]{},predH[2][3]{},sumL[2][3]{},predL[2][3]{};long long nh[2][3]{},nl[2][3]{};
  for(Long64_t i=0;i<r->GetEntries();++i){need(r->GetEntry(i)>0&&t->GetEntry(i)>0&&event==i,"Truth/raw event identity differs");need(rid->size()==tag->size()&&rid->size()==h->size()&&rid->size()==l->size()&&tid->size()==e->size(),"Truth/raw vectors differ");std::map<int,size_t> index;for(size_t j=0;j<rid->size();++j)if(tag->at(j)==1)need(index.emplace(rid->at(j),j).second,"Duplicate truth/raw key");
   for(size_t j=0;j<tid->size();++j){int id=physical(tid->at(j));if(!connected(id)||!index.count(tid->at(j)))continue;double en=e->at(j);int b=en>=1&&en<2?0:en>=30&&en<100?1:en>=100?2:-1;if(b<0)continue;auto&p=pars.at(id);size_t k=index.at(tid->at(j));int g=group(id);double pred=en/.305*p.mip;if(l->at(k)<3000){sumL[g][b]+=(l->at(k)-p.pl)*p.gain+p.intercept;predL[g][b]+=pred;++nl[g][b];}if(h->at(k)<2400){sumH[g][b]+=h->at(k)-p.ph;predH[g][b]+=pred;++nh[g][b];}}
  }
  for(int g=0;g<2;++g)for(int b=0;b<3;++b)result<<energy<<'\t'<<(g?10:15)<<'\t'<<std::array<std::string,3>{"1-2","30-100",">=100"}[b]<<'\t'<<nh[g][b]<<'\t'<<(predH[g][b]?sumH[g][b]/predH[g][b]:NAN)<<'\t'<<nl[g][b]<<'\t'<<(predL[g][b]?sumL[g][b]/predL[g][b]:NAN)<<'\n';
  need(beforeRaw==fingerprint(rawpath)&&beforeTruth==fingerprint(truthpath),"Truth/raw input changed while reading");for(const auto&path:{rawpath,truthpath})metadata<<fs::canonical(path).string()<<'\t'<<fs::file_size(path)<<'\t'<<fs::last_write_time(path).time_since_epoch().count()<<'\n';
  r->ResetBranchAddresses();t->ResetBranchAddresses();delete rid;delete tid;delete tag;delete h;delete l;delete e;
 }
}
