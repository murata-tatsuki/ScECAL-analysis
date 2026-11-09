#include "../../analysis/resolution/HLTailSelection.hh"
#include "StudyCommon.hh"
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <unordered_map>

struct Event {
 int sample=0,run=0,bin=0; double e=0,central=0,hits=0,chits=0,cx=0,cy=0,noTemp=0,digiMip=0,cleanMip=0,validMip=0,lg=0,temp=0,negative=0;
 bool selected=false,tail=false; std::array<double,30> layers{}; std::array<double,4> cuts{};
};
struct Weighted {
 double w=0,w2=0,e=0,e2=0,c=0,n=0,cn=0,t=0,t2=0,m=0,m2=0,clean=0,clean2=0,v=0,lg=0,temp=0,negative=0;
 long long events=0; std::array<double,30> layers{}; std::array<std::array<double,30>,30> cross{};
 void add(const Event&a,double weight){
  ++events;w+=weight;w2+=weight*weight;e+=weight*a.e;e2+=weight*a.e*a.e;c+=weight*a.central;n+=weight*a.hits;cn+=weight*a.chits;
  t+=weight*a.noTemp;t2+=weight*a.noTemp*a.noTemp;m+=weight*a.digiMip;m2+=weight*a.digiMip*a.digiMip;clean+=weight*a.cleanMip;clean2+=weight*a.cleanMip*a.cleanMip;
  v+=weight*a.validMip;lg+=weight*a.lg;temp+=weight*a.temp;negative+=weight*a.negative;
  for(int i=0;i<30;++i){layers[i]+=weight*a.layers[i];for(int j=0;j<30;++j)cross[i][j]+=weight*a.layers[i]*a.layers[j];}
 }
};
struct Line {long long n=0;double h=0,l=0,hh=0,hl=0;Stats residual,difference;
 void add(double x,double y,double r,double d){++n;h+=x;l+=y;hh+=x*x;hl+=x*y;residual.add(r);difference.add(d);}
 double slope()const{return n>1&&hh-h*h/n>0?(hl-h*l/n)/(hh-h*h/n):0;}
};
struct Channel {Stats energy,hg,lg,temp;long long selectedLG=0,hg4000=0,lg3000=0,earlySwitch=0;std::array<Line,4> bands;};
const std::string header="sample\tselection\tevents\tsum_weight\teffective_events\tmean_GeV\tRMS_over_mean\tcentral_GeV\thits\tcentral_hits\ttemperature20_GeV\tdigiMip_GeV\tnonfallback_GeV\tcleanMip_GeV\ttemperature20_RMS_over_mean\tdigiMip_RMS_over_mean\tcleanMip_RMS_over_mean\tLG_GeV\tenergy_weighted_temperature\tnegative_all_connected_GeV\tSEM_GeV\n";
double resolution(double sum,double sum2,double w){double mean=sum/w;return mean>0?std::sqrt(std::max(0.,sum2/w-mean*mean))/mean:0;}
void write(std::ostream&f,int s,const std::string&name,const Weighted&v){
 if(!v.w)return;
 f<<std::setprecision(12)<<(s?"mc":"data")<<'\t'<<name<<'\t'<<v.events<<'\t'<<v.w<<'\t'<<v.w*v.w/v.w2<<'\t'<<v.e/v.w<<'\t'<<resolution(v.e,v.e2,v.w)<<'\t'<<v.c/v.w<<'\t'<<v.n/v.w<<'\t'<<v.cn/v.w<<'\t'<<v.t/v.w<<'\t'<<v.m/v.w<<'\t'<<v.v/v.w<<'\t'<<v.clean/v.w<<'\t'<<resolution(v.t,v.t2,v.w)<<'\t'<<resolution(v.m,v.m2,v.w)<<'\t'<<resolution(v.clean,v.clean2,v.w)<<'\t'<<v.lg/v.w<<'\t'<<v.temp/v.e<<'\t'<<v.negative/v.w<<'\t'<<std::sqrt(std::max(0.,v.e2/v.w-std::pow(v.e/v.w,2))*v.w2/(v.w*v.w))<<'\n';
}
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"input",""},{"out",""},{"audit",""},{"hl",""},{"mip",""},{"pedestal",""},{"mask",""},{"tail-results",""}});
 auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);
 std::ifstream mf(o["mask"]);std::string source((std::istreambuf_iterator<char>(mf)),{});auto a=source.find("const unordered_set<int> bad_channels = {");require(a!=std::string::npos,"Missing mask");a=source.find('{',a);source=source.substr(a,source.find("};",a)-a);std::regex numbers("[0-9]+");std::set<int> bad;for(std::sregex_iterator i(source.begin(),source.end(),numbers),end;i!=end;++i)bad.insert(std::stoi(i->str()));
 std::map<int,double>digi,cleanMip;std::ifstream af(o["audit"]);require(bool(af),"Missing MIP audit");std::string line;std::getline(af,line);while(std::getline(af,line)){std::istringstream row(line);int id,g,nd;double mp,chi,cal,d;require(bool(row>>id>>g>>mp>>chi>>nd>>cal>>d),"Bad MIP row");digi[id]=d;}
 {auto f=openFile(o["mip"]);auto*t=object<TTree>(*f,"MIP_Fit");int id,nd;double mp,chi;bind(t,"CellID",&id);bind(t,"NDF",&nd);bind(t,"LandauMPV",&mp);bind(t,"ChiSquare",&chi);double sum[2]={};int count[2]={};
  for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id<0||id/100000>=30||nd<=0||!std::isfinite(mp)||mp<=0||!std::isfinite(chi)||chi<0||chi/nd>(group(id)?2.:1.7))continue;sum[group(id)]+=mp;++count[group(id)];if(connected(id))cleanMip[id]=mp;}
  for(auto&[id,p]:pars)if(!cleanMip.count(id)||cleanMip[id]==1){require(count[group(id)]>0,"Empty clean MIP group");cleanMip[id]=sum[group(id)]/count[group(id)];}
 }
 fs::path input=o["input"],out=o["out"];auto referenceFile=openFile((input/"study.root").string());double beam=object<TParameter<double>>(*referenceFile,"beam_energy_GeV")->GetVal();std::string energy=std::to_string(int(beam))+"GeV";
 fs::path joint=out/"comparison"/"20mm"/energy;fs::create_directories(joint/"figures");
 auto* paired=object<TTree>(*referenceFile,"paired_events");double en;bind(paired,"nominal_GeV",&en);std::vector<double> energies;for(Long64_t i=0;i<paired->GetEntries();++i){paired->GetEntry(i);energies.push_back(en);}std::sort(energies.begin(),energies.end());require(!energies.empty(),"Empty MC");double ref=energies[energies.size()/2];
 std::unique_ptr<hl_tail::Selection> tail;if(!o["tail-results"].empty())tail=std::make_unique<hl_tail::Selection>(o["tail-results"]);
 std::vector<Event> events;std::map<int,std::array<int,2>> count;std::array<std::array<Weighted,4>,2> totals;std::map<std::pair<int,int>,Weighted> runs;
 std::array<std::map<int,Channel>,2> channel;std::array<std::array<Stats,4>,2> cutStats;std::array<std::array<Stats,8>,2> bandStats;std::array<std::map<int,Stats>,2> fixedHits;
 std::array<long long,2> sampled{},rejected{},checked{};std::array<double,2> maxH{},maxL{};
 std::ofstream detail(joint/"events.tsv");detail<<std::setprecision(12)<<"sample\tfile_index\tentry\trun\tselected\ttail_rich\tcog9x_mm\tcog10y_mm\tenergy_GeV\tcentral_GeV\thits\tcentral_hits\tenergy_temperature20_GeV\tenergy_digiMip_GeV\tenergy_nonfallback_GeV\tenergy_cleanMip_GeV\tLG_GeV\n";
 for(int s=0;s<2;++s){std::ifstream files(input/"inputs"/(s?"mc.tsv":"data.tsv"));require(bool(files),"Missing input manifest");int fi=-1;
  while(std::getline(files,line)){++fi;auto tab=line.find('\t');require(tab!=std::string::npos,"Bad file row");auto path=line.substr(tab+1);std::vector<unsigned char> flags;if(tail)flags=tail->mask(path);
   auto f=openFile(path),rf=openFile(line.substr(0,tab));auto*t=object<TTree>(*f,"Calib_Hit"),*raw=object<TTree>(*rf,"Raw_Hit");require(t->GetEntries()==raw->GetEntries(),"Entry count mismatch");t->SetBranchStatus("*",0);raw->SetBranchStatus("*",0);
   int run,time,trigger,rr,rt,re;std::vector<int>*ids=nullptr,*rid=nullptr,*tag=nullptr;std::vector<double>*signal=nullptr,*hg=nullptr,*lg=nullptr,*x=nullptr,*y=nullptr,*temp=nullptr,*rh=nullptr,*rl=nullptr;
   bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,"Event_Num",&trigger);bind(t,"CellID",&ids);bind(t,"Hit_Energy",&signal);bind(t,"Hit_HG_Energy",&hg);bind(t,"Hit_LG_Energy",&lg);bind(t,"Hit_X",&x);bind(t,"Hit_Y",&y);bind(t,"NewTemperature",&temp);
   bind(raw,"Run_Num",&rr);bind(raw,"Event_Time",&rt);bind(raw,"TriggerID",&re);bind(raw,"CellID",&rid);bind(raw,"HitTag",&tag);bind(raw,"HG_Charge",&rh);bind(raw,"LG_Charge",&rl);
   Long64_t n=t->GetEntries(),take=std::min<Long64_t>(n,s?1000:2000);int blocks=std::min<Long64_t>(10,take);
   for(int b=0;b<blocks;++b){Long64_t begin=n*b/blocks,end=n*(b+1)/blocks,len=take/blocks+(b<take%blocks);begin+=(end-begin-len)/2;
    for(Long64_t entry=begin;entry<begin+len;++entry){require(t->GetEntry(entry)>0&&raw->GetEntry(entry)>0,"Read failed");require(run==rr&&time==rt&&trigger==re,"Event identity mismatch");++sampled[s];Event e;e.sample=s;e.run=run;e.tail=tail&&flags.at(entry);std::array<double,30> lx{},ly{};double pos=0;
     require(ids->size()==signal->size()&&ids->size()==hg->size()&&ids->size()==lg->size()&&ids->size()==x->size()&&ids->size()==y->size()&&ids->size()==temp->size(),"Calib vector mismatch");require(rid->size()==rh->size()&&rid->size()==rl->size()&&rid->size()==tag->size(),"Raw vector mismatch");
     std::unordered_map<int,size_t> indices;indices.reserve(rid->size());for(size_t j=0;j<rid->size();++j)if(tag->at(j)==1)require(indices.emplace(rid->at(j),j).second,"Duplicate raw full CellID");
     for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!connected(id))continue;require(indices.count(ids->at(j)),"Calib hit absent from decode");size_t k=indices.at(ids->at(j));auto&p=pars.at(id);double C=1-(temp->at(j)-20)*(group(id)?1.6/135:3.5/230),qh=rh->at(k)-p.ph,ql=rl->at(k)-p.pl;
      double eh=.305/p.mip*C*qh,el=.305/p.mip*(C*ql*p.gain+p.intercept);maxH[s]=std::max(maxH[s],std::abs(hg->at(j)-eh));maxL[s]=std::max(maxL[s],std::abs(lg->at(j)-el));++checked[s];require(std::abs(hg->at(j)-eh)<1e-6&&std::abs(lg->at(j)-el)<1e-6,"Saved calibration differs from parameter model");bool isLG=rh->at(k)>=2600;require(std::abs(signal->at(j)-(isLG?el:eh))<1e-6,"Gain selection differs");
      double v=signal->at(j)/1000;require(std::isfinite(v),"Nonfinite energy");if(v<0){e.negative+=v;continue;}pos+=v;if(bad.count(id))continue;require(p.mip>0,"Unexpected nonpositive good MIP");++e.hits;e.e+=v;double offset=isLG?p.offset()/1000:0;
      require(C>0&&digi.at(id)>0,"Invalid sensitivity denominator");e.noTemp+=(v-offset)/C+offset;e.digiMip+=v*p.mip/digi.at(id);e.cleanMip+=v*p.mip/cleanMip.at(id);e.validMip+=p.mipFallback?0:v;e.lg+=isLG?v:0;e.temp+=v*temp->at(j);
      int layer=id/100000;e.layers[layer]+=v;lx[layer]+=x->at(j)*v;ly[layer]+=y->at(j)*v;if(std::abs(x->at(j))<22.5&&std::abs(y->at(j))<22.5){e.central+=v;++e.chits;}
      const std::array<double,4> cuts={0,.3,1,3};for(int c=0;c<4;++c)if(v*1000>=cuts[c])e.cuts[c]+=v;
     }
     e.selected=pos>.5*ref&&pos<1.5*ref;for(int l:{9,10})e.selected=e.selected&&e.layers[l]>0&&std::abs(lx[l]/e.layers[l])<20&&std::abs(ly[l]/e.layers[l])<20;
     e.cx=e.layers[9]>0?lx[9]/e.layers[9]:999;e.cy=e.layers[10]>0?ly[10]/e.layers[10]:999;
     detail<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<entry<<'\t'<<run<<'\t'<<e.selected<<'\t'<<e.tail<<'\t'<<e.cx<<'\t'<<e.cy<<'\t'<<e.e<<'\t'<<e.central<<'\t'<<e.hits<<'\t'<<e.chits<<'\t'<<e.noTemp<<'\t'<<e.digiMip<<'\t'<<e.validMip<<'\t'<<e.cleanMip<<'\t'<<e.lg<<'\n';
     if(!e.selected)continue;
     e.bin=int((e.cx+20)/5)*8+int((e.cy+20)/5);events.push_back(e);totals[s][0].add(e,1);if(e.tail){++rejected[s];continue;}
     ++count[e.bin][s];totals[s][1].add(e,1);runs[{s,run}].add(e,1);if(std::abs(e.cx)<5&&std::abs(e.cy)<5)totals[s][3].add(e,1);
     for(int c=0;c<4;++c){cutStats[s][c].add(e.cuts[c]);}fixedHits[s][int(e.chits)/2].add(e.central);
     for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!connected(id)||bad.count(id)||signal->at(j)<0)continue;size_t k=indices.at(ids->at(j));auto&p=pars.at(id);auto&ch=channel[s][id];double h=rh->at(k),l=rl->at(k),qh=h-p.ph,ql=l-p.pl;
      ch.energy.add(signal->at(j));ch.hg.add(qh);ch.lg.add(ql);ch.temp.add(temp->at(j));ch.selectedLG+=h>=2600;ch.hg4000+=h==4000;ch.lg3000+=l==3000;ch.earlySwitch+=h==4000&&ql*p.gain+p.intercept+p.ph<2600;
      int band=h<800?0:h<1400?1:h<2200?2:3;if(h>=300&&h<2600)ch.bands[band].add(qh,ql,ql-qh/p.gain+p.intercept/p.gain,lg->at(j)-hg->at(j));
      const std::array<double,7> edges={.1,.3,1,3,10,30,100};int eb=std::upper_bound(edges.begin(),edges.end(),signal->at(j))-edges.begin();bandStats[s][eb].add(signal->at(j)/1000);
     }
    }
   }t->ResetBranchAddresses();raw->ResetBranchAddresses();delete ids;delete signal;delete hg;delete lg;delete x;delete y;delete temp;delete rid;delete tag;delete rh;delete rl;
  }
 }
 std::ofstream weights(joint/"cog_weights.tsv");weights<<"bin\tdata_events\tmc_events\tcommon_weight\n";for(auto&[bin,n]:count)weights<<bin<<'\t'<<n[0]<<'\t'<<n[1]<<'\t'<<(std::min(n[0],n[1])>=20?std::min(n[0],n[1]):0)<<'\n';
 for(const auto&e:events)if(!e.tail){auto n=count.at(e.bin);if(std::min(n[0],n[1])>=20)totals[e.sample][2].add(e,double(std::min(n[0],n[1]))/n[e.sample]);}
 std::ofstream summary(joint/"summary.tsv"),runTable(joint/"runs.tsv");summary<<header;runTable<<header;
 const std::array<std::string,4> names={"keep","exclude","exclude_matched_cog5mm_bins_min20","exclude_cog_abs5mm"};
 for(int s=0;s<2;++s)for(int k=0;k<4;++k){write(summary,s,names[k],totals[s][k]);fs::path dest=out/(s?"mc":"data")/names[k]/(k==3?"5mm":"20mm")/energy;fs::create_directories(dest);std::ofstream f(dest/"summary.tsv");f<<header;write(f,s,names[k],totals[s][k]);}
 for(auto&[key,v]:runs)write(runTable,key.first,std::to_string(key.second),v);
 std::ofstream channels(joint/"channels.tsv");channels<<std::setprecision(12)<<"sample\tcell\thits\tmean_hit_MeV\tenergy_GeV\tmean_HG_sub\tmean_LG_sub\tmean_temperature\tLG_hits\tHG4000\tLG3000\tHG4000_implied_HG_lt2600\tcal_mip\tdigi_mip\tclean_mip\tfallback";for(int b=0;b<4;++b)channels<<"\tband"<<b<<"_n\tband"<<b<<"_residual_mean\tband"<<b<<"_residual_RMS\tband"<<b<<"_ELminusEH\tband"<<b<<"_slope";channels<<'\n';
 for(int s=0;s<2;++s)for(auto&[id,ch]:channel[s]){channels<<(s?"mc":"data")<<'\t'<<id<<'\t'<<ch.energy.n<<'\t'<<ch.energy.mean<<'\t'<<ch.energy.mean*ch.energy.n/1000<<'\t'<<ch.hg.mean<<'\t'<<ch.lg.mean<<'\t'<<ch.temp.mean<<'\t'<<ch.selectedLG<<'\t'<<ch.hg4000<<'\t'<<ch.lg3000<<'\t'<<ch.earlySwitch<<'\t'<<pars.at(id).mip<<'\t'<<digi.at(id)<<'\t'<<cleanMip.at(id)<<'\t'<<pars.at(id).mipFallback;for(auto&b:ch.bands)channels<<'\t'<<b.n<<'\t'<<b.residual.mean<<'\t'<<b.residual.rms()<<'\t'<<b.difference.mean<<'\t'<<b.slope();channels<<'\n';}
 std::ofstream cuts(joint/"hit_energy_cuts.tsv"),bands(joint/"hit_energy_bands.tsv"),fixed(joint/"fixed_central_hits.tsv"),cov(joint/"layer_covariance.tsv"),validation(joint/"validation.tsv");
 cuts<<"sample\tminimum_hit_MeV\tmean_GeV\tRMS_over_mean\n";bands<<"sample\tband\thits\tenergy_GeV\n";fixed<<"central_hits_min\tcentral_hits_max\tdata_n\tmc_n\tdata_mean_GeV\tmc_mean_GeV\n";cov<<"sample\tlayer_i\tlayer_j\tmean_i_GeV\tmean_j_GeV\tcovariance_GeV2\n";
 validation<<"sample\tsampled\tselected_keep\tselected_exclude\ttail_rejected\tchecked_hits\tmax_HG_error_MeV\tmax_LG_error_MeV\n";
 for(int s=0;s<2;++s){for(int c=0;c<4;++c){auto&v=cutStats[s][c];cuts<<(s?"mc":"data")<<'\t'<<std::array<double,4>{0,.3,1,3}[c]<<'\t'<<v.mean<<'\t'<<(v.mean?v.rms()/v.mean:0)<<'\n';}for(int b=0;b<8;++b)bands<<(s?"mc":"data")<<'\t'<<b<<'\t'<<bandStats[s][b].n<<'\t'<<bandStats[s][b].n*bandStats[s][b].mean<<'\n';
  auto&v=totals[s][1];for(int i=0;i<30;++i)for(int j=0;j<30;++j)cov<<(s?"mc":"data")<<'\t'<<i<<'\t'<<j<<'\t'<<v.layers[i]/v.w<<'\t'<<v.layers[j]/v.w<<'\t'<<v.cross[i][j]/v.w-v.layers[i]*v.layers[j]/(v.w*v.w)<<'\n';
  validation<<std::setprecision(12)<<(s?"mc":"data")<<'\t'<<sampled[s]<<'\t'<<totals[s][0].events<<'\t'<<v.events<<'\t'<<rejected[s]<<'\t'<<checked[s]<<'\t'<<maxH[s]<<'\t'<<maxL[s]<<'\n';
 }
 std::ofstream covarianceSummary(joint/"covariance_summary.tsv");covarianceSummary<<"sample\tselection\tevents\tmean_GeV\tRMS_over_mean\tnormalized_diagonal_variance\tnormalized_offdiagonal_covariance\n";
 auto writeCov=[&](int s,const std::string&name,const Weighted&v){if(!v.w)return;double sum=0,diag=0,total=0;for(int i=0;i<30;++i){sum+=v.layers[i]/v.w;for(int j=0;j<30;++j){double cv=v.cross[i][j]/v.w-v.layers[i]*v.layers[j]/(v.w*v.w);total+=cv;if(i==j)diag+=cv;}}double mean=v.e/v.w;require(std::abs(sum-mean)<1e-9&&std::abs(total-(v.e2/v.w-mean*mean))<1e-9,"Layer covariance does not reproduce event variance");covarianceSummary<<std::setprecision(12)<<(s?"mc":"data")<<'\t'<<name<<'\t'<<v.events<<'\t'<<mean<<'\t'<<std::sqrt(std::max(0.,total))/mean<<'\t'<<diag/(mean*mean)<<'\t'<<(total-diag)/(mean*mean)<<'\n';};
 for(int s=0;s<2;++s)for(int k=0;k<4;++k)writeCov(s,names[k],totals[s][k]);
 for(auto&[key,v]:runs)writeCov(key.first,"run"+std::to_string(key.second),v);
 for(auto&[n,d]:fixedHits[0]){auto m=fixedHits[1][n];if(d.n>=30&&m.n>=30)fixed<<n*2<<'\t'<<n*2+1<<'\t'<<d.n<<'\t'<<m.n<<'\t'<<d.mean<<'\t'<<m.mean<<'\n';}
 gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);TFile root((joint/"response.root").c_str(),"RECREATE");require(!root.IsZombie(),"Cannot create output ROOT");
 TCanvas canvas("response","",1400,550);canvas.Divide(2,1);
 for(int plot=0;plot<2;++plot){canvas.cd(plot+1);std::array<std::unique_ptr<TH1D>,2> hist;double lo=1e99,hi=-1e99;for(const auto&e:events)if(!e.tail){double v=plot?e.central:e.e;lo=std::min(lo,v);hi=std::max(hi,v);}double margin=(hi-lo)*.04;lo-=margin;hi+=margin;require(hi>lo,"Empty plot range");double ymax=0;
  for(int s=0;s<2;++s){std::string name=std::string(s?"mc":"data")+(plot?"_central":"_total");hist[s]=std::make_unique<TH1D>(name.c_str(),(energy+", tail excluded;"+(plot?"Central":"Total")+" positive good energy [GeV];Probability / bin").c_str(),120,lo,hi);hist[s]->Sumw2();for(const auto&e:events)if(e.sample==s&&!e.tail)hist[s]->Fill(plot?e.central:e.e);hist[s]->SetLineColor(s?kRed+1:kBlack);hist[s]->SetLineWidth(2);hist[s]->Scale(1/hist[s]->Integral());ymax=std::max(ymax,hist[s]->GetMaximum());hist[s]->Write();}
  hist[0]->SetMaximum(ymax*1.2);auto* drawnData=hist[0]->DrawCopy("HIST");auto* drawnMC=hist[1]->DrawCopy("HIST SAME");TLegend legend(.68,.74,.89,.88);legend.AddEntry(drawnData,"data","l");legend.AddEntry(drawnMC,"MC e_ssa","l");legend.DrawClone();
 }
 canvas.SaveAs((joint/"figures/event_energy.png").c_str());canvas.Clear();canvas.Divide(2,1);
 for(int plot=0;plot<2;++plot){canvas.cd(plot+1);std::array<TGraphErrors,2> graphs;double ymax=0;for(int s=0;s<2;++s){auto&v=totals[s][1];for(int i=0;i<30;++i){double mean=v.layers[i]/v.w,variance=v.cross[i][i]/v.w-mean*mean,value=plot?std::sqrt(std::max(0.,variance))/mean:mean;graphs[s].SetPoint(i,i,value);graphs[s].SetPointError(i,0,plot?0:std::sqrt(std::max(0.,variance)/v.w));ymax=std::max(ymax,value);}graphs[s].SetLineColor(s?kRed+1:kBlack);graphs[s].SetMarkerColor(s?kRed+1:kBlack);graphs[s].SetMarkerStyle(s?24:20);}
  gPad->DrawFrame(-.5,0,29.5,ymax*1.2,(energy+", tail excluded;Layer;"+(plot?"Layer RMS / mean":"Mean energy [GeV]")).c_str());TLegend legend(.65,.76,.88,.89);for(int s=0;s<2;++s){auto*obj=graphs[s].DrawClone("LP SAME");legend.AddEntry(obj,s?"MC e_ssa":"data","lp");graphs[s].Write((std::string(s?"mc":"data")+(plot?"_layer_relative_width":"_layer_mean")).c_str());}legend.DrawClone();
 }
 canvas.SaveAs((joint/"figures/layer_response.png").c_str());TNamed("tail_results",o["tail-results"].c_str()).Write();TNamed("selection","All-connected nonnegative energy in [0.5,1.5] MC signed median; good-hit CoG x/y in layers9,10 within20mm; original full-sample tail flags; fixed nominal nonnegative good hits for variants").Write();TParameter<double>("reference_GeV",ref).Write();root.Close();
 std::cout<<"Sensitivity complete: tail-excluded data="<<totals[0][1].events<<" MC="<<totals[1][1].events<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
