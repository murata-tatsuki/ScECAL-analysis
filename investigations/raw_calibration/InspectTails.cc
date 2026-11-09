#include "StudyCommon.hh"
#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <set>
#include <regex>
#include <unordered_map>
struct Row {std::map<std::string,double> v;};
std::vector<std::string> split(const std::string&s){std::istringstream i(s);std::vector<std::string>v;std::string x;while(std::getline(i,x,'\t'))v.push_back(x);return v;}
struct Acc {std::vector<double> residual;long long zero=0,tail=0;double sum=0;};
struct Event {Stats e,alternative,rejected,tailNhits,earlyMC;};
double quant(const std::vector<double>&x,double q){if(x.empty())return 0;double p=q*(x.size()-1);size_t k=p;return x[k]+(p-k)*(x[std::min(k+1,x.size()-1)]-x[k]);}
std::set<int> mask(const std::string&path){std::ifstream f(path);std::string t((std::istreambuf_iterator<char>(f)),{});auto p=t.find("const unordered_set<int> bad_channels = {");require(p!=std::string::npos,"Missing mask");p=t.find('{',p);t=t.substr(p,t.find("};",p)-p);std::regex re("[0-9]+");std::set<int>v;for(std::sregex_iterator it(t.begin(),t.end(),re),end;it!=end;++it)v.insert(std::stoi(it->str()));return v;}
int main(int argc,char**argv)try{
 require(argc==4,"InspectTails INPUT_DIAGNOSTIC_DIRECTORY COMPARISON_SOURCE IMPACT_STUDY_ROOT");fs::path out=argv[1];auto bad=mask(argv[2]);std::ifstream input(out/"channels.tsv");std::string line;std::getline(input,line);auto header=split(line);std::map<int,Row> rows;
 while(std::getline(input,line)){auto v=split(line);if(v[0]!="data")continue;Row row;for(size_t j=1;j<v.size();++j)row.v[header[j]]=std::stod(v[j]);rows[int(row.v.at("cell"))]=row;}
 std::map<int,std::array<Acc,2>> acc;std::map<std::tuple<int,int,int>,Acc> byrun;std::array<std::array<Event,2>,2> events;std::array<long long,2> ngain0{},nlogical{},nlogicalGain0{};std::array<double,2> logicalDiff{},logicalGain0Diff{},earlyDelta{};
 const std::array<int,5> targets={240008,420008,820008,920012,920018};
 std::map<int,std::array<std::unique_ptr<TH1D>,2>> plots;
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 for(int id:targets)for(int s=0;s<2;++s)plots[id][s]=std::make_unique<TH1D>((std::to_string(id)+(s?"mc":"data")).c_str(),("CellID "+std::to_string(id)+";LG residual from calibrated HL line [ADC];Probability / ADC").c_str(),240,-120,120);
 std::ofstream example(out/"tail_examples.tsv");example<<std::setprecision(12)<<"sample\tfile_index\tentry\trun\tcell\tgain_tag\tHG_raw\tLG_raw\tEH_MeV\tEL_MeV\tselected_MeV\tBCID\tmodal_event_BCID\n";
 std::map<int,long long> countExamples;std::array<long long,2> bcidMismatch{},bcidTailMismatch{};struct Drop {long long hits=0,zero=0,negative=0,lowHG=0,selectedLG=0,selectedNegative=0;double selectedEnergy=0;};
 std::array<Drop,2> drops;std::map<std::tuple<int,int,int>,Drop> dropChannels;
 std::ofstream dropped(out/"lg_dropout_examples.tsv");dropped<<"sample\tfile_index\tentry\tcell\tHG_raw\tLG_raw\tEH_MeV\tEL_MeV\tselected_MeV\tgain_tag\tBCID\tmodal_BCID\n";
 std::array<long long,2> nDropExamples{};
 std::ofstream timing(out/"event_readout.tsv");timing<<std::setprecision(12)<<"sample\tfile_index\tentry\trun\tmodal_BCID\tgood_energy_GeV\tshower_cog20\tgood_positive_mip_hits\toverlap_hits\toverlap_negative_tail\toverlap_positive_tail\tLG_selected_hits\tLG_selected_negative\tLG_selected_negative_energy_MeV\tmodal_LG_ADC\tLG_selected_near_ADC_mode\n";
 double reference=0; {auto f=openFile(argv[3]); auto*t=object<TTree>(*f,"paired_events");double e;bind(t,"nominal_GeV",&e);std::vector<double> v;for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);v.push_back(e);}std::sort(v.begin(),v.end());require(!v.empty(),"Empty MC reference");reference=v[v.size()/2];}
 std::ifstream manifest(out/"inputs/sampled_files.tsv");std::getline(manifest,line);
 while(std::getline(manifest,line)){auto fields=split(line);int s=fields[0]=="mc",fi=std::stoi(fields[1]);Long64_t count=std::stoll(fields[2]),take=std::stoll(fields[3]);auto raw=openFile(fields[4]),cal=openFile(fields[5]);auto*r=object<TTree>(*raw,"Raw_Hit"),*c=object<TTree>(*cal,"Calib_Hit");require(count==r->GetEntries()&&count==c->GetEntries(),"Entry mismatch");r->SetBranchStatus("*",0);c->SetBranchStatus("*",0);
  int rr,rt,re,cr,ct,ce;std::vector<int>*ri=nullptr,*tag=nullptr,*gain=nullptr,*ci=nullptr,*bcid=nullptr;std::vector<double>*rh=nullptr,*rl=nullptr,*eh=nullptr,*el=nullptr,*en=nullptr,*temp=nullptr,*x=nullptr,*y=nullptr;
  bind(r,"Run_Num",&rr);bind(r,"Event_Time",&rt);bind(r,"TriggerID",&re);bind(r,"CellID",&ri);bind(r,"HitTag",&tag);bind(r,"GainTag",&gain);bind(r,"HG_Charge",&rh);bind(r,"LG_Charge",&rl);bind(r,"BCID",&bcid);
  bind(c,"Run_Num",&cr);bind(c,"Event_Time",&ct);bind(c,"Event_Num",&ce);bind(c,"CellID",&ci);bind(c,"Hit_HG_Energy",&eh);bind(c,"Hit_LG_Energy",&el);bind(c,"Hit_Energy",&en);bind(c,"NewTemperature",&temp);bind(c,"Hit_X",&x);bind(c,"Hit_Y",&y);
  for(auto*t:{r,c}){t->SetCacheSize(16*1024*1024);t->AddBranchToCache("*",true);t->StopCacheLearningPhase();}
  int blocks=std::min<Long64_t>(10,take);
  for(int b=0;b<blocks;++b){Long64_t start=count*b/blocks,end=count*(b+1)/blocks,n=take/blocks+(b<take%blocks);start+=(end-start-n)/2;
   for(Long64_t entry=start;entry<start+n;++entry){r->GetEntry(entry);c->GetEntry(entry);require(rr==cr&&rt==ct&&re==ce,"Event key mismatch");std::unordered_map<int,size_t> index;std::map<int,int> bc;for(size_t k=0;k<ri->size();++k)if(tag->at(k)==1){require(index.emplace(ri->at(k),k).second,"Duplicate hit");++bc[bcid->at(k)];}int mode=-1,num=0;for(auto&[id,cnt]:bc)if(cnt>num){mode=id;num=cnt;}
    double goodE=0,pos=0,replacement=0,removed=0,early=0;int ntail=0;
    long long validHits=0,overlapN=0,negativeN=0,positiveN=0,lgN=0,lgNegative=0;double lgNegativeE=0;std::map<int,int> lgValues;std::array<double,30> le{},lx{},ly{};
    for(size_t j=0;j<ci->size();++j){int id=physical(ci->at(j));if(!connected(id))continue;require(index.count(ci->at(j)),"Missing raw hit");size_t k=index.at(ci->at(j));auto&p=rows.at(id).v;double h=rh->at(k),l=rl->at(k),qh=h-p.at("ped_HG"),ql=l-p.at("ped_LG"),mip=p.at("mip"),a=p.at("calib_HL_A"),offset=.305*p.at("calib_HL_B")/mip;
     if(en->at(j)>=0)pos+=en->at(j)/1000;
     if(mip>0&&!bad.count(id)){
      bool drop=ql-qh/a-(s?0:p.at("raw_HL_b")) < -std::max(15.,5*p.at("ped_sigma_LG"));
      ++validHits;if(h>=800&&h<2200){++overlapN;negativeN+=drop;positiveN+=ql-qh/a-(s?0:p.at("raw_HL_b"))>std::max(15.,5*p.at("ped_sigma_LG"));}
      if(h>=2600){++lgN;++lgValues[int(l)];if(en->at(j)<0){++lgNegative;lgNegativeE+=en->at(j);}}
      auto record=[&](Drop&v){++v.hits;v.zero+=l==0;v.negative+=el->at(j)<0;v.lowHG+=drop&&h<2600;v.selectedLG+=drop&&h>=2600;v.selectedNegative+=h>=2600&&en->at(j)<0;if(drop&&h>=2600)v.selectedEnergy+=en->at(j);};
      record(drops[s]);record(dropChannels[{s,fi,id}]);
      if(drop&&nDropExamples[s]++<100)dropped<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<entry<<'\t'<<id<<'\t'<<h<<'\t'<<l<<'\t'<<eh->at(j)<<'\t'<<el->at(j)<<'\t'<<en->at(j)<<'\t'<<gain->at(k)<<'\t'<<bcid->at(k)<<'\t'<<mode<<'\n';
     }
     if(mip<=0||bad.count(id)||en->at(j)<0)continue;
     goodE+=en->at(j)/1000;int layer=id/100000;le[layer]+=en->at(j);lx[layer]+=en->at(j)*x->at(j);ly[layer]+=en->at(j)*y->at(j);bool gain0=gain->at(k)==0;ngain0[s]+=gain0;bcidMismatch[s]+=bcid->at(k)!=mode;
     // Subtract each sample's linear ADC intercept. This isolates the tail shape from the known MC intercept omission.
     double residual=ql-qh/a-(s?0:p.at("raw_HL_b"));
     if(h>=800&&h<2200){auto&v=acc[id][s];v.residual.push_back(residual);v.zero+=gain0;if(!gain0)byrun[{s,fi,id}].residual.push_back(residual);if(plots.count(id))plots[id][s]->Fill(residual);}
     // Diagnostic flag only: strong positive HG/LG disagreement while the saved energy selects HG.
     double correctedLG=el->at(j)-(s?offset:0);bool suspect=h<2600&&correctedLG>2*eh->at(j)&&correctedLG-eh->at(j)>3;
     if(suspect){++nlogical[s];nlogicalGain0[s]+=gain0;double diff=(correctedLG-en->at(j))/1000;logicalDiff[s]+=diff;if(gain0)logicalGain0Diff[s]+=diff;bcidTailMismatch[s]+=bcid->at(k)!=mode;replacement+=diff;removed+=en->at(j)/1000;++ntail;if(countExamples[s]++<100)example<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<entry<<'\t'<<rr<<'\t'<<id<<'\t'<<gain->at(k)<<'\t'<<h<<'\t'<<l<<'\t'<<eh->at(j)<<'\t'<<el->at(j)<<'\t'<<en->at(j)<<'\t'<<bcid->at(k)<<'\t'<<mode<<'\n';}
     if(s&&h==4000&&ql*p.at("digi_HL_A")+p.at("ped_HG")<2600){early+=offset/1000;earlyDelta[s]+=offset/1000;}
    }
    bool central=pos>.5*reference&&pos<1.5*reference;for(int layer:{9,10})central=central&&le[layer]>0&&std::abs(lx[layer]/le[layer])<20&&std::abs(ly[layer]/le[layer])<20;
    int lgMode=-1,lgBest=0;for(auto&[a,n]:lgValues)if(n>lgBest){lgMode=a;lgBest=n;}long long lgNear=0;for(int a=lgMode-2;a<=lgMode+2;++a)lgNear+=lgValues[a];
    timing<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<entry<<'\t'<<rr<<'\t'<<mode<<'\t'<<goodE<<'\t'<<central<<'\t'<<validHits<<'\t'<<overlapN<<'\t'<<negativeN<<'\t'<<positiveN<<'\t'<<lgN<<'\t'<<lgNegative<<'\t'<<lgNegativeE<<'\t'<<lgMode<<'\t'<<lgNear<<'\n';
    for(int selection=0;selection<2;++selection)if(selection==0||central){auto&v=events[s][selection];v.e.add(goodE);v.alternative.add(goodE+replacement);v.rejected.add(goodE-removed);v.tailNhits.add(ntail);v.earlyMC.add(early);}
   }
  }
  r->ResetBranchAddresses();c->ResetBranchAddresses();delete ri;delete tag;delete gain;delete ci;delete bcid;delete rh;delete rl;delete eh;delete el;delete en;delete temp;delete x;delete y;
  std::cout<<(s?"mc":"data")<<" file "<<fi<<" tail diagnosis complete"<<std::endl;
 }
 std::ofstream outt(out/"robust_residuals.tsv");outt<<std::setprecision(12)<<"sample\tcell\tinner10um\tn\tgain0\tmedian_LG_ADC\tcentral68_sigma_LG_ADC\tRMS_LG_ADC\tpredicted_noise_LG_ADC\ttail5sigma\ttail_positive\ttail_negative\tcore3sigma_n\tcore3sigma_RMS\n";
 for(auto&[id,both]:acc)for(int s=0;s<2;++s){auto&v=both[s];auto&p=rows.at(id).v;auto&x=v.residual;std::sort(x.begin(),x.end());double med=quant(x,.5),sigma=std::sqrt(std::pow(p.at("ped_sigma_LG"),2)+std::pow(p.at("ped_sigma_HG")/p.at("calib_HL_A"),2)+1./12);Stats all,core;long long plus=0,minus=0;for(double z:x){all.add(z);if(std::abs(z-med)<3*sigma)core.add(z);plus+=z-med>5*sigma;minus+=z-med<-5*sigma;}outt<<(s?"mc":"data")<<'\t'<<id<<'\t'<<p.at("inner10um")<<'\t'<<x.size()<<'\t'<<v.zero<<'\t'<<med<<'\t'<<(quant(x,.84)-quant(x,.16))/2<<'\t'<<all.rms()<<'\t'<<sigma<<'\t'<<plus+minus<<'\t'<<plus<<'\t'<<minus<<'\t'<<core.n<<'\t'<<core.rms()<<'\n';}
 std::ofstream runs(out/"residuals_by_run.tsv");runs<<std::setprecision(12)<<"sample\tfile_index\tcell\tgain1_hits\tmedian_residual\tcentral68_sigma\tRMS\n";for(auto&[key,v]:byrun){auto [s,fi,id]=key;auto&x=v.residual;if(x.size()<20)continue;std::sort(x.begin(),x.end());Stats all;for(double a:x)all.add(a);runs<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<id<<'\t'<<x.size()<<'\t'<<quant(x,.5)<<'\t'<<(quant(x,.84)-quant(x,.16))/2<<'\t'<<all.rms()<<'\n';}
 std::ofstream ev(out/"tail_event_effect.tsv");ev<<std::setprecision(12)<<"sample\tselection\tevents\tmean_GeV\trms_GeV\tmean_LG_substitution_GeV\trms_LG_substitution_GeV\tmean_reject_suspect_GeV\trms_reject_suspect_GeV\tsuspect_hits_per_event\tearly_HG_switch_intercept_GeV\n";for(int s=0;s<2;++s)for(int k=0;k<2;++k){auto&v=events[s][k];ev<<(s?"mc":"data")<<'\t'<<(k?"shower_cog20":"all")<<'\t'<<v.e.n<<'\t'<<v.e.mean<<'\t'<<v.e.rms()<<'\t'<<v.alternative.mean<<'\t'<<v.alternative.rms()<<'\t'<<v.rejected.mean<<'\t'<<v.rejected.rms()<<'\t'<<v.tailNhits.mean<<'\t'<<v.earlyMC.mean<<'\n';}
 std::ofstream summary(out/"tail_counts.tsv");summary<<"sample\tgood_gain0_hits\tsuspect_hits\tsuspect_gain0\tsum_LG_substitution_GeV\tgood_BCID_not_event_mode\tsuspect_BCID_not_event_mode\n";for(int s=0;s<2;++s)summary<<(s?"mc":"data")<<'\t'<<ngain0[s]<<'\t'<<nlogical[s]<<'\t'<<nlogicalGain0[s]<<'\t'<<logicalDiff[s]<<'\t'<<bcidMismatch[s]<<'\t'<<bcidTailMismatch[s]<<'\n';
 std::ofstream dropfile(out/"lg_dropout_counts.tsv");dropfile<<"sample\tfile_index\tcell\tgood_channel_hits\tLG_raw_zero\tEL_negative\tnegative_residual_HG_selected\tnegative_residual_LG_selected\tselected_negative_LG\tnegative_residual_selected_energy_MeV\n";
 auto writeDrop=[&](int s,int fi,int id,const Drop&d){dropfile<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<id<<'\t'<<d.hits<<'\t'<<d.zero<<'\t'<<d.negative<<'\t'<<d.lowHG<<'\t'<<d.selectedLG<<'\t'<<d.selectedNegative<<'\t'<<d.selectedEnergy<<'\n';};
 for(int s=0;s<2;++s)writeDrop(s,-1,-1,drops[s]);
 for(auto&[key,v]:dropChannels){auto [s,fi,id]=key;writeDrop(s,fi,id,v);}
 TFile rf((out/"residual_tails.root").c_str(),"RECREATE");TCanvas can("tails","",1600,900);can.Divide(3,2);int pad=0;for(auto&[id,both]:plots){can.cd(++pad);gPad->SetLogy();double max=0;for(auto&h:both){h->Write();h->Scale(1./std::max(1.,h->GetEntries()));max=std::max(max,h->GetMaximum());}TLegend leg(.6,.75,.89,.89);for(int s=0;s<2;++s){auto&h=both[s];h->SetLineColor(s?kBlue+1:kBlack);h->SetMaximum(max*3);h->SetMinimum(1e-5);h->Draw(s?"HIST SAME":"HIST");leg.AddEntry(h.get(),s?"MC":"data","l");}leg.DrawClone();}can.SaveAs((out/"figures/hg_lg_residual_tails.png").c_str());return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
