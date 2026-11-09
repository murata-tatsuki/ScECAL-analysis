#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TNamed.h>
#include <TH1.h>
#include <TProfile.h>
#include <TKey.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iomanip>
namespace fs=std::filesystem;
struct Stats {double n=0, temp=0, temp2=0, gain=0, eh=0,el=0, energy=0, raw4000=0, nraw=0;};
int main(int argc,char** argv) {
 if(argc!=3)return 1;TH1::AddDirectory(false);
 std::ofstream out(fs::path(argv[2])/"temperature_sample.tsv"), totals(fs::path(argv[2])/"exact_channel_energy.tsv");
 out<<std::setprecision(10)<<"energy\tsample\tcell\tn\tmean_temp\ttemp_rms\tmean_eh_over_qhg\traw_n\traw_hg4000\n";
 totals<<std::setprecision(10)<<"energy\tsample\tcell\tcal_hits\tmean_energy\tenergy_per_event\tmean_lg_adc\n";
 for(std::string energy:{"0.5","5","30","100","250"})for(std::string sample:{"data","mc"}) {
  auto name=fs::path(argv[1])/(sample=="data"?"data":"simulation/threshold")/(energy+"GeV.root");TFile histfile(name.c_str());
  auto* pairs=dynamic_cast<TNamed*>(histfile.Get("input_pairs"));if(!pairs)return 2;
  auto* nev=dynamic_cast<TParameter<Long64_t>*>(histfile.Get("processed_events"));
  std::map<int,Stats> stats;std::map<int,double> pedestal;
  for(int cell:{20008,220008,420008,620008,820008,920012,1020008,1220008,1420010,1820008,2820008}) {
   auto prefix="channel_"+std::to_string(cell)+"/";auto* ped=dynamic_cast<TParameter<double>*>(histfile.Get((prefix+"pedestal_hg").c_str()));if(!ped)continue;pedestal[cell]=ped->GetVal();
   auto* pe=dynamic_cast<TProfile*>(histfile.Get((prefix+"hg_selected_energy_profile").c_str()));auto* pl=dynamic_cast<TProfile*>(histfile.Get((prefix+"hg_lg_profile").c_str()));
   double n=0,sum=0,nl=0,sl=0;for(int b=0;b<=pe->GetNbinsX()+1;++b){n+=pe->GetBinEntries(b);sum+=pe->GetBinEntries(b)*pe->GetBinContent(b);}
   for(int b=0;b<=pl->GetNbinsX()+1;++b){nl+=pl->GetBinEntries(b);sl+=pl->GetBinEntries(b)*pl->GetBinContent(b);}
   totals<<energy<<'\t'<<sample<<'\t'<<cell<<'\t'<<n<<'\t'<<sum/n<<'\t'<<sum/nev->GetVal()<<'\t'<<sl/nl<<'\n';pe->SetDirectory(nullptr);pl->SetDirectory(nullptr);delete pe;delete pl;
  }
  std::istringstream manifest(pairs->GetTitle());std::string line;int files=0;
  while(files++<3&&std::getline(manifest,line)) {
   auto sep=line.find('\t');if(sep==std::string::npos)return 3;TFile rawf(line.substr(0,sep).c_str()), calf(line.substr(sep+1).c_str());
   auto* r=dynamic_cast<TTree*>(rawf.Get("Raw_Hit"));auto* c=dynamic_cast<TTree*>(calf.Get("Calib_Hit"));if(!r||!c||r->GetEntries()!=c->GetEntries())return 4;
   r->SetBranchStatus("*",0);c->SetBranchStatus("*",0);
   std::vector<int>* rid=nullptr,*tag=nullptr,*cid=nullptr;std::vector<double>*hg=nullptr,*eh=nullptr,*temp=nullptr;
   auto bind=[](TTree* t,const char* name,auto* ptr){t->SetBranchStatus(name,1);if(t->SetBranchAddress(name,ptr)<0)throw std::runtime_error(name);};
   int rr,rt,ri,cr,ct,ci;
   bind(r,"Run_Num",&rr);bind(r,"Event_Time",&rt);bind(r,"TriggerID",&ri);bind(c,"Run_Num",&cr);bind(c,"Event_Time",&ct);bind(c,"Event_Num",&ci);
   bind(r,"CellID",&rid);bind(r,"HitTag",&tag);bind(r,"HG_Charge",&hg);bind(c,"CellID",&cid);bind(c,"Hit_HG_Energy",&eh);bind(c,"NewTemperature",&temp);
   for(Long64_t e=0;e<std::min<Long64_t>(r->GetEntries(),1000);++e){r->GetEntry(e);c->GetEntry(e);if(rr!=cr||rt!=ct||ri!=ci)return 5;
    std::map<int,size_t> match;for(size_t j=0;j<cid->size();++j)match[(*cid)[j]]=j;
    for(size_t j=0;j<rid->size();++j){int id=(*rid)[j];int cell=id/100000*100000+(id%100000/10000)*10000+id%100;if(!pedestal.count(cell)||(*tag)[j]!=1)continue;
     auto& s=stats[cell];++s.nraw;if((*hg)[j]==4000)++s.raw4000;
     double q=(*hg)[j]-pedestal[cell];auto it=match.find(id);if(q<100||q>1500||it==match.end())continue;auto k=it->second;double t=temp->at(k),g=eh->at(k)/q;
     ++s.n;s.temp+=t;s.temp2+=t*t;s.gain+=g;
    }
   }
   r->ResetBranchAddresses();c->ResetBranchAddresses();delete rid;delete tag;delete cid;delete hg;delete eh;delete temp;
  }
  for(const auto&[cell,s]:stats)out<<energy<<'\t'<<sample<<'\t'<<cell<<'\t'<<s.n<<'\t'<<s.temp/s.n<<'\t'<<std::sqrt(std::max(0.,s.temp2/s.n-std::pow(s.temp/s.n,2)))<<'\t'<<s.gain/s.n<<'\t'<<s.nraw<<'\t'<<s.raw4000<<'\n';
  std::cout<<energy<<" "<<sample<<" checked"<<std::endl;
 }
}
