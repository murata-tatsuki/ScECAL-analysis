// Raw HG spectra; no production calibration or gain constants are applied.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TNamed.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "../../channel_response/temperature_calibration/TemperatureGeometry.hh"
using Key=std::tuple<int,int,int,int,int,int>;
struct Spectrum {
  std::array<unsigned long long,4096> counts{};
  unsigned long long n=0,nt=0,overflow=0;
  double t=0,tt=0;
};
void need(bool ok,const std::string&s){if(!ok)throw std::runtime_error(s);}
template<class T> void bind(TTree*t,const char*n,T*p){
 need(t->GetBranch(n),std::string("Missing branch ")+n); t->SetBranchStatus(n,1);
 need(t->SetBranchAddress(n,p)>=0,std::string("Invalid branch type ")+n);
}
#include "LEDGroups.hh"
int main(int argc,char**argv)try{
 need(argc==8,"Collect RAW OUTPUT CHANNELS|all MAX_EVENTS TEMP_BIN_WIDTH LED_GROUP|-1 TIME_BIN_SECONDS");
 std::string raw=argv[1],out=argv[2],selection=argv[3];
 long long limit=std::stoll(argv[4]);double tw=std::stod(argv[5]),timeWidth=std::stod(argv[7]);
 int group=std::stoi(argv[6]);
 need(limit>=0&&tw>0&&timeWidth>=0,"Invalid collection options");
 need(group>=-1&&group<14,"Invalid LED group");
 std::set<int> selected;
 if(selection!="all"){std::istringstream s(selection);std::string v;while(std::getline(s,v,','))selected.insert(std::stoi(v));need(!selected.empty(),"Empty channels");}
 std::unique_ptr<TFile>f(TFile::Open(raw.c_str(),"READ"));
 need(f&&!f->IsZombie()&&!f->TestBit(TFile::kRecovered),"Unreadable/recovered ROOT");
 auto*t=dynamic_cast<TTree*>(f->Get("Raw_Hit"));need(t,"Missing Raw_Hit");
 t->SetBranchStatus("*",0);
 int run=0,time=0;std::vector<int>*ids=nullptr,*hits=nullptr,*gains=nullptr;
 std::vector<double>*hg=nullptr;std::vector<std::vector<double>>*temps=nullptr;
 bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,"CellID",&ids);
 bind(t,"HitTag",&hits);bind(t,"GainTag",&gains);bind(t,"HG_Charge",&hg);bind(t,"Temperature",&temps);
 using FullKey=std::tuple<int,int,int,int,int,int,int>;
 std::map<FullKey,Spectrum> spectra;
 // Memory diagnostics remain separate even when signal spectra are pooled.
 using MemKey=std::tuple<int,int,int,int,int>;
 std::map<MemKey,std::array<double,3>> diagnostics;
 Long64_t total=t->GetEntries(),n=limit?std::min<Long64_t>(total,limit):total;
 unsigned long long disconnected=0,other=0; std::set<int>runs;
 for(Long64_t j=0;j<n;++j){
  Long64_t entry=limit?j*total/n:j; need(t->GetEntry(entry)>0,"Read failure");
  need(ids&&hits&&gains&&hg&&temps&&ids->size()==hits->size()&&ids->size()==gains->size()&&ids->size()==hg->size(),"Raw vector size mismatch");
  int timeBin=timeWidth?static_cast<int>(std::floor(time/timeWidth)):0;
  runs.insert(run);
  for(size_t h=0;h<ids->size();++h){
   int full=ids->at(h);need(full>=0,"Negative CellID");
   int l=full/100000,p=full%100000/10000,c=full%100,mem=full%10000/100;
   need(l<32&&p<6&&c<36&&mem<16,"Malformed CellID");
   if(l>=30||(p==5&&c>=30)){++disconnected;continue;}
   if(mem!=0){++other;continue;}
   int id=l*100000+p*10000+c;
   if(!selected.empty()&&!selected.count(id))continue;
   if(group>=0&&!gain::ledChannel(group,p*36+c)){++other;continue;}
   int ht=hits->at(h),gt=gains->at(h);need((ht==0||ht==1)&&(gt==0||gt==1),"Unexpected tag");
   double q=hg->at(h);need(std::isfinite(q),"Nonfinite HG");
   int bin=std::lround(q);need(std::abs(q-bin)<1e-6,"Noninteger HG: ADC-bin convention must be updated");
   double temp=NAN;
   if(static_cast<int>(temps->size())>l&&temps->at(l).size()==16){
    const auto&v=temps->at(l);bool valid=true;for(double x:v)valid &= std::isfinite(x)&&x>0&&x<60;
    if(valid){auto*pos=EBUdecode(l,p,c);temp=tempReconstruction(l,pos,v);}
   }
   int tb=std::isfinite(temp)?static_cast<int>(std::floor(temp/tw)):-1;
   auto&a=spectra[{run,id,ht,gt,tb,mem,timeBin}];++a.n;
   if(std::isfinite(temp)){++a.nt;a.t+=temp;a.tt+=temp*temp;}
   if(bin>=0&&bin<4096){++a.counts[bin];}else ++a.overflow;
   if(ht==0){auto&d=diagnostics[{run,id,gt,mem,tb}];d[0]++;d[1]+=q;d[2]+=q*q;}
  }
 }
 std::unique_ptr<TFile>o(TFile::Open(out.c_str(),"CREATE"));need(o&&!o->IsZombie(),"Output already exists or cannot be created");
 int sid=0,id,ht,gt,tb,mem,tbin;unsigned long long entries,nt,overflow;double tm,trms;
 TTree meta("Spectra","Uncorrected spectra; physical CellID=layer*100000+chip*10000+channel");
 meta.Branch("sid",&sid);meta.Branch("run",&run);meta.Branch("cellid",&id);meta.Branch("hit_tag",&ht);meta.Branch("gain_tag",&gt);
 meta.Branch("temperature_bin",&tb);meta.Branch("memory_cell",&mem);meta.Branch("time_bin",&tbin);
 meta.Branch("entries",&entries);meta.Branch("valid_temperature_entries",&nt);meta.Branch("overflow",&overflow);
 meta.Branch("temperature_C",&tm);meta.Branch("temperature_RMS_C",&trms);
 for(auto&[key,a]:spectra){
  std::tie(run,id,ht,gt,tb,mem,tbin)=key;entries=a.n;nt=a.nt;overflow=a.overflow;
  tm=nt?a.t/nt:NAN;trms=nt?std::sqrt(std::max(0.,a.tt/nt-tm*tm)):NAN;meta.Fill();
  TH1D hh(("s"+std::to_string(sid)).c_str(),"Raw HG;HG ADC;entries",4096,-.5,4095.5);
  for(int b=0;b<4096;++b){hh.SetBinContent(b+1,a.counts[b]);hh.SetBinError(b+1,std::sqrt(a.counts[b]));}
  hh.SetEntries(entries);hh.Write();
  ++sid;
 }
 meta.Write();
 TNamed mode("bootstrap_storage","disabled; memory cell 0 only");mode.Write();
 double dn,dm,ds;TTree md("MemoryDiagnostics","HitTag0 conditional moments, not an unbiased pedestal measurement");
 md.Branch("run",&run);md.Branch("cellid",&id);md.Branch("gain_tag",&gt);md.Branch("memory_cell",&mem);md.Branch("temperature_bin",&tb);
 md.Branch("entries",&dn);md.Branch("mean_ADC",&dm);md.Branch("RMS_ADC",&ds);
 for(auto&[k,a]:diagnostics){std::tie(run,id,gt,mem,tb)=k;dn=a[0];dm=a[1]/dn;ds=std::sqrt(std::max(0.,a[2]/dn-dm*dm));md.Fill();}md.Write();
 TNamed("sampling",("raw="+raw+"; total_events="+std::to_string(total)+"; sampled_events="+std::to_string(n)+"; memory_cell=0; bootstrap=0").c_str()).Write();
 o->Write();o->Close();
 std::cout<<"events="<<n<<" spectra="<<spectra.size()<<" disconnected_hits="<<disconnected<<" outside_LED_group="<<other<<"\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
