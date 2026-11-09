#include <TFile.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TNamed.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <iomanip>
#include <stdexcept>
namespace fs=std::filesystem;
template<class T> T* get(TFile& f,const std::string& n) {auto* p=dynamic_cast<T*>(f.Get(n.c_str()));if(!p)throw std::runtime_error("Missing "+n);return p;}
template<class T> std::unique_ptr<T> hist(TFile& f,const std::string& n) {auto* p=get<T>(f,n);p->SetDirectory(nullptr);return std::unique_ptr<T>(p);}
int main(int argc,char** argv) {
  try {
    if(argc!=3)throw std::runtime_error("review INPUT_DIR OUTPUT_DIR");
    TH1::AddDirectory(false);fs::create_directories(argv[2]);
    std::ofstream points(fs::path(argv[2])/"profile_points.tsv"), meta(fs::path(argv[2])/"samples.tsv"), summary(fs::path(argv[2])/"channel_summary.tsv");
    points<<std::setprecision(10);summary<<std::setprecision(10);
    points<<"energy\tsample\tcell\tplot\tx\txlow\txhigh\tn\tmean\tsem\n";
    meta<<"energy\tsample\tevents\tchannels\tpedestal\n";
    summary<<"energy\tsample\tcell\thits\thits_per_event\tmean_hg\tmean_lg\tmean_energy\tmean_hg_energy\tmean_lg_energy\tped_hg\tped_lg\n";
    for(const std::string energy:{"0.5","1","5","10","30","100","250"}) {
      for(const std::string sample:{"data","mc"}) {
        fs::path path=fs::path(argv[1])/(sample=="data"?"data":"simulation/threshold")/(energy+"GeV.root");
        TFile file(path.c_str());if(file.IsZombie())throw std::runtime_error("Cannot open "+path.string());
        auto events=get<TParameter<Long64_t>>(file,"processed_events")->GetVal();
        int channels=0;TIter next(file.GetListOfKeys());while(auto* k=dynamic_cast<TKey*>(next()))if(std::string(k->GetName()).rfind("channel_",0)==0)++channels;
        meta<<energy<<'\t'<<sample<<'\t'<<events<<'\t'<<channels<<'\t'<<get<TNamed>(file,"pedestal_file")->GetTitle()<<'\n';
        for(int cell:{20008,120012,220008,320012,420008,520012,620008,720012,820008,920012,1020008,1120012,1220008,1420010,1820008,2220008,2820008}) {
          const auto prefix="channel_"+std::to_string(cell)+"/";
          if(!file.GetDirectory(("channel_"+std::to_string(cell)).c_str()))continue;
          auto density=hist<TH2D>(file,prefix+"hg_lg");auto eh=hist<TH2D>(file,prefix+"hg_energy");auto el=hist<TH2D>(file,prefix+"lg_energy");auto en=hist<TH2D>(file,prefix+"hg_selected_energy");
          summary<<energy<<'\t'<<sample<<'\t'<<cell<<'\t'<<density->GetEntries()<<'\t'<<density->GetEntries()/events<<'\t'<<density->GetMean(1)<<'\t'<<density->GetMean(2)<<'\t'<<en->GetMean(2)<<'\t'<<eh->GetMean(2)<<'\t'<<el->GetMean(2)<<'\t'<<get<TParameter<double>>(file,prefix+"pedestal_hg")->GetVal()<<'\t'<<get<TParameter<double>>(file,prefix+"pedestal_lg")->GetVal()<<'\n';
          for(const std::string name:{"hg_lg_zoom","hg_energy","lg_energy","hg_selected_energy","lg_selected_energy"}) {
            auto p=hist<TProfile>(file,prefix+name+"_profile");
            for(int b=1;b<=p->GetNbinsX();++b)if(p->GetBinEntries(b)>0)points<<energy<<'\t'<<sample<<'\t'<<cell<<'\t'<<name<<'\t'<<p->GetBinCenter(b)<<'\t'<<p->GetXaxis()->GetBinLowEdge(b)<<'\t'<<p->GetXaxis()->GetBinUpEdge(b)<<'\t'<<p->GetBinEntries(b)<<'\t'<<p->GetBinContent(b)<<'\t'<<p->GetBinError(b)<<'\n';
          }
        }
        std::cout<<energy<<" "<<sample<<" events="<<events<<std::endl;
      }
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
