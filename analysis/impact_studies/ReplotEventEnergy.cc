// Redraw completed studies without repeating raw/calib matching or channel plots.
#include "EventHistograms.hh"
#include <TParameter.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>

int main(int argc,char**argv)try {
 auto o=options(argc,argv,{{"input",""},{"out",""}});fs::path input=o["input"],out=o["out"];
 require(!input.empty()&&!out.empty(),"--input and --out are required");
 require(fs::weakly_canonical(input)!=fs::weakly_canonical(out),"Use a separate staging output directory");
 auto source=openFile((input/"study.root").string());
 double beam=object<TParameter<double>>(*source,"beam_energy_GeV")->GetVal();
 std::ifstream summary(input/"summary.tsv");require(bool(summary),"Missing summary.tsv");std::string line;std::getline(summary,line);
 std::array<std::vector<std::string>,3> rows;
 for(int s=0;s<3;++s){require(bool(std::getline(summary,line)),"Incomplete summary");std::istringstream row(line);std::string value;while(std::getline(row,value,'\t'))rows[s].push_back(value);require(rows[s].size()==8&&rows[s][0]==eventplots::samples[s],"Unexpected summary format");}
 eventplots::Values values;
 auto*paired=object<TTree>(*source,"paired_events");double nominal,zero;paired->SetBranchStatus("*",0);bind(paired,"nominal_GeV",&nominal);bind(paired,"intercept_zero_GeV",&zero);
 for(Long64_t i=0;i<paired->GetEntries();++i){require(paired->GetEntry(i)>0,"Cannot read paired event");values[1].push_back(nominal);values[2].push_back(zero);}paired->ResetBranchAddresses();
 auto*data=dynamic_cast<TTree*>(source->Get("data_events"));std::unique_ptr<TFile> cache;
 if(!data&&fs::exists(input/"event_histograms.root")){cache=openFile((input/"event_histograms.root").string());data=dynamic_cast<TTree*>(cache->Get("data_events"));}
 if(data){double energy;data->SetBranchStatus("*",0);bind(data,"energy_GeV",&energy);for(Long64_t i=0;i<data->GetEntries();++i){require(data->GetEntry(i)>0,"Cannot read data event");values[0].push_back(energy);}data->ResetBranchAddresses();}
 else {
  // Old study.root stores MC event totals only. Read just two data branches,
  // preserving the manifest order, connected-channel selection and event limit.
  const auto count=std::stoll(rows[0][1]);std::ifstream manifest(input/"inputs/data.tsv");require(bool(manifest),"Missing data manifest");
  while(static_cast<long long>(values[0].size())<count&&std::getline(manifest,line)){
   auto tab=line.find('\t');require(tab!=std::string::npos,"Invalid data manifest");auto file=openFile(line.substr(tab+1));auto*t=object<TTree>(*file,"Calib_Hit");t->SetBranchStatus("*",0);
   std::vector<int>*ids=nullptr;std::vector<double>*energies=nullptr;bind(t,"CellID",&ids);bind(t,"Hit_Energy",&energies);t->SetCacheSize(32*1024*1024);t->AddBranchToCache("CellID",true);t->AddBranchToCache("Hit_Energy",true);t->StopCacheLearningPhase();
   for(Long64_t i=0;i<t->GetEntries()&&static_cast<long long>(values[0].size())<count;++i){require(t->GetEntry(i)>0,"Cannot read calibrated data");require(ids->size()==energies->size(),"Data hit vector mismatch");double energy=0;for(size_t j=0;j<ids->size();++j)if(connected(physical(ids->at(j))))energy+=energies->at(j)/1000.;values[0].push_back(energy);}
   t->ResetBranchAddresses();delete ids;delete energies;std::cout<<"Data totals: "<<values[0].size()<<" / "<<count<<std::endl;
  }
 }
 for(int s=0;s<3;++s){Stats stats;for(double x:values[s])stats.add(x);require(stats.n==std::stoll(rows[s][1]),"Event count differs from completed study");for(auto pair:{std::pair<double,double>{stats.mean,std::stod(rows[s][2])},{stats.rms(),std::stod(rows[s][4])}})require(std::abs(pair.first-pair.second)<1e-9*std::max(1.,std::abs(pair.second)),"Event statistics differ from completed study");}
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetTextFont(42);gStyle->SetTitleFont(42,"XYZ");gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleSize(.04,"XYZ");gStyle->SetLabelSize(.035,"XYZ");
 fs::create_directories(out);eventplots::Histograms histograms(values);
 {TFile file((out/"event_histograms.root").c_str(),"RECREATE");require(!file.IsZombie(),"Cannot create histograms");histograms.write();TParameter<double>("beam_energy_GeV",beam).Write();TNamed("source_study",fs::canonical(input/"study.root").c_str()).Write();TTree events("data_events","Data event totals with the original study selection");double energy;events.Branch("energy_GeV",&energy);for(double x:values[0]){energy=x;events.Fill();}events.Write();}
 histograms.metadata(out);histograms.draw(out,beam);
 std::ofstream updated(out/"summary.tsv");updated<<"sample\tevents\tmean_GeV\tsem_GeV\trms_GeV\trms_over_mean\tunderflow\toverflow\n";
 for(int s=0;s<3;++s){for(int col=0;col<6;++col)updated<<rows[s][col]<<'\t';auto&h=*histograms.totals[s];updated<<h.GetBinContent(0)<<'\t'<<h.GetBinContent(h.GetNbinsX()+1)<<'\n';}
 std::cout<<"Replotted "<<beam<<" GeV: "<<out<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
