// Small MC reference for the shared diagnostics; no channel histograms are needed.
#include "StudyCommon.hh"
#include <TMemFile.h>
#include <TParameter.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"mc-manifest",""},{"hl",""},{"mip",""},{"pedestal",""},{"out",""},{"energy","100"}});auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);fs::path out=o["out"];fs::create_directories(out);
 TMemFile result("study.root","RECREATE");TTree events("paired_events","Current MC calibrated event sums; intercept-zero branch only for compatibility with legacy diagnostic");int fileIndex=-1;Long64_t entry;double nominal,zero;
 events.Branch("file_index",&fileIndex);events.Branch("entry",&entry);events.Branch("nominal_GeV",&nominal);events.Branch("intercept_zero_GeV",&zero);
 std::ifstream manifest(o["mc-manifest"]);require(bool(manifest),"Missing manifest");std::string line;Stats stats;std::vector<double> sums;
 while(std::getline(manifest,line)){++fileIndex;auto f=openFile(line.substr(line.find('\t')+1));auto*t=object<TTree>(*f,"Calib_Hit");t->SetBranchStatus("*",0);std::vector<int>*id=nullptr;std::vector<double>*e=nullptr,*h=nullptr,*l=nullptr;bind(t,"CellID",&id);bind(t,"Hit_Energy",&e);bind(t,"Hit_HG_Energy",&h);bind(t,"Hit_LG_Energy",&l);
  for(entry=0;entry<t->GetEntries();++entry){require(t->GetEntry(entry)>0,"Read failure");require(id->size()==e->size()&&id->size()==h->size()&&id->size()==l->size(),"Vector sizes differ");nominal=zero=0;for(size_t j=0;j<id->size();++j){int cell=physical(id->at(j));if(!connected(cell))continue;double v=e->at(j);require(std::isfinite(v),"Nonfinite energy");bool lg=v==l->at(j)&&v!=h->at(j);nominal+=v/1000;zero+=(v-(lg?pars.at(cell).offset():0))/1000;}stats.add(nominal);sums.push_back(nominal);events.Fill();}
  t->ResetBranchAddresses();delete id;delete e;delete h;delete l;
 }
 require(stats.n>0,"No MC events");result.cd();events.Write();TParameter<double>("beam_energy_GeV",std::stod(o["energy"])).Write();result.Write();Long64_t size=result.GetEND();std::vector<char> bytes(size);require(result.CopyTo(bytes.data(),size)==size,"Incomplete in-memory ROOT copy");std::ofstream output(out/"study.root.buffer",std::ios::binary);output.write(bytes.data(),size);output.close();require(bool(output),"ROOT write failed");fs::rename(out/"study.root.buffer",out/"study.root");events.SetDirectory(nullptr);result.Close();
 std::sort(sums.begin(),sums.end());std::ofstream summary(out/"reference.tsv");summary<<std::setprecision(12)<<"mc_events\tsigned_mean_GeV\tsigned_median_GeV\n"<<stats.n<<'\t'<<stats.mean<<'\t'<<sums[sums.size()/2]<<'\n';std::cout<<"MC reference: "<<stats.n<<" events, median="<<sums[sums.size()/2]<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
