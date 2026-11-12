#pragma once
// All bins, including flow bins, are retained. The compact table stores nonzero
// bin indices/counts, the original axes, entries and ROOT statistics for each
// distribution; absent indices have exactly zero content.
class Histograms {
 TFile&file;std::unique_ptr<TDirectory>dir;std::unique_ptr<TTree>rows;
 size_t n=0;std::string prefix,name;int nbins=0;double xmin=0,xmax=0,entries=0,stats[4]={};
 std::vector<int>indices;std::vector<double>counts;bool compact;
 public:explicit Histograms(TFile&f,const std::string&mode="compact"):file(f),compact(mode=="compact"){
  ts::require(mode=="compact"||mode=="legacy","Histogram storage must be compact or legacy");
  if(compact){file.cd();rows=std::make_unique<TTree>("histograms","Exact absolute distributions; sparse nonzero bins");rows->SetAutoFlush(-16*1024*1024);rows->SetMaxVirtualSize(16*1024*1024);rows->Branch("name",&name);rows->Branch("nbins",&nbins);rows->Branch("xmin",&xmin);rows->Branch("xmax",&xmax);rows->Branch("entries",&entries);rows->Branch("stats",stats,"stats[4]/D");rows->Branch("bin_indices",&indices);rows->Branch("bin_counts",&counts);}
 }
 std::string write(TH1D&h){
  if(compact){name=h.GetName();nbins=h.GetNbinsX();xmin=h.GetXaxis()->GetXmin();xmax=h.GetXaxis()->GetXmax();entries=h.GetEntries();h.GetStats(stats);indices.clear();counts.clear();for(int i=0;i<=nbins+1;++i){double v=h.GetBinContent(i);if(v!=0){indices.push_back(i);counts.push_back(v);}}ts::require(rows->Fill()>0,"Compact histogram write failed");return "histograms/entry_"+std::to_string(n++);}
  if(n%1024==0){close();prefix="hist_"+std::to_string(n/1024);file.cd();dir.reset(file.mkdir(prefix.c_str()));ts::require(bool(dir),"Cannot create histogram directory");}dir->cd();ts::require(h.Write()>0,"Histogram write failed");++n;return prefix+"/"+h.GetName();
 }
 void close(){if(rows){file.cd();ts::require(rows->Write()>0,"Compact histogram table write failed");rows->SetDirectory(nullptr);rows.reset();}if(dir){dir->SaveSelf(true);file.cd();dir.reset();}}
 ~Histograms(){if(rows)rows->SetDirectory(nullptr);rows.reset();dir.reset();}
};
namespace ts {
// Reconstruct either layout. Plotting itself continues to use numerical TSVs.
inline std::unique_ptr<TH1D> readHistogram(TFile&f,const std::string&reference){
 const std::string prefix="histograms/entry_";
 if(reference.compare(0,prefix.size(),prefix)!=0){auto*h=dynamic_cast<TH1D*>(f.Get(reference.c_str()));require(h,"Missing histogram "+reference);auto result=std::unique_ptr<TH1D>(static_cast<TH1D*>(h->Clone()));result->SetDirectory(nullptr);return result;}
 auto*t=dynamic_cast<TTree*>(f.Get("histograms"));require(t,"Missing compact histogram table");Long64_t entry=std::stoll(reference.substr(prefix.size()));require(entry>=0&&entry<t->GetEntries(),"Invalid histogram table entry");
 std::string*name=nullptr;std::vector<int>*indices=nullptr;std::vector<double>*counts=nullptr;int nbins;double xmin,xmax,entries,stats[4];
 t->SetBranchAddress("name",&name);t->SetBranchAddress("nbins",&nbins);t->SetBranchAddress("xmin",&xmin);t->SetBranchAddress("xmax",&xmax);t->SetBranchAddress("entries",&entries);t->SetBranchAddress("stats",stats);t->SetBranchAddress("bin_indices",&indices);t->SetBranchAddress("bin_counts",&counts);require(t->GetEntry(entry)>0,"Cannot read compact histogram entry");
 auto result=std::make_unique<TH1D>(name->c_str(),"Absolute saved distribution",nbins,xmin,xmax);result->SetDirectory(nullptr);require(indices->size()==counts->size(),"Invalid sparse histogram");for(size_t i=0;i<indices->size();++i){require((*indices)[i]>=0&&(*indices)[i]<=nbins+1,"Invalid saved bin index");result->SetBinContent((*indices)[i],(*counts)[i]);}result->PutStats(stats);result->SetEntries(entries);t->ResetBranchAddresses();delete name;delete indices;delete counts;return result;
}
}
