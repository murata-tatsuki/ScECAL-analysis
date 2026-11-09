#include "EventHistograms.hh"
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <set>
#include <regex>
std::set<int> loadMask(const std::string&path){std::ifstream f(path);std::string txt((std::istreambuf_iterator<char>(f)),{});auto start=txt.find("const unordered_set<int> bad_channels = {");require(start!=std::string::npos,"Mask not found");start=txt.find('{',start);txt=txt.substr(start,txt.find("};",start)-start);std::regex number("[0-9]+");std::set<int> result;for(std::sregex_iterator it(txt.begin(),txt.end(),number),end;it!=end;++it)result.insert(std::stoi(it->str()));return result;}
struct Group {Stats n,active,energy,x0,y0,x9,y9,x10,y10;long long empty=0,track=0,cog=0;};
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"input",""},{"out",""},{"comparison-source",""}});fs::path input=o["input"],out=o["out"];fs::create_directories(out/"figures");auto bad=loadMask(o["comparison-source"]);auto source=openFile((input/"study.root").string());auto*t=object<TTree>(*source,"paired_events");double en;t->SetBranchStatus("*",0);bind(t,"nominal_GeV",&en);std::vector<double> mc;for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);mc.push_back(en);}t->ResetBranchAddresses();std::sort(mc.begin(),mc.end());double ref=eventplots::quantile(mc,.5);gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);std::array<std::array<Group,3>,2> group{};std::array<std::unique_ptr<TH2D>,2> xy;
 for(int s=0;s<2;++s){xy[s]=std::make_unique<TH2D>((s?"mc_cog":"data_cog"),";Layer 9 CoG X [mm];Layer 10 CoG Y [mm]",90,-135,135,90,-135,135);std::ifstream mf(input/"inputs"/(s?"mc.tsv":"data.tsv"));std::string line;std::vector<std::string> paths;while(std::getline(mf,line))paths.push_back(line.substr(line.find('\t')+1));long long perFile=s?std::max<size_t>(200,20000/paths.size()):2000;
  for(const auto&path:paths){auto file=openFile(path);auto*tree=object<TTree>(*file,"Calib_Hit");tree->SetBranchStatus("*",0);std::vector<int>*ids=nullptr;std::vector<double>*energy=nullptr,*x=nullptr,*y=nullptr;bind(tree,"CellID",&ids);bind(tree,"Hit_Energy",&energy);bind(tree,"Hit_X",&x);bind(tree,"Hit_Y",&y);tree->SetCacheSize(16*1024*1024);for(auto key:{"CellID","Hit_Energy","Hit_X","Hit_Y"})tree->AddBranchToCache(key,true);tree->StopCacheLearningPhase();Long64_t count=tree->GetEntries(),take=std::min(count,perFile);int blocks=std::min<Long64_t>(10,take);
   for(int b=0;b<blocks;++b){Long64_t begin=count*b/blocks,end=count*(b+1)/blocks,n=take/blocks+(b<take%blocks);begin+=(end-begin-n)/2;
    for(Long64_t i=begin;i<begin+n;++i){tree->GetEntry(i);std::array<double,30> sums{},xx{},yy{};int hits=0,positiveHits=0;double signedE=0,posE=0;for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!connected(id))continue;++hits;double e=energy->at(j)/1000;signedE+=e;if(e<0)continue;posE+=e;if(bad.count(id))continue;++positiveHits;int l=id/100000;sums[l]+=e;xx[l]+=x->at(j)*e;yy[l]+=y->at(j)*e;}
     bool shower=posE>.5*ref&&posE<1.5*ref,low=signedE<.1*ref;int active=0;for(double e:sums)active+=e>0;bool cog=true;for(int l:{9,10})cog&=sums[l]>0&&std::abs(xx[l]/sums[l])<20&&std::abs(yy[l]/sums[l])<20;
     for(int k=0;k<3;++k)if(k==0||(k==1&&low)||(k==2&&shower)){auto&g=group[s][k];g.n.add(hits);g.energy.add(signedE);g.active.add(active);g.empty+=hits==0;g.track+=active>=20&&positiveHits<=2*active;g.cog+=cog;for(int l:{0,9,10})if(sums[l]>0){auto&sx=l==0?g.x0:(l==9?g.x9:g.x10);auto&sy=l==0?g.y0:(l==9?g.y9:g.y10);sx.add(xx[l]/sums[l]);sy.add(yy[l]/sums[l]);}}
     if(shower&&sums[9]>0&&sums[10]>0)xy[s]->Fill(xx[9]/sums[9],yy[10]/sums[10]);
    }
   }tree->ResetBranchAddresses();delete ids;delete energy;delete x;delete y;
  }std::cout<<input.filename()<<" topology "<<(s?"mc":"data")<<" sampled "<<group[s][0].n.n<<std::endl;
 }
 std::ofstream table(out/"topology.tsv");table<<std::setprecision(12)<<"sample\tselection\tevents\tmean_hits\tmean_active_layers\tmean_signed_GeV\tempty\ttracklike\tcog20\tcog0_n\tcog0_x\tcog0_y\tcog9_n\tcog9_x\tcog9_y\tcog10_n\tcog10_x\tcog10_y\n";
 for(int s=0;s<2;++s)for(int k=0;k<3;++k){auto&g=group[s][k];table<<(s?"mc":"data")<<'\t'<<std::array<const char*,3>{"all","low_signed","shower"}[k]<<'\t'<<g.n.n<<'\t'<<g.n.mean<<'\t'<<g.active.mean<<'\t'<<g.energy.mean<<'\t'<<g.empty<<'\t'<<g.track<<'\t'<<g.cog<<'\t'<<g.x0.n<<'\t'<<g.x0.mean<<'\t'<<g.y0.mean<<'\t'<<g.x9.n<<'\t'<<g.x9.mean<<'\t'<<g.y9.mean<<'\t'<<g.x10.n<<'\t'<<g.x10.mean<<'\t'<<g.y10.mean<<'\n';}
 TFile f((out/"topology.root").c_str(),"RECREATE");TCanvas c("position","",1300,600);c.Divide(2,1);for(int s=0;s<2;++s){xy[s]->Write();c.cd(s+1);gPad->SetRightMargin(.16);xy[s]->SetTitle((input.filename().string()+": "+(s?"MC":"data")+", shower selection").c_str());xy[s]->Draw("COLZ");}c.SaveAs((out/"figures/shower_position.png").c_str());return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
