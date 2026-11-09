#include "StudyCommon.hh"
#include <array>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
struct Bin {long long truthN=0,detected=0,hgN=0;double truthE=0,detectedE=0,hg=0,hgPred=0,lg=0,pred=0,calPred=0;};
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"manifest",""},{"audit",""},{"hl",""},{"mip",""},{"pedestal",""},{"out",""}});
 auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);std::map<int,double>digi;
 std::ifstream audit(o["audit"]);require(bool(audit),"Missing audit");std::string line;std::getline(audit,line);
 while(std::getline(audit,line)){std::istringstream row(line);int id,g,nd;double mp,chi,cal,d;require(bool(row>>id>>g>>mp>>chi>>nd>>cal>>d),"Invalid audit row");digi[id]=d;}
 std::array<double,12> edges={0,.05,.1,.2,.305,.5,1,2,5,10,30,100};
 std::map<std::pair<int,int>,Bin> bins;Stats posx,posy,truthSum,rawHits;double maxPrimaryDifference=0;
 std::ifstream manifest(o["manifest"]);require(bool(manifest),"Missing manifest");fs::path out=o["out"];fs::create_directories(out);
 std::ofstream files(out/"files.tsv");files<<"decode\ttruth\tevents\n";
 while(std::getline(manifest,line)){
  auto path=line.substr(0,line.find('\t'));auto truth=path;const std::string from="/decode/e-/e_ssa/threshold/";auto at=truth.find(from);require(at!=std::string::npos,"Expected e_ssa/threshold");truth.replace(at,from.size(),"/generate/e-/e_ssa/");
  auto rf=openFile(path),tf=openFile(truth);auto*r=object<TTree>(*rf,"Raw_Hit"),*t=object<TTree>(*tf,"MC_Truth");require(r->GetEntries()==t->GetEntries(),"Truth/decode entry count differs");
  r->SetBranchStatus("*",0);t->SetBranchStatus("*",0);
  std::vector<int>*rid=nullptr,*tid=nullptr,*tags=nullptr;std::vector<double>*h=nullptr,*l=nullptr,*e=nullptr;double x,y,primary;int event;
  bind(r,"CellID",&rid);bind(r,"HitTag",&tags);bind(r,"HG_Charge",&h);bind(r,"LG_Charge",&l);bind(r,"Event_Time",&event);
  bind(t,"CellID",&tid);bind(t,"Hit_Energy",&e);bind(t,"PrimaryPosX",&x);bind(t,"PrimaryPosY",&y);bind(t,"PrimaryEnergy",&primary);
  for(Long64_t i=0;i<t->GetEntries();++i){require(t->GetEntry(i)>0&&r->GetEntry(i)>0,"Read failed");require(event==i,"Unexpected MC event identity");require(tid->size()==e->size()&&rid->size()==h->size()&&rid->size()==l->size()&&rid->size()==tags->size(),"Vector size mismatch");
   std::map<int,size_t> indices;for(size_t j=0;j<rid->size();++j)if(tags->at(j)==1)require(indices.emplace(rid->at(j),j).second,"Duplicate raw ID");
   std::set<int> seen;double sum=0;for(size_t j=0;j<tid->size();++j){const int full=tid->at(j),id=physical(full);require(seen.insert(full).second,"Duplicate truth ID");if(!connected(id))continue;double en=e->at(j);sum+=en;auto&p=pars.at(id);if(p.mip<=0)continue;
    int band=std::upper_bound(edges.begin(),edges.end(),en)-edges.begin()-1;if(band<0)continue;int category=group(id)*2+(std::abs(digi.at(id)-p.mip)>1e-8);auto&b=bins[{category,band}];++b.truthN;b.truthE+=en;
    if(!indices.count(full))continue;
    size_t k=indices.at(full);++b.detected;b.detectedE+=en;double prediction=en/.305*digi.at(id);b.pred+=prediction;b.calPred+=en/.305*p.mip;
    // Undo the raw LG intercept before testing the digitized signal scale.
    b.lg+=(l->at(k)-p.pl)*p.gain+p.intercept;
    if(h->at(k)<2400){++b.hgN;b.hg+=h->at(k)-p.ph;b.hgPred+=prediction;}
   }
   for(const auto& kv:indices)require(seen.count(kv.first),"Decoded ID absent from truth");
   posx.add(x);posy.add(y);truthSum.add(sum/1000);rawHits.add(indices.size());maxPrimaryDifference=std::max(maxPrimaryDifference,std::abs(primary/1000-std::stod(fs::path(path).parent_path().filename().string())));
  }
  files<<path<<'\t'<<truth<<'\t'<<t->GetEntries()<<'\n';r->ResetBranchAddresses();t->ResetBranchAddresses();delete rid;delete tid;delete tags;delete h;delete l;delete e;
 }
 std::ofstream table(out/"response.tsv");table<<std::setprecision(12)<<"category\ttruth_min_MeV\ttruth_max_MeV\ttruth_hits\tdecoded_hits\tefficiency\tHG_hits\tHG_over_digi_prediction\tLG_corrected_over_digi_prediction\tLG_corrected_over_calib_prediction\n";
 for(auto&[key,b]:bins){auto[category,band]=key;table<<category<<'\t'<<edges[band]<<'\t'<<(band+1<int(edges.size())?edges[band+1]:1e9)<<'\t'<<b.truthN<<'\t'<<b.detected<<'\t'<<double(b.detected)/b.truthN<<'\t'<<b.hgN<<'\t'<<(b.hgPred?b.hg/b.hgPred:0)<<'\t'<<(b.pred?b.lg/b.pred:0)<<'\t'<<(b.calPred?b.lg/b.calPred:0)<<'\n';}
 std::ofstream summary(out/"summary.tsv");summary<<std::setprecision(12)<<"events\tprimary_x_mean_mm\tprimary_x_RMS_mm\tprimary_y_mean_mm\tprimary_y_RMS_mm\ttruth_energy_mean_GeV\ttruth_energy_RMS_GeV\tdecoded_hits_mean\tmax_primary_energy_difference_GeV\n"<<posx.n<<'\t'<<posx.mean<<'\t'<<posx.rms()<<'\t'<<posy.mean<<'\t'<<posy.rms()<<'\t'<<truthSum.mean<<'\t'<<truthSum.rms()<<'\t'<<rawHits.mean<<'\t'<<maxPrimaryDifference<<'\n';
 std::cout<<"Truth/decode matched "<<posx.n<<" events"<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
