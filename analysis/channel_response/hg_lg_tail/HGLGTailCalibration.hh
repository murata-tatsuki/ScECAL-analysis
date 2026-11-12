#pragma once
#include "StudyCommon.hh"
#include <array>
#include <set>

// Match the current Calibration.cxx factor extraction. Keep this local to tail
// analysis: older impact studies intentionally reproduce the legacy calibration.
inline std::map<int,Parameter> tailParameters(const std::string& hl,
                                             const std::string& mip,
                                             const std::string& ped) {
 std::map<int,Parameter> a;
 for(int l=0;l<30;++l)for(int c=0;c<6;++c)for(int ch=0;ch<(c==5?30:36);++ch)
  a[l*100000+c*10000+ch]=Parameter{};
 std::set<int> validHL, validMIP;
 {
  auto f=openFile(hl);auto*t=object<TTree>(*f,"InterCalib");
  int id=0;double slope=0,intercept=0;
  bind(t,"CellID",&id);bind(t,"Slope",&slope);bind(t,"Intercept",&intercept);
  for(Long64_t i=0;i<t->GetEntries();++i){
   require(t->GetEntry(i)>0,"Cannot read HL factor entry");
   if(id<0||!connected(physical(id))||!std::isfinite(slope)||slope<.02||slope>.05||!std::isfinite(intercept))continue;
   double gain=1/slope,offset=-intercept/slope;
   if(!std::isfinite(gain)||!std::isfinite(offset))continue;
   int cell=physical(id);a.at(cell).gain=gain;a.at(cell).intercept=offset;validHL.insert(cell);
  }
 }
 {
  auto f=openFile(mip);auto*t=object<TTree>(*f,"MIP_Fit");
  int id=0,ndf=0;double mpv=0,chi=0;
  bind(t,"CellID",&id);bind(t,"LandauMPV",&mpv);bind(t,"ChiSquare",&chi);bind(t,"NDF",&ndf);
  for(Long64_t i=0;i<t->GetEntries();++i){
   require(t->GetEntry(i)>0,"Cannot read MIP factor entry");
   if(id<0||!connected(physical(id))||!std::isfinite(mpv)||mpv<=0||ndf<=0||!std::isfinite(chi)||chi<0||chi/ndf>(group(id)?2.:1.7))continue;
   int cell=physical(id);a.at(cell).mip=mpv;validMIP.insert(cell);
  }
 }
 // Count each accepted physical channel once, and exclude layers 30/31.
 std::array<long double,2> gains{},offsets{},mips{};
 std::array<int,2> nh{},nm{};
 for(const auto& [id,p]:a){int g=group(id);
  if(validHL.count(id)){gains[g]+=p.gain;offsets[g]+=p.intercept;++nh[g];}
  if(validMIP.count(id)){mips[g]+=p.mip;++nm[g];}
 }
 for(int g=0;g<2;++g)require(nh[g]>0&&nm[g]>0,"No accepted constants in a SiPM group");
 for(auto& [id,p]:a){int g=group(id);
  p.hlFallback=!validHL.count(id);p.mipFallback=!validMIP.count(id);
  if(p.hlFallback){p.gain=static_cast<double>(gains[g]/nh[g]);p.intercept=static_cast<double>(offsets[g]/nh[g]);}
  if(p.mipFallback)p.mip=static_cast<double>(mips[g]/nm[g]);
 }
 {
  auto f=openFile(ped);auto*t=object<TTree>(*f,"ChnLevel");
  std::vector<int>*ids=nullptr;std::vector<double>*ph=nullptr,*pl=nullptr;
  bind(t,"CellID",&ids);bind(t,"PedHighMean",&ph);bind(t,"PedLowMean",&pl);
  for(Long64_t i=0;i<t->GetEntries();++i){
   require(t->GetEntry(i)>0,"Cannot read pedestal entry");
   require(ids->size()==ph->size()&&ids->size()==pl->size(),"Pedestal size mismatch");
   for(size_t j=0;j<ids->size();++j){int id=physical(ids->at(j));if(!a.count(id))continue;
    auto&p=a.at(id);p.ph=static_cast<float>(ph->at(j));p.pl=static_cast<float>(pl->at(j));p.ped=true;
   }
  }
  t->ResetBranchAddresses();delete ids;delete ph;delete pl;
 }
 return a;
}
