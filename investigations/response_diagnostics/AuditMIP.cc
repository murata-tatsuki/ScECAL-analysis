#include "StudyCommon.hh"
#include <iomanip>
#include <iostream>
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"hl",""},{"mip",""},{"pedestal",""},{"out",""}});auto cal=parameters(o["hl"],o["mip"],o["pedestal"]);std::map<int,double> digi,mpv,chi;std::map<int,int> ndf;for(auto&[id,p]:cal)digi[id]=1;
 double sums[2]={};long long n[2]={};auto f=openFile(o["mip"]);auto*t=object<TTree>(*f,"MIP_Fit");int id,nd;double m,c;bind(t,"CellID",&id);bind(t,"LandauMPV",&m);bind(t,"ChiSquare",&c);bind(t,"NDF",&nd);
 for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id<0||id/100000>=30)continue;mpv[id]=m;chi[id]=c;ndf[id]=nd;if(m<0)continue;sums[group(id)]+=m;++n[group(id)];if(digi.count(id))digi[id]=m;}
 std::ofstream out(o["out"]);out<<std::setprecision(12)<<"cell\tgroup10um\traw_mpv\tchi2\tndf\tcalibration_mip\tdigitization_mip\tdigi_over_calib\tcalib_fallback\tdigi_fallback\n";
 long long negative=0,different=0,fallback=0;
 for(auto&[id,p]:cal){double d=digi[id];bool missing=d==1;if(missing)d=sums[group(id)]/n[group(id)];if(d<=5)d=group(id)?90:400;negative+=p.mip<0;different+=std::abs(d-p.mip)>1e-8;fallback+=p.mipFallback;out<<id<<'\t'<<group(id)<<'\t'<<(mpv.count(id)?mpv[id]:0)<<'\t'<<chi[id]<<'\t'<<ndf[id]<<'\t'<<p.mip<<'\t'<<d<<'\t'<<d/p.mip<<'\t'<<p.mipFallback<<'\t'<<missing<<'\n';}
 std::cout<<"channels="<<cal.size()<<" negative_calibration="<<negative<<" different_digi_vs_calib="<<different<<" calib_fallback="<<fallback<<" digi_mean_outer="<<sums[0]/n[0]<<" digi_mean_inner="<<sums[1]/n[1]<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
