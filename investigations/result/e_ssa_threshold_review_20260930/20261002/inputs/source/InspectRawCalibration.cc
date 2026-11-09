#include "StudyCommon.hh"
#include <TROOT.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <array>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <regex>
#include <set>
#include <unordered_map>

struct Line {
 long long n=0; double x=0,y=0,xx=0,xy=0,yy=0;
 void add(double a,double b){++n;x+=a;y+=b;xx+=a*a;xy+=a*b;yy+=b*b;}
 double slope()const{return n>1 && xx-x*x/n>0?(xy-x*y/n)/(xx-x*x/n):0;}
 double intercept()const{return n?y/n-slope()*x/n:0;}
 double width()const{return n>2?std::sqrt(std::max(0.,(yy-y*y/n-slope()*(xy-x*y/n))/(n-2))):0;}
 double se()const{return n>2 && xx-x*x/n>0?width()/std::sqrt(xx-x*x/n):0;}
};
struct Extra {double sh=0,sl=0,slope=0,b=0,xmax=0,digiA=0,sw=0,threshold=0,thresholdSigma=0;bool hl=false;};
struct Channel {
 Stats temperature,hgError,lgError,mipInferred,offH,offL,offCoreH,offCoreL;
 Line overlap,low,high,offLine;
 std::array<Stats,5> residual,energy;
 long long raw=0,off=0,hits=0,lg=0,negative=0,clipH=0,clipL=0,overL=0,belowSwitch=0,lowHG=0,gain0=0;
 double e=0,lgE=0,clipLE=0,tempAfterDelta=0;
};
struct EventSum {Stats energy,hits,lgEnergy,lgHits,thermalOrder,clipEnergy,negative;};
int band(double h){return h<800?0:h<1400?1:h<2200?2:h<2600?3:4;}
std::set<int> mask(const std::string&path){std::ifstream f(path);std::string t((std::istreambuf_iterator<char>(f)),{});auto p=t.find("const unordered_set<int> bad_channels = {");require(p!=std::string::npos,"Missing mask");p=t.find('{',p);t=t.substr(p,t.find("};",p)-p);std::regex re("[0-9]+");std::set<int> v;for(std::sregex_iterator it(t.begin(),t.end(),re),end;it!=end;++it)v.insert(std::stoi(it->str()));return v;}
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"input",""},{"out",""},{"hl",""},{"mip",""},{"pedestal",""},{"threshold",""},{"comparison-source",""},{"data-per-file","2000"},{"mc-per-file","1000"},{"mc-files","6"}});
 fs::path out=o["out"];require(!fs::exists(out/"channels.tsv"),"Results already exist; choose another output directory");fs::create_directories(out/"figures");fs::create_directories(out/"inputs");
 auto pars=parameters(o["hl"],o["mip"],o["pedestal"]);auto bad=mask(o["comparison-source"]);std::map<int,Extra> extra;for(auto&[id,p]:pars)extra[id]=Extra{};
 {auto f=openFile(o["pedestal"]);auto*t=object<TTree>(*f,"ChnLevel");std::vector<int>*ids=nullptr;std::vector<double>*sh=nullptr,*sl=nullptr;bind(t,"CellID",&ids);bind(t,"PedHighSig",&sh);bind(t,"PedLowSig",&sl);for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);for(size_t j=0;j<ids->size();++j)if(extra.count(physical(ids->at(j)))){auto&p=extra.at(physical(ids->at(j)));p.sh=sh->at(j);p.sl=sl->at(j);}}t->ResetBranchAddresses();delete ids;delete sh;delete sl;}
 {auto f=openFile(o["hl"]);auto*t=object<TTree>(*f,"InterCalib");int id;double a,b,x;bind(t,"CellID",&id);bind(t,"Slope",&a);bind(t,"Intercept",&b);bind(t,"XMax",&x);double sum=0,minsw=4000;int n=0;for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id<0||id/100000>=30)continue;sum+=1/a;++n;minsw=std::min(minsw,x-600);if(extra.count(id)){auto&p=extra.at(id);p.hl=true;p.slope=a;p.b=b;p.xmax=x;p.digiA=1/a;p.sw=x-600;}}for(auto&[id,p]:extra)if(!p.hl){p.digiA=sum/n;p.sw=minsw;}}
 {auto f=openFile(o["threshold"]);auto*t=object<TTree>(*f,"expErfThre");int id;double a,b;bind(t,"CellID",&id);bind(t,"threshold",&a);bind(t,"sigma",&b);for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(extra.count(id)){extra.at(id).threshold=a;extra.at(id).thresholdSigma=b;}}}
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 std::map<int,std::array<Channel,2>> channels;for(auto&[id,p]:pars)channels[id]={};
 std::array<std::array<EventSum,2>,2> sums; // all, shower + CoG20, positive good hits
 std::array<std::array<std::unique_ptr<TH1D>,3>,2> hist;
 for(int s=0;s<2;++s){auto n=std::string(s?"mc":"data");hist[s][0]=std::make_unique<TH1D>((n+"_lg_adc").c_str(),";Raw LG ADC;Hits (unit area)",205,0,4100);hist[s][1]=std::make_unique<TH1D>((n+"_hg_adc").c_str(),";Raw HG ADC;Hits (unit area)",205,0,4100);hist[s][2]=std::make_unique<TH1D>((n+"_lg_energy").c_str(),";Selected LG hit energy [MeV];Hits (unit area)",200,0,400);}
 std::ofstream files(out/"inputs/sampled_files.tsv");files<<"sample\tfile_index\tentries\tsampled\traw\tcalib\n";
 std::ofstream runTable(out/"runs.tsv");runTable<<"sample\tfile_index\trun\thits\tmax_hg_error\tmax_lg_error\tmean_temperature\n";
 long long mismatch[2]={},ncheck[2]={},gainTag0[2]={};double maxH[2]={},maxL[2]={};long long nOff[2]={};
 double ref=0;{auto f=openFile((fs::path(o["input"])/"study.root").string());auto*t=object<TTree>(*f,"paired_events");double e;bind(t,"nominal_GeV",&e);std::vector<double> values;for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);values.push_back(e);}std::sort(values.begin(),values.end());require(!values.empty(),"Empty reference MC");ref=values[values.size()/2];}
 for(int s=0;s<2;++s){std::ifstream manifest{fs::path(o.at("input"))/"inputs"/(s?"mc.tsv":"data.tsv")};require(bool(manifest),"Missing manifest");std::string line;int fi=0;
  while(std::getline(manifest,line)){if(s&&fi>=std::stoi(o["mc-files"]))break;auto tab=line.find('\t');require(tab!=std::string::npos,"Bad manifest");auto raw=openFile(line.substr(0,tab)),cal=openFile(line.substr(tab+1));auto*r=object<TTree>(*raw,"Raw_Hit"),*c=object<TTree>(*cal,"Calib_Hit");require(r->GetEntries()==c->GetEntries(),"Entry mismatch");r->SetBranchStatus("*",0);c->SetBranchStatus("*",0);
   int rr,rt,re,cr,ct,ce;std::vector<int>*ri=nullptr,*tag=nullptr,*gain=nullptr,*ci=nullptr;std::vector<double>*rh=nullptr,*rl=nullptr,*eh=nullptr,*el=nullptr,*en=nullptr,*temp=nullptr,*x=nullptr,*y=nullptr;
   bind(r,"Run_Num",&rr);bind(r,"Event_Time",&rt);bind(r,"TriggerID",&re);bind(r,"CellID",&ri);bind(r,"HitTag",&tag);bind(r,"GainTag",&gain);bind(r,"HG_Charge",&rh);bind(r,"LG_Charge",&rl);
   bind(c,"Run_Num",&cr);bind(c,"Event_Time",&ct);bind(c,"Event_Num",&ce);bind(c,"CellID",&ci);bind(c,"Hit_HG_Energy",&eh);bind(c,"Hit_LG_Energy",&el);bind(c,"Hit_Energy",&en);bind(c,"NewTemperature",&temp);bind(c,"Hit_X",&x);bind(c,"Hit_Y",&y);
   for(auto*t:{r,c}){t->SetCacheSize(16*1024*1024);t->AddBranchToCache("*",true);t->StopCacheLearningPhase();}
   Long64_t count=c->GetEntries(),take=std::min<Long64_t>(count,std::stoll(o[s?"mc-per-file":"data-per-file"]));int blocks=std::min<Long64_t>(10,take);Stats runT;long long runHits=0;double runH=0,runL=0;
   for(int b=0;b<blocks;++b){Long64_t start=count*b/blocks,end=count*(b+1)/blocks,n=take/blocks+(b<take%blocks);start+=(end-start-n)/2;
    for(Long64_t entry=start;entry<start+n;++entry){require(r->GetEntry(entry)>0&&c->GetEntry(entry)>0,"Read failed");require(rr==cr&&rt==ct&&re==ce,"Event key mismatch");require(ri->size()==rh->size()&&ri->size()==rl->size()&&ri->size()==tag->size()&&ri->size()==gain->size(),"Raw sizes mismatch");require(ci->size()==en->size()&&ci->size()==eh->size()&&ci->size()==el->size()&&ci->size()==temp->size()&&ci->size()==x->size()&&ci->size()==y->size(),"Cal sizes mismatch");
     std::unordered_map<int,size_t> index;index.reserve(ri->size());
     for(size_t k=0;k<ri->size();++k){int id=physical(ri->at(k));if(!connected(id))continue;auto&p=pars.at(id);auto&ex=extra.at(id);auto&ch=channels.at(id)[s];++ch.raw;if(tag->at(k)!=1){++ch.off;++nOff[s];double qh=rh->at(k)-p.ph,ql=rl->at(k)-p.pl;ch.offH.add(qh);ch.offL.add(ql);if(ex.sh>0&&ex.sl>0&&std::abs(qh)<5*ex.sh&&std::abs(ql)<5*ex.sl){ch.offCoreH.add(qh);ch.offCoreL.add(ql);ch.offLine.add(qh,ql);}}else require(index.emplace(ri->at(k),k).second,"Duplicate tagged full CellID");}
     double positive=0,good=0,negative=0,lgE=0,clipE=0,order=0;int hits=0,nlg=0;std::array<double,30> le{},lx{},ly{};
     for(size_t j=0;j<ci->size();++j){int id=physical(ci->at(j));if(!connected(id))continue;require(index.count(ci->at(j)),"Cal hit missing raw match");size_t k=index.at(ci->at(j));auto&p=pars.at(id);auto&ex=extra.at(id);auto&ch=channels.at(id)[s];double h=rh->at(k),l=rl->at(k),qh=h-p.ph,ql=l-p.pl,t=temp->at(j),C=1-(t-20)*(group(id)?1.6/135:3.5/230),f=.305/p.mip;
      double expectedH=f*C*qh,expectedL=f*(C*ql*p.gain+p.intercept),dh=eh->at(j)-expectedH,dl=el->at(j)-expectedL;
      ch.hgError.add(dh);ch.lgError.add(dl);maxH[s]=std::max(maxH[s],std::abs(dh));maxL[s]=std::max(maxL[s],std::abs(dl));runH=std::max(runH,std::abs(dh));runL=std::max(runL,std::abs(dl));++ncheck[s];mismatch[s]+=std::abs(dh)>1e-5+1e-6*std::abs(expectedH)||std::abs(dl)>1e-5+1e-6*std::abs(expectedL);
      bool lg=h>=2600;require(std::abs(en->at(j)-(lg?el->at(j):eh->at(j)))<1e-9,"Saved gain selection differs");++ch.hits;ch.lg+=lg;ch.negative+=en->at(j)<0;ch.clipH+=h==4000;ch.clipL+=l==3000;ch.overL+=l>3000;ch.belowSwitch+=h==4000&&ql*ex.digiA+p.ph<2600;ch.lowHG+=qh<100;ch.gain0+=gain->at(k)==0;gainTag0[s]+=gain->at(k)==0;ch.e+=en->at(j);if(lg)ch.lgE+=en->at(j);if(l==3000)ch.clipLE+=en->at(j);ch.temperature.add(t);runT.add(t);++runHits;if(std::abs(eh->at(j))>1e-6)ch.mipInferred.add(.305*C*qh/eh->at(j));
      if(h>=800&&h<2200)ch.overlap.add(qh,ql);
      if(h>=800&&h<1400)ch.low.add(qh,ql);
      if(h>=1400&&h<2200)ch.high.add(qh,ql);
      if(h>=300&&h<3000){ch.residual[band(h)].add(ql-qh/p.gain);ch.energy[band(h)].add(el->at(j)-eh->at(j));}
      double shift=lg?(C-1)*p.offset():0;ch.tempAfterDelta+=shift;
      if(en->at(j)<0){negative+=en->at(j)/1000;continue;}positive+=en->at(j)/1000;if(bad.count(id))continue;good+=en->at(j)/1000;++hits;if(lg){lgE+=en->at(j)/1000;++nlg;}if(l==3000)clipE+=en->at(j)/1000;order+=shift/1000;int layer=id/100000;le[layer]+=en->at(j);lx[layer]+=x->at(j)*en->at(j);ly[layer]+=y->at(j)*en->at(j);
     }
     bool selected=positive>.5*ref&&positive<1.5*ref;for(int layer:{9,10})selected=selected&&le[layer]>0&&std::abs(lx[layer]/le[layer])<20&&std::abs(ly[layer]/le[layer])<20;
     for(int selection=0;selection<2;++selection)if(selection==0||selected){auto&a=sums[s][selection];a.energy.add(good);a.hits.add(hits);a.lgEnergy.add(lgE);a.lgHits.add(nlg);a.thermalOrder.add(order);a.clipEnergy.add(clipE);a.negative.add(negative);}
     if(selected)for(size_t j=0;j<ci->size();++j){int id=physical(ci->at(j));if(!connected(id)||bad.count(id)||en->at(j)<0)continue;size_t k=index.at(ci->at(j));hist[s][0]->Fill(rl->at(k));hist[s][1]->Fill(rh->at(k));if(rh->at(k)>=2600)hist[s][2]->Fill(en->at(j));}
    }
   }
   files<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<count<<'\t'<<take<<'\t'<<line<<'\n';runTable<<(s?"mc":"data")<<'\t'<<fi<<'\t'<<rr<<'\t'<<runHits<<'\t'<<runH<<'\t'<<runL<<'\t'<<runT.mean<<'\n';std::cout<<(s?"mc":"data")<<" file "<<fi<<" sampled "<<take<<" max calibration errors "<<runH<<","<<runL<<std::endl;
   r->ResetBranchAddresses();c->ResetBranchAddresses();delete ri;delete tag;delete gain;delete ci;delete rh;delete rl;delete eh;delete el;delete en;delete temp;delete x;delete y;++fi;
  }
 }
 std::ofstream table(out/"channels.tsv");table<<std::setprecision(12)<<"sample\tcell\tinner10um\tbad\tmip\tmip_fallback\tcalib_HL_A\tcalib_HL_B\tdigi_HL_A\traw_HL_b\tHL_present\tMC_HG_switch\tped_HG\tped_LG\tped_sigma_HG\tped_sigma_LG\tthreshold\tthreshold_sigma\traw_entries\thittag0\tcal_hits\tLG_hits\tnegative_hits\tHG_eq4000\tLG_eq3000\tLG_gt3000\tHG4000_impliedHG_lt2600\tHG_sub_lt100\tgain_tag0\tenergy_MeV\tLG_energy_MeV\tLG3000_energy_MeV\ttemperature\tinferred_MIP\tHG_error_mean\tHG_error_RMS\tLG_error_mean\tLG_error_RMS\tcore_ped_n\tcore_ped_HG_mean\tcore_ped_HG_sigma\tcore_ped_LG_mean\tcore_ped_LG_sigma\toverlap_n\tfit_slope\tfit_intercept\tfit_sigma_LG\tfit_slope_error\tlow_n\tlow_slope\thigh_n\thigh_slope\tthermal_order_delta_MeV";
 for(int b=0;b<5;++b)table<<"\tband"<<b<<"_n\tband"<<b<<"_raw_residual_mean\tband"<<b<<"_raw_residual_sigma\tband"<<b<<"_ELminusEH_mean";
 table<<'\n';
 for(auto&[id,ch]:channels)for(int s=0;s<2;++s){auto&a=ch[s];auto&p=pars.at(id);auto&ex=extra.at(id);table<<(s?"mc":"data")<<'\t'<<id<<'\t'<<group(id)<<'\t'<<bad.count(id)<<'\t'<<p.mip<<'\t'<<p.mipFallback<<'\t'<<p.gain<<'\t'<<p.intercept<<'\t'<<ex.digiA<<'\t'<<ex.b<<'\t'<<ex.hl<<'\t'<<ex.sw<<'\t'<<p.ph<<'\t'<<p.pl<<'\t'<<ex.sh<<'\t'<<ex.sl<<'\t'<<ex.threshold<<'\t'<<ex.thresholdSigma<<'\t'<<a.raw<<'\t'<<a.off<<'\t'<<a.hits<<'\t'<<a.lg<<'\t'<<a.negative<<'\t'<<a.clipH<<'\t'<<a.clipL<<'\t'<<a.overL<<'\t'<<a.belowSwitch<<'\t'<<a.lowHG<<'\t'<<a.gain0<<'\t'<<a.e<<'\t'<<a.lgE<<'\t'<<a.clipLE<<'\t'<<a.temperature.mean<<'\t'<<a.mipInferred.mean<<'\t'<<a.hgError.mean<<'\t'<<a.hgError.rms()<<'\t'<<a.lgError.mean<<'\t'<<a.lgError.rms()<<'\t'<<a.offCoreH.n<<'\t'<<a.offCoreH.mean<<'\t'<<a.offCoreH.rms()<<'\t'<<a.offCoreL.mean<<'\t'<<a.offCoreL.rms()<<'\t'<<a.overlap.n<<'\t'<<a.overlap.slope()<<'\t'<<a.overlap.intercept()<<'\t'<<a.overlap.width()<<'\t'<<a.overlap.se()<<'\t'<<a.low.n<<'\t'<<a.low.slope()<<'\t'<<a.high.n<<'\t'<<a.high.slope()<<'\t'<<a.tempAfterDelta;for(int b=0;b<5;++b)table<<'\t'<<a.residual[b].n<<'\t'<<a.residual[b].mean<<'\t'<<a.residual[b].rms()<<'\t'<<a.energy[b].mean;table<<'\n';}
 std::ofstream val(out/"validation.tsv");val<<std::setprecision(12)<<"sample\tvalidated_hits\tformula_mismatch_hits\tmax_HG_error_MeV\tmax_LG_error_MeV\traw_hittag0\tgain_tag0\n";for(int s=0;s<2;++s)val<<(s?"mc":"data")<<'\t'<<ncheck[s]<<'\t'<<mismatch[s]<<'\t'<<maxH[s]<<'\t'<<maxL[s]<<'\t'<<nOff[s]<<'\t'<<gainTag0[s]<<'\n';
 std::ofstream ev(out/"event_summary.tsv");ev<<std::setprecision(12)<<"sample\tselection\tevents\tgood_energy_GeV\tgood_hits\tLG_energy_GeV\tLG_hits\tLG3000_energy_GeV\ttemperature_after_HL_delta_GeV\tnegative_energy_GeV\n";for(int s=0;s<2;++s)for(int k=0;k<2;++k){auto&a=sums[s][k];ev<<(s?"mc":"data")<<'\t'<<(k?"shower_cog20":"all")<<'\t'<<a.energy.n<<'\t'<<a.energy.mean<<'\t'<<a.hits.mean<<'\t'<<a.lgEnergy.mean<<'\t'<<a.lgHits.mean<<'\t'<<a.clipEnergy.mean<<'\t'<<a.thermalOrder.mean<<'\t'<<a.negative.mean<<'\n';}
 TFile rf((out/"raw_calibration.root").c_str(),"RECREATE");TCanvas canvas("adc","",1500,500);canvas.Divide(3,1);for(int j=0;j<3;++j){canvas.cd(j+1);gPad->SetLogy();double ymax=0;for(int s=0;s<2;++s){hist[s][j]->Write();if(hist[s][j]->Integral()>0)hist[s][j]->Scale(1./hist[s][j]->Integral());ymax=std::max(ymax,hist[s][j]->GetMaximum());}TLegend leg(.55,.75,.89,.89);for(int s=0;s<2;++s){hist[s][j]->SetLineColor(s?kBlue+1:kBlack);hist[s][j]->SetMinimum(1e-7);hist[s][j]->SetMaximum(ymax*2);hist[s][j]->Draw(s?"HIST SAME":"HIST");leg.AddEntry(hist[s][j].get(),s?"MC":"data","l");}leg.DrawClone();}canvas.SaveAs((out/"figures/raw_adc_and_lg_energy.png").c_str());
 std::cout<<"Finished: "<<out<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
