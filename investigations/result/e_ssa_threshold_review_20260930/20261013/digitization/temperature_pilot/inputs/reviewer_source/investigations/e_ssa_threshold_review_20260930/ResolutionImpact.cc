// Resolution diagnostics from validated, paired raw/calibrated event tables.
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>
namespace fs=std::filesystem;
void need(bool ok,const std::string&m){if(!ok)throw std::runtime_error(m);}
struct Event{int file,run;long long entry;double x,y,e,central,hits,chits,temp,digi,valid,clean,lg;};
using Events=std::vector<Event>;
using Values=std::vector<double>;
struct Width{double mean,rms,core,gaus;int fit;};
double quantile(Values a,double p){std::sort(a.begin(),a.end());double at=p*(a.size()-1);size_t i=at;return a[i]+(at-i)*(a[std::min(i+1,a.size()-1)]-a[i]);}
double mean(const Values&a){return std::accumulate(a.begin(),a.end(),0.)/a.size();}
double rms(const Values&a){double m=mean(a),v=0;for(double x:a)v+=(x-m)*(x-m);return std::sqrt(v/a.size())/m;}
Width width(const Values&a,TH1D*hist=nullptr){need(a.size()>20,"Too few selected events");Width w{mean(a),rms(a),(quantile(a,.84)-quantile(a,.16))/(2*quantile(a,.5)),0,-1};if(hist){TF1 f("diagnostic_gaussian","gaus",.015,10);f.SetParameters(hist->GetMaximum(),w.mean,w.mean*w.rms);auto r=hist->Fit(&f,"QSN","",std::max(.015,w.mean-3*w.mean*w.rms),w.mean+3*w.mean*w.rms);double m=f.GetParameter(1),s=std::abs(f.GetParameter(2));if(int(r)==0&&m>0&&s>0){r=hist->Fit(&f,"QSN","",std::max(.015,m-1.5*s),m+1.5*s);w.fit=int(r);w.gaus=std::abs(f.GetParameter(2))/f.GetParameter(1);}}return w;}
unsigned long long block(const Event&e){return (static_cast<unsigned long long>(e.file)<<32)+e.entry/100;}
int fold(const Event&e){auto z=block(e)+0x9e3779b97f4a7c15ULL;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return (z^(z>>31))%5;}
Values features(const Event&e,int mode,const std::vector<int>&runs){Values x{1};if(mode!=0){double a=e.x/20,b=e.y/20;x.insert(x.end(),{a,b,a*a,a*b,b*b});}if(mode!=1)for(size_t i=1;i<runs.size();++i)x.push_back(e.run==runs[i]);return x;}
Values corrected(const Events&ev,int mode,std::ostream&models,const std::string&sample,int energy){std::vector<int>runs;for(const auto&e:ev)runs.push_back(e.run);std::sort(runs.begin(),runs.end());runs.erase(std::unique(runs.begin(),runs.end()),runs.end());Values result(ev.size());int p=features(ev.front(),mode,runs).size();for(int held=0;held<5;++held){int n=0;double m=0;for(const auto&e:ev)if(fold(e)!=held){++n;m+=e.e;}m/=n;TMatrixD design(n,p);TVectorD response(n);int row=0;for(const auto&e:ev)if(fold(e)!=held){auto x=features(e,mode,runs);for(int j=0;j<p;++j)design(row,j)=x[j];response[row++]=e.e/m;}TDecompSVD solver(design);Bool_t ok;auto beta=solver.Solve(response,ok);need(ok,"Position/run regression singular");for(int j=0;j<p;++j)models<<energy<<'\t'<<sample<<'\t'<<mode<<'\t'<<held<<'\t'<<j<<'\t'<<beta[j]<<'\n';for(size_t i=0;i<ev.size();++i)if(fold(ev[i])==held){auto x=features(ev[i],mode,runs);double pred=0;for(int j=0;j<p;++j)pred+=x[j]*beta[j];need(pred>.4&&pred<1.6,"Invalid held-out response prediction");result[i]=ev[i].e/pred;}}return result;}
std::array<Events,2> read(const fs::path&path){std::ifstream f(path);need(bool(f),"Missing "+path.string());std::string line;std::getline(f,line);std::array<Events,2> all;while(std::getline(f,line)){std::istringstream r(line);std::string sample;Event e;int selected,tail;r>>sample>>e.file>>e.entry>>e.run>>selected>>tail>>e.x>>e.y>>e.e>>e.central>>e.hits>>e.chits>>e.temp>>e.digi>>e.valid>>e.clean>>e.lg;need(bool(r),"Malformed event row");if(selected&&!tail){need(e.e>0,"Invalid energy");all[sample=="mc"].push_back(e);}}return all;}
std::vector<std::vector<size_t>> blocks(const Events&ev){std::map<unsigned long long,std::vector<size_t>> m;for(size_t i=0;i<ev.size();++i)m[block(ev[i])].push_back(i);std::vector<std::vector<size_t>> out;for(auto&kv:m)out.push_back(kv.second);return out;}
void components(const Events&ev,int energy,int sample,std::ostream&out){double m=0;for(const auto&e:ev)m+=e.e;m/=ev.size();for(int mode=0;mode<2;++mode){double a=0,b=0,aa=0,bb=0,ab=0;for(const auto&e:ev){double x=mode?e.lg:e.central,y=e.e-x;a+=x;b+=y;aa+=x*x;bb+=y*y;ab+=x*y;}double n=ev.size();a/=n;b/=n;double va=(aa/n-a*a)/(m*m),vb=(bb/n-b*b)/(m*m),cov=2*(ab/n-a*b)/(m*m);out<<energy<<'\t'<<(sample?"mc":"data")<<'\t'<<(mode?"LG_HG":"central_peripheral")<<'\t'<<a/m<<'\t'<<va<<'\t'<<vb<<'\t'<<cov<<'\t'<<va+vb+cov<<'\n';}}
void fitScan(const fs::path&,const fs::path&);
int main(int argc,char**argv)try{need(argc==3||argc==4,"Usage: ResolutionImpact INPUT_BASE OUTPUT_BASE");fs::path input=argv[1],out=argv[2];if(argc==4){need(std::string(argv[3])=="--fit-scan","Invalid option");fitScan(input,out);return 0;}fs::create_directories(out);gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);std::ofstream summary(out/"summary.tsv"),impact(out/"impact.tsv"),models(out/"models.tsv"),parts(out/"components.tsv"),validation(out/"validation.tsv");for(auto*f:{&summary,&impact,&models,&parts,&validation})*f<<std::setprecision(12);
summary<<"energy_GeV\tsample\tcondition\tevents\tmean_GeV\tRMS_over_mean\tsigma68_over_median\tgaussian_sigma_over_mean\tfit_status\n";
impact<<"energy_GeV\tcondition\tdata_RMS\tmc_RMS\textra_RMS\tvariance_gap_reduction\treduction_q16\treduction_q84\tdata_RMS_delta_pp\tmc_RMS_delta_pp\n";
models<<"energy_GeV\tsample\tmode_0_run_1_position_2_both\theld_out_fold\tcoefficient_index\tcoefficient\n";
parts<<"energy_GeV\tsample\tpartition\tfirst_energy_fraction\tfirst_normalized_variance\tsecond_normalized_variance\ttwice_normalized_covariance\tsummed_relative_variance\n";
validation<<"energy_GeV\tsample\tevents\tmax_mean_difference_GeV\tmax_RMS_difference\n";
for(int energy:{5,20,100}){auto ev=read(input/"sensitivity/comparison/20mm"/(std::to_string(energy)+"GeV")/"events.tsv");std::array<std::map<std::string,Values>,2> vals;std::array<std::map<std::string,Width>,2> widths;
for(int s=0;s<2;++s){for(const auto&e:ev[s]){auto&v=vals[s];v["baseline"].push_back(e.e);v["no_temperature"].push_back(e.temp);v["digi_mip"].push_back(e.digi);v["clean_mip"].push_back(e.clean);if(std::abs(e.x)<5&&std::abs(e.y)<5)v["center5"].push_back(e.e);}for(int mode=0;mode<3;++mode)vals[s][std::array<std::string,3>{"run_cv","position_cv","position_run_cv"}[mode]]=corrected(ev[s],mode,models,s?"mc":"data",energy);components(ev[s],energy,s,parts);
// Independently reproduce cached summary before accepting any impact result.
std::ifstream old(input/"sensitivity"/(s?"mc":"data")/"exclude/20mm"/(std::to_string(energy)+"GeV")/"summary.tsv");need(bool(old),"Missing reference summary");std::string line;std::getline(old,line);std::getline(old,line);std::istringstream r(line);std::string name,selection;long long n;double sw,ne,m,w;r>>name>>selection>>n>>sw>>ne>>m>>w;double dm=std::abs(mean(vals[s]["baseline"])-m),dw=std::abs(rms(vals[s]["baseline"])-w);need(n==static_cast<long long>(ev[s].size())&&dm<1e-9&&dw<1e-9,"Cached event closure failed");validation<<energy<<'\t'<<name<<'\t'<<n<<'\t'<<dm<<'\t'<<dw<<'\n';
for(auto&[condition,v]:vals[s]){auto dest=out/(s?"mc":"data")/condition/(condition=="center5"?"5mm":"20mm")/(std::to_string(energy)+"GeV");fs::create_directories(dest);double mu=mean(v);TH1D h("energy",";Energy [GeV];Events",240,.3*mu,1.7*mu);for(double x:v)h.Fill(x);auto w=width(v,&h);widths[s][condition]=w;summary<<energy<<'\t'<<(s?"mc":"data")<<'\t'<<condition<<'\t'<<v.size()<<'\t'<<w.mean<<'\t'<<w.rms<<'\t'<<w.core<<'\t'<<w.gaus<<'\t'<<w.fit<<'\n';std::ofstream f(dest/"width.tsv");f<<std::setprecision(12)<<"events\tmean_GeV\tRMS_over_mean\tsigma68_over_median\tgaussian_sigma_over_mean\tfit_status\n"<<v.size()<<'\t'<<w.mean<<'\t'<<w.rms<<'\t'<<w.core<<'\t'<<w.gaus<<'\t'<<w.fit<<'\n';}
fs::path dest=out/(s?"mc":"data")/"position_run_cv/20mm"/(std::to_string(energy)+"GeV");std::ofstream f(dest/"events.tsv");f<<std::setprecision(12)<<"file\tentry\trun\tfold\tx_mm\ty_mm\tnominal_GeV\tposition_cv_GeV\tposition_run_cv_GeV\n";for(size_t i=0;i<ev[s].size();++i){auto&e=ev[s][i];f<<e.file<<'\t'<<e.entry<<'\t'<<e.run<<'\t'<<fold(e)<<'\t'<<e.x<<'\t'<<e.y<<'\t'<<e.e<<'\t'<<vals[s]["position_cv"][i]<<'\t'<<vals[s]["position_run_cv"][i]<<'\n';}}
double basegap=std::pow(widths[0]["baseline"].rms,2)-std::pow(widths[1]["baseline"].rms,2);auto b0=blocks(ev[0]),b1=blocks(ev[1]);
for(auto&[condition,v]:vals[0]){double d=widths[0][condition].rms,m=widths[1][condition].rms,gap=d*d-m*m;Values reductions;if(condition!="center5"){std::mt19937_64 rng(energy+4242);for(int rep=0;rep<400;++rep){std::array<double,2> nominal{},changed{};for(int s=0;s<2;++s){auto&bb=s?b1:b0;Values a,c;for(size_t k=0;k<bb.size();++k){auto&group=bb[std::uniform_int_distribution<size_t>(0,bb.size()-1)(rng)];for(size_t j:group){a.push_back(vals[s]["baseline"][j]);c.push_back(vals[s][condition][j]);}}nominal[s]=std::pow(rms(a),2);changed[s]=std::pow(rms(c),2);}double denom=nominal[0]-nominal[1];if(denom>0)reductions.push_back(1-(changed[0]-changed[1])/denom);}}impact<<energy<<'\t'<<condition<<'\t'<<d<<'\t'<<m<<'\t'<<std::sqrt(std::max(0.,gap))<<'\t'<<1-gap/basegap<<'\t'<<(reductions.empty()?NAN:quantile(reductions,.16))<<'\t'<<(reductions.empty()?NAN:quantile(reductions,.84))<<'\t'<<100*(d-widths[0]["baseline"].rms)<<'\t'<<100*(m-widths[1]["baseline"].rms)<<'\n';}
fs::path plots=out/"comparison/20mm"/(std::to_string(energy)+"GeV")/"figures";fs::create_directories(plots);TCanvas c("resolution_impact","",1400,600);c.Divide(2,1);for(int p=0;p<2;++p){c.cd(p+1);const std::string cond=p?"position_run_cv":"baseline";std::array<TH1D,2> hs={TH1D("data","",100,.5,1.5),TH1D("mc","",100,.5,1.5)};double ymax=0;for(int s=0;s<2;++s){auto&v=vals[s][cond];double mu=mean(v);for(double x:v)hs[s].Fill(x/mu);hs[s].Scale(1./v.size());hs[s].SetLineColor(s?kRed+1:kBlack);hs[s].SetLineWidth(2);ymax=std::max(ymax,hs[s].GetMaximum());}hs[0].SetTitle((std::to_string(energy)+" GeV, "+cond+";Energy / sample mean;Probability / bin").c_str());hs[0].SetMaximum(1.15*ymax);auto*a=hs[0].DrawCopy("HIST");auto*b=hs[1].DrawCopy("HIST SAME");TLegend l(.62,.74,.88,.88);l.AddEntry(a,"data, LG tail excluded","l");l.AddEntry(b,"MC e_ssa threshold","l");l.DrawClone();}c.SaveAs((plots/"resolution_position_impact.png").c_str());
c.Clear();c.Divide(2,1);for(int axis=0;axis<2;++axis){c.cd(axis+1);std::array<TProfile,2> ps={TProfile("data_pos","",16,-20,20),TProfile("mc_pos","",16,-20,20)};for(int s=0;s<2;++s){double mu=mean(vals[s]["baseline"]);for(auto&e:ev[s])ps[s].Fill(axis?e.y:e.x,e.e/mu);ps[s].SetLineColor(s?kRed+1:kBlack);ps[s].SetMarkerColor(s?kRed+1:kBlack);ps[s].SetMarkerStyle(s?24:20);}ps[0].SetTitle((std::to_string(energy)+" GeV;"+(axis?"Layer10 CoG y [mm]":"Layer9 CoG x [mm]")+";Mean energy / sample mean").c_str());ps[0].SetMinimum(.8);ps[0].SetMaximum(1.2);auto*a=ps[0].DrawCopy();auto*b=ps[1].DrawCopy("SAME");TLegend l(.63,.77,.88,.89);l.AddEntry(a,"data","lp");l.AddEntry(b,"MC","lp");l.DrawClone();}c.SaveAs((plots/"energy_position_profile.png").c_str());std::cout<<energy<<" GeV complete\n";
}return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}

// Bin/range sensitivity of Gaussian core width; no re-training of position model.
void fitScan(const fs::path&input,const fs::path&out){
 gROOT->SetBatch(true);TH1::AddDirectory(false);std::ofstream table(out/"gaussian_scan.tsv");table<<std::setprecision(12)<<"energy_GeV\tsample\tcondition\tbins\tfit_halfwidth_sigma\tstatus\tmean_GeV\tsigma_over_mean\tfit_error\tlow_3sigma_fraction\thigh_3sigma_fraction\n";
 for(int energy:{5,20,100}){auto ev=read(input/"sensitivity/comparison/20mm"/(std::to_string(energy)+"GeV")/"events.tsv");for(int s=0;s<2;++s){std::map<std::string,Values> variants;for(auto&e:ev[s])variants["baseline"].push_back(e.e);
  std::ifstream saved(out/(s?"mc":"data")/"position_run_cv/20mm"/(std::to_string(energy)+"GeV")/"events.tsv");need(bool(saved),"Missing CV events");std::string line;std::getline(saved,line);size_t index=0;while(std::getline(saved,line)){std::istringstream r(line);int file,run,fold;long long entry;double x,y,e,p,b;r>>file>>entry>>run>>fold>>x>>y>>e>>p>>b;need(bool(r)&&index<ev[s].size()&&file==ev[s][index].file&&entry==ev[s][index].entry,"CV event identity mismatch");variants["position_run_cv"].push_back(b);++index;}need(index==ev[s].size(),"Incomplete CV cache");
  for(auto&[condition,v]:variants)for(int bins:{120,240,360})for(double window:{1.5,2.}){double mu=mean(v),sd=mu*rms(v);TH1D h("fit_scan","",bins,.3*mu,1.7*mu);for(double x:v)h.Fill(x);TF1 f("scan_gaussian","gaus",.015,10);f.SetParameters(h.GetMaximum(),mu,sd);auto result=h.Fit(&f,"QSN","",std::max(.015,mu-3*sd),mu+3*sd);int status=result;mu=f.GetParameter(1);sd=std::abs(f.GetParameter(2));if(status==0){result=h.Fit(&f,"QSN","",std::max(.015,mu-window*sd),mu+window*sd);status=result;}mu=f.GetParameter(1);sd=std::abs(f.GetParameter(2));double ratio=sd/mu,error=NAN;if(status==0){double variance=result->CovMatrix(2,2)/(mu*mu)+sd*sd*result->CovMatrix(1,1)/std::pow(mu,4)-2*sd*result->CovMatrix(1,2)/std::pow(mu,3);error=std::sqrt(std::max(0.,variance));}int low=0,high=0;for(double x:v){low+=x<mu-3*sd;high+=x>mu+3*sd;}table<<energy<<'\t'<<(s?"mc":"data")<<'\t'<<condition<<'\t'<<bins<<'\t'<<window<<'\t'<<status<<'\t'<<mu<<'\t'<<ratio<<'\t'<<error<<'\t'<<double(low)/v.size()<<'\t'<<double(high)/v.size()<<'\n';}
 }}
}
