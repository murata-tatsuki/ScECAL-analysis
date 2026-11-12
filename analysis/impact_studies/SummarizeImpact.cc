#include "StudyCommon.hh"
#include <TParameter.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <array>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
struct Row {double beam;std::array<double,3> mean,sem,rms,resolution;double delta,deltaSem;};
int main(int argc,char**argv)try{
 require(argc>=3,"Usage: SummarizeImpact OUTPUT_DIRECTORY COMPLETED_STUDY_ROOT...");
 std::vector<Row> rows;std::set<double> seen;
 for(int i=2;i<argc;++i){auto f=openFile(argv[i]);Row r{};r.beam=object<TParameter<double>>(*f,"beam_energy_GeV")->GetVal();require(seen.insert(r.beam).second,"Duplicate beam energy");
  std::ifstream s(fs::path(argv[i]).parent_path()/"summary.tsv");std::string line;std::getline(s,line);
  for(int j=0;j<3;++j){require(bool(std::getline(s,line)),"Incomplete summary");std::istringstream in(line);std::string name;long long n;in>>name>>n>>r.mean[j]>>r.sem[j]>>r.rms[j]>>r.resolution[j];require(!in.fail()&&n>0&&r.mean[j]>0,"Invalid summary row");require(name==std::array<std::string,3>{"data","mc_nominal","mc_intercept_zero"}[j],"Unexpected sample order");}
  std::ifstream p(fs::path(argv[i]).parent_path()/"paired_effect.tsv");std::getline(p,line);require(bool(std::getline(p,line)),"Missing paired effect");std::istringstream in(line);long long n;in>>n>>r.delta>>r.deltaSem;require(!in.fail()&&n>0,"Invalid paired effect");rows.push_back(r);
 }
 std::sort(rows.begin(),rows.end(),[](const Row&a,const Row&b){return a.beam<b.beam;});fs::path out=argv[1];fs::create_directories(out);
 gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetTitleFont(42,"XYZ");gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleSize(.04,"XYZ");gStyle->SetLabelSize(.035,"XYZ");
 std::ofstream table(out/"energy_scan.tsv");table<<std::setprecision(12)<<"beam_GeV\tdata_mean_GeV\tmc_nominal_mean_GeV\tmc_zero_mean_GeV\tnominal_over_data\tzero_over_data\tpaired_delta_GeV\tpaired_delta_sem_GeV\n";
 TGraphErrors means[3],ratios[2],res[3],effect;
 for(size_t i=0;i<rows.size();++i){const auto&r=rows[i];for(int s=0;s<3;++s){means[s].SetPoint(i,r.beam,r.mean[s]);means[s].SetPointError(i,0,r.sem[s]);res[s].SetPoint(i,r.beam,r.resolution[s]);}
  for(int s=0;s<2;++s){double ratio=r.mean[s+1]/r.mean[0];ratios[s].SetPoint(i,r.beam,ratio);ratios[s].SetPointError(i,0,std::hypot(r.sem[s+1]/r.mean[0],ratio*r.sem[0]/r.mean[0]));}
  effect.SetPoint(i,r.beam,r.delta);effect.SetPointError(i,0,r.deltaSem);table<<r.beam<<'\t'<<r.mean[0]<<'\t'<<r.mean[1]<<'\t'<<r.mean[2]<<'\t'<<r.mean[1]/r.mean[0]<<'\t'<<r.mean[2]/r.mean[0]<<'\t'<<r.delta<<'\t'<<r.deltaSem<<'\n';
 }
 std::array<int,3> colors={kBlack,kBlue+1,kRed+1};std::array<std::string,3> names={"data","mc_nominal","mc_intercept_zero"};
 auto plot=[&](TGraphErrors*g,int count,int offset,const char*title,const char*file){TCanvas c("scan","",1000,650);c.SetLogx();TLegend leg(.14,.72,.54,.88);leg.SetTextFont(42);leg.SetTextSize(.03);double lo=0,hi=0;for(int s=0;s<count;++s)for(int j=0;j<g[s].GetN();++j){lo=std::min(lo,g[s].GetPointY(j)-g[s].GetErrorY(j));hi=std::max(hi,g[s].GetPointY(j)+g[s].GetErrorY(j));}double range=std::max(.001,hi-lo);
  for(int s=0;s<count;++s){g[s].SetTitle(title);g[s].SetLineColor(colors[s+offset]);g[s].SetMarkerColor(colors[s+offset]);g[s].SetMarkerStyle(20+s);g[s].SetMinimum(lo-.1*range);g[s].SetMaximum(hi+.25*range);g[s].Draw(s?"LP SAME":"ALP");leg.AddEntry(&g[s],names[s+offset].c_str(),"lp");}leg.Draw();c.SaveAs((out/file).c_str());};
 plot(means,3,0,";Beam energy [GeV];Mean event energy [GeV]","mean_energy_vs_beam.png");plot(ratios,2,1,";Beam energy [GeV];MC mean / data mean","mc_over_data_vs_beam.png");plot(res,3,0,";Beam energy [GeV];RMS / mean (no fit)","rms_over_mean_vs_beam.png");
 TCanvas c("effect","",1000,650);c.SetLogx();effect.SetTitle(";Beam energy [GeV];Paired mean E_{zero}-E_{nominal} [GeV]");effect.SetMarkerStyle(20);effect.Draw("ALP");c.SaveAs((out/"paired_energy_change_vs_beam.png").c_str());
 std::cout<<"Summarized "<<rows.size()<<" completed energies in "<<out<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
