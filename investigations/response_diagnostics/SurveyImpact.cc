#include "EventHistograms.hh"
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TParameter.h>
#include <iostream>
struct Summary {Stats all,high;double low=0,negative=0,median=0,peak=0,error=0,sigma=0,chi=0;int status=-1;};
Summary summarize(const std::vector<double>&v,double ref){
 Summary s;std::vector<double> high;TH1D h("peak_counts","",640,0,1.6*ref);h.SetDirectory(nullptr);
 for(double x:v){s.all.add(x);s.low+=x<.1*ref;s.negative+=x<0;if(x>.5*ref&&x<1.5*ref){s.high.add(x);high.push_back(x);}h.Fill(x);}
 std::sort(high.begin(),high.end());s.median=eventplots::quantile(high,.5);s.low/=v.size();s.negative/=v.size();
 int bmax=0;double best=0;for(int b=h.FindBin(.5*ref);b<=h.FindBin(1.4*ref);++b){double n=h.GetBinContent(b-1)+h.GetBinContent(b)+h.GetBinContent(b+1);if(n>best){best=n;bmax=b;}}
 if(s.high.n<100||!bmax)return s;
 double mode=h.GetBinCenter(bmax);int l=bmax,r=bmax;while(l>1&&h.GetBinContent(l)>best/6)--l;while(r<h.GetNbinsX()&&h.GetBinContent(r)>best/6)++r;
 double width=std::max(.01*ref,(r-l)*h.GetBinWidth(1)/2.355);TF1 fit("core_gaussian","gaus",mode-1.5*width,mode+1.5*width);fit.SetParameters(h.GetBinContent(bmax),mode,width);fit.SetParLimits(0,0,v.size());fit.SetParLimits(1,.4*ref,1.5*ref);fit.SetParLimits(2,.002*ref,.4*ref);
 for(int k=0;k<2;++k){fit.SetRange(mode-1.5*width,mode+1.5*width);s.status=int(h.Fit(&fit,"Q0RSN"));mode=fit.GetParameter(1);width=fit.GetParameter(2);}
 s.peak=mode;s.error=fit.GetParError(1);s.sigma=width;s.chi=fit.GetNDF()>0?fit.GetChisquare()/fit.GetNDF():0;return s;
}
void row(std::ostream&o,double beam,const std::string&sample,const Summary&s){o<<beam<<'\t'<<sample<<'\t'<<s.all.n<<'\t'<<s.all.mean<<'\t'<<s.low<<'\t'<<s.negative<<'\t'<<s.high.n<<'\t'<<s.high.mean<<'\t'<<s.median<<'\t'<<s.peak<<'\t'<<s.error<<'\t'<<s.sigma<<'\t'<<s.status<<'\t'<<s.chi;}
int main(int argc,char**argv)try{
 auto o=options(argc,argv,{{"base",""},{"out",""}});fs::path base=o["base"],out=o["out"];fs::create_directories(out/"figures");gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 std::ofstream table(out/"event_populations.tsv"),runs(out/"data_runs.tsv"),ratios(out/"peak_ratios.tsv");table<<std::setprecision(12);runs<<std::setprecision(12);ratios<<std::setprecision(12);
 std::string header="beam_GeV\tsample\tevents\tmean_GeV\tfraction_below_0.1_mc_median\tfraction_negative\tcore_events\tcore_mean_GeV\tcore_median_GeV\tgaussian_mean_GeV\tgaussian_mean_error\tgaussian_sigma_GeV\tfit_status\tchi2_ndf";
 table<<header<<'\n';runs<<header<<"\tfile\n";ratios<<"beam_GeV\tdata_peak\tmc_peak\tzero_peak\tdata_over_mc_peak\tdata_over_zero_peak\tpaired_mean_relative_shift\tdata_low_fraction\tmc_low_fraction\tdata_over_mc_core_median\n";
 std::vector<std::pair<double,fs::path>> dirs;for(auto&entry:fs::directory_iterator(base))if(entry.is_directory()&&fs::exists(entry.path()/"event_histograms.root")){auto name=entry.path().filename().string();dirs.emplace_back(std::stod(name),entry.path());}std::sort(dirs.begin(),dirs.end());
 TGraphErrors nominalRatio,zeroRatio;TGraph lowData,lowMC;
 for(const auto&[beam,dir]:dirs){auto cache=openFile((dir/"event_histograms.root").string()),study=openFile((dir/"study.root").string());eventplots::Values v;auto*d=object<TTree>(*cache,"data_events");double energy;bind(d,"energy_GeV",&energy);for(Long64_t i=0;i<d->GetEntries();++i){d->GetEntry(i);v[0].push_back(energy);}d->ResetBranchAddresses();auto*m=object<TTree>(*study,"paired_events");double a,b;m->SetBranchStatus("*",0);bind(m,"nominal_GeV",&a);bind(m,"intercept_zero_GeV",&b);Stats shift;for(Long64_t i=0;i<m->GetEntries();++i){m->GetEntry(i);v[1].push_back(a);v[2].push_back(b);if(a!=0)shift.add((b-a)/a);}m->ResetBranchAddresses();auto sorted=v[1];std::sort(sorted.begin(),sorted.end());double ref=eventplots::quantile(sorted,.5);require(ref>0,"Nonpositive MC median");std::array<Summary,3>s;
  for(int k=0;k<3;++k){s[k]=summarize(v[k],ref);row(table,beam,eventplots::samples[k],s[k]);table<<'\n';}
  std::ifstream manifest(dir/"inputs/data.tsv");std::string line;size_t cursor=0;
  while(std::getline(manifest,line)&&cursor<v[0].size()){auto path=line.substr(line.find('\t')+1);auto f=openFile(path);size_t count=std::min<size_t>(object<TTree>(*f,"Calib_Hit")->GetEntries(),v[0].size()-cursor);std::vector<double>part(v[0].begin()+cursor,v[0].begin()+cursor+count);row(runs,beam,"data",summarize(part,ref));runs<<'\t'<<fs::path(path).filename().string()<<'\n';cursor+=count;}require(cursor==v[0].size(),"Data manifest count mismatch");
  double nr=s[0].peak/s[1].peak,zr=s[0].peak/s[2].peak;ratios<<beam<<'\t'<<s[0].peak<<'\t'<<s[1].peak<<'\t'<<s[2].peak<<'\t'<<nr<<'\t'<<zr<<'\t'<<shift.mean<<'\t'<<s[0].low<<'\t'<<s[1].low<<'\t'<<s[0].median/s[1].median<<'\n';
  int n=nominalRatio.GetN();nominalRatio.SetPoint(n,beam,nr);nominalRatio.SetPointError(n,0,nr*std::hypot(s[0].error/s[0].peak,s[1].error/s[1].peak));zeroRatio.SetPoint(n,beam,zr);zeroRatio.SetPointError(n,0,zr*std::hypot(s[0].error/s[0].peak,s[2].error/s[2].peak));lowData.SetPoint(n,beam,s[0].low);lowMC.SetPoint(n,beam,s[1].low);
  std::cout<<beam<<" GeV: peak data/MC="<<nr<<", data/zero="<<zr<<", low fraction="<<s[0].low<<std::endl;
 }
 TCanvas c("scan","",1100,650);c.SetLogx();c.SetGrid();nominalRatio.SetTitle("Gaussian core: diagnostic fit, no event identification;Beam energy [GeV];Data / MC peak");nominalRatio.SetMarkerStyle(20);nominalRatio.SetLineColor(kBlue+1);nominalRatio.SetMarkerColor(kBlue+1);nominalRatio.Draw("APL");nominalRatio.GetYaxis()->SetRangeUser(.4,1.4);zeroRatio.SetMarkerStyle(21);zeroRatio.SetMarkerColor(kRed+1);zeroRatio.SetLineColor(kRed+1);zeroRatio.Draw("PL SAME");TLegend leg(.15,.72,.44,.86);leg.AddEntry(&nominalRatio,"MC nominal","lp");leg.AddEntry(&zeroRatio,"MC intercept zero","lp");leg.Draw();c.SaveAs((out/"figures/peak_ratio_vs_beam.png").c_str());c.Clear();c.SetLogx();c.SetGrid();lowData.SetTitle("E_{sum} < 0.1 #times median(E_{MC});Beam energy [GeV];Fraction of all events");lowData.SetMarkerStyle(20);lowData.Draw("APL");lowData.GetYaxis()->SetRangeUser(0,.8);lowMC.SetMarkerStyle(21);lowMC.SetLineColor(kBlue+1);lowMC.SetMarkerColor(kBlue+1);lowMC.Draw("PL SAME");TLegend lowLegend(.15,.72,.4,.86);lowLegend.AddEntry(&lowData,"Data","lp");lowLegend.AddEntry(&lowMC,"MC nominal","lp");lowLegend.Draw();c.SaveAs((out/"figures/low_energy_fraction_vs_beam.png").c_str());return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
