#pragma once
#include "StudyCommon.hh"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <array>
#include <algorithm>
#include <iomanip>

namespace eventplots {
inline const std::array<std::string,3> samples={"data","mc_nominal","mc_intercept_zero"};
inline const std::array<int,3> colors={kBlack,kBlue+1,kRed+1};
using Values=std::array<std::vector<double>,3>;
struct Axis { int bins; double low,high; };
inline double quantile(const std::vector<double>& sorted,double p) {
 if(sorted.empty())return 0;
 double index=p*(sorted.size()-1);size_t i=static_cast<size_t>(index);
 return sorted[i]+(index-i)*(sorted[std::min(i+1,sorted.size()-1)]-sorted[i]);
}
inline double niceWidth(double x) {
 double power=std::pow(10.,std::floor(std::log10(x)));
 for(double m:{1.,2.,2.5,5.,10.})if(m*power>=x)return m*power;
 return 10*power;
}
// Keep 99.8% of EACH sample in the main view, including zero as a reference.
// Use ~10 bins per robust MC sigma where possible, with a practical bin cap.
// A separate full-range view preserves all tails, including negative energies.
inline Axis axis(const std::vector<std::vector<double>>& sorted,size_t reference,bool full) {
 double lo=0,hi=0;
 for(const auto&v:sorted)if(!v.empty()){
  lo=std::min(lo,full?v.front():quantile(v,.001));
  hi=std::max(hi,full?v.back():quantile(v,.999));
 }
 double span=hi-lo;
 if(!(span>0)){lo=-.5e-6;hi=.5e-6;span=hi-lo;}
 lo-=.04*span;hi+=.04*span;span=hi-lo;
 const auto&ref=sorted.at(reference);
 double sigma=(quantile(ref,.84)-quantile(ref,.16))/2;
 double target=full?span/400.:std::min(span/160.,sigma>0?sigma/10:span/160.);
 double width=niceWidth(std::max(target,span/(full?1600.:800.)));
 lo=std::floor(lo/width)*width;hi=(std::floor(hi/width)+1)*width;
 int bins=std::max(1,static_cast<int>(std::lround((hi-lo)/width)));
 return {bins,lo,lo+bins*width};
}
inline std::unique_ptr<TH1D> histogram(const std::string&name,const std::string&title,const Axis&a,const std::vector<double>&v) {
 auto h=std::make_unique<TH1D>(name.c_str(),title.c_str(),a.bins,a.low,a.high);h->SetDirectory(nullptr);
 for(double x:v){require(std::isfinite(x),"Nonfinite event energy");h->Fill(x);}return h;
}
struct Histograms {
 std::array<std::unique_ptr<TH1D>,3> totals,fullTotals;
 std::unique_ptr<TH1D> change,relative,fullChange,fullRelative;
 explicit Histograms(const Values&values) {
  require(!values[0].empty()&&!values[1].empty()&&values[1].size()==values[2].size(),"Invalid event sample sizes");
  std::vector<std::vector<double>> sorted(values.begin(),values.end());
  for(auto&v:sorted)std::sort(v.begin(),v.end());
  auto central=axis(sorted,1,false),full=axis(sorted,1,true);
  for(int s=0;s<3;++s){
   totals[s]=histogram(samples[s]+"_event_energy",";Sum of calibrated hit energies [GeV];Events",central,values[s]);
   fullTotals[s]=histogram(samples[s]+"_event_energy_full_range",";Sum of calibrated hit energies [GeV];Events",full,values[s]);
  }
  std::vector<double> d,r;
  for(size_t i=0;i<values[1].size();++i){double x=values[2][i]-values[1][i];d.push_back(x);if(values[1][i]!=0)r.push_back(x/values[1][i]);}
  std::sort(d.begin(),d.end());std::sort(r.begin(),r.end());
  change=histogram("mc_event_energy_change",";E_{zero}-E_{nominal} [GeV];Events",axis({d},0,false),d);
  fullChange=histogram("mc_event_energy_change_full_range",change->GetTitle(),axis({d},0,true),d);
  relative=histogram("mc_event_relative_change",";(E_{zero}-E_{nominal})/E_{nominal};Events",axis({r},0,false),r);
  fullRelative=histogram("mc_event_relative_change_full_range",relative->GetTitle(),axis({r},0,true),r);
 }
 std::vector<TH1D*> all()const {
  std::vector<TH1D*> h;for(const auto&a:{&totals,&fullTotals})for(const auto&p:*a)h.push_back(p.get());
  for(auto*p:{change.get(),relative.get(),fullChange.get(),fullRelative.get()})h.push_back(p);
  return h;
 }
 void write()const { for(auto*h:all())h->Write(); }
 void metadata(const fs::path&out)const {
  std::ofstream f(out/"event_histogram_binning.tsv");require(bool(f),"Cannot write binning table");
  f<<std::setprecision(12)<<"histogram\tbins\txmin\txmax\tbin_width\tentries\tunderflow\toverflow\n";
  for(auto*h:all())f<<h->GetName()<<'\t'<<h->GetNbinsX()<<'\t'<<h->GetXaxis()->GetXmin()<<'\t'<<h->GetXaxis()->GetXmax()<<'\t'<<h->GetBinWidth(1)<<'\t'<<h->GetEntries()<<'\t'<<h->GetBinContent(0)<<'\t'<<h->GetBinContent(h->GetNbinsX()+1)<<'\n';
 }
 static void flowLabel(const TH1D&h,double y,int color=kBlack) {
  std::ostringstream s;s<<"outside: "<<h.GetBinContent(0)<<" below, "<<h.GetBinContent(h.GetNbinsX()+1)<<" above";
  TLatex text;text.SetNDC();text.SetTextColor(color);text.SetTextSize(.023);text.DrawLatex(.54,y,s.str().c_str());
 }
 void draw(const fs::path&out,double beam)const {
  fs::create_directories(out/"figures");std::ostringstream title;title<<beam<<" GeV beam";
  for(bool full:{false,true}){
   TCanvas canvas("event_overview","",1100,650);canvas.SetBottomMargin(.13);canvas.SetLogy(full);
   TLegend leg(.48,.72,.90,.90);leg.SetTextFont(42);leg.SetTextSize(.027);leg.SetBorderSize(0);leg.SetHeader("Sample (below / above range)");
   std::array<std::unique_ptr<TH1D>,3> curves;double peak=0,smallest=1;
   const auto&source=full?fullTotals:totals;
   for(int s=0;s<3;++s){
    curves[s].reset(static_cast<TH1D*>(source[s]->Clone(("draw_"+samples[s]).c_str())));auto&h=*curves[s];
    h.Scale(1./h.GetEntries());h.SetLineColor(colors[s]);h.SetLineWidth(2);h.GetYaxis()->SetTitle("Probability / bin");h.SetTitle(title.str().c_str());
    peak=std::max(peak,h.GetMaximum());smallest=std::min(smallest,1./source[s]->GetEntries());
    std::ostringstream label;label<<samples[s]<<" ("<<source[s]->GetBinContent(0)<<" / "<<source[s]->GetBinContent(source[s]->GetNbinsX()+1)<<")";leg.AddEntry(&h,label.str().c_str(),"l");
   }
   for(int s=0;s<3;++s){curves[s]->SetMaximum(peak*(full?5:1.45));curves[s]->SetMinimum(full?smallest*.5:0);curves[s]->Draw(s?"HIST SAME":"HIST");}
   leg.Draw();
   canvas.SaveAs((out/"figures"/(full?"event_energy_full_range.png":"event_energy.png")).c_str());
   canvas.Clear();canvas.SetLogy(false);canvas.Divide(2,1);
   std::array<TH1D*,2> effects={full?fullChange.get():change.get(),full?fullRelative.get():relative.get()};
   for(int p=0;p<2;++p){canvas.cd(p+1);gPad->SetLeftMargin(.16);gPad->SetRightMargin(.04);gPad->SetBottomMargin(.13);gPad->SetLogy(full&&effects[p]->GetEntries()>0);effects[p]->GetYaxis()->SetTitleOffset(1.7);effects[p]->GetYaxis()->CenterTitle();effects[p]->Draw("HIST");flowLabel(*effects[p],.86);}
   canvas.SaveAs((out/"figures"/(full?"event_energy_change_full_range.png":"event_energy_change.png")).c_str());
  }
 }
};
} // namespace eventplots
