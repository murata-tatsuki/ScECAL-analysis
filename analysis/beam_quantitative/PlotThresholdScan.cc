#include "ThresholdPlotCommon.hh"

struct Candidate {
  std::string name;double scale;
  std::unique_ptr<TFile> data,mc;
  std::map<std::string,std::array<double,3>> total;
};
inline double graphHigh(const TGraphAsymmErrors& g) {
  double value=0;for(int i=0;i<g.GetN();++i)value=std::max(value,g.GetPointY(i)+g.GetErrorYhigh(i));return value;
}
int main(int argc,char** argv) {
  try {
    check(argc==7,"PlotThresholdScan conditions.tsv plots comparison.root energy name cut_keys");
    gROOT->SetBatch();TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetOptTitle(1);
    gStyle->SetTitleSize(.05,"XYZ");gStyle->SetLabelSize(.045,"XYZ");
    std::vector<Candidate> candidates;std::ifstream manifest(argv[1]);check(bool(manifest),"Cannot open condition manifest");std::string line;
    while(std::getline(manifest,line)) {
      std::istringstream ss(line);std::string name,scale,dh,mh,table;
      check(bool(std::getline(ss,name,'\t'))&&bool(std::getline(ss,scale,'\t'))&&bool(std::getline(ss,dh,'\t'))&&bool(std::getline(ss,mh,'\t'))&&bool(std::getline(ss,table)),"Invalid condition manifest");
      Candidate c{name,std::stod(scale),std::unique_ptr<TFile>(TFile::Open(dh.c_str())),std::unique_ptr<TFile>(TFile::Open(mh.c_str())),{}};
      check(c.data&&!c.data->IsZombie()&&c.mc&&!c.mc->IsZombie(),"Cannot open diagnostic histograms");
      std::ifstream metrics(table);check(bool(metrics),"Cannot open diagnostic metric table");std::string row;std::getline(metrics,row);
      std::vector<std::string> columns;{std::istringstream in(row);while(std::getline(in,row,'\t'))columns.push_back(row);}
      while(std::getline(metrics,row)) {
        std::istringstream in(row);std::map<std::string,std::string> r;
        for(const auto& column:columns){std::string v;std::getline(in,v,'\t');r[column]=v;}
        if(r.at("layer")=="-1"&&r.at("metric")=="hits_relative")c.total[r.at("cut_key")]={std::stod(r.at("difference")),std::stod(r.at("ci_low")),std::stod(r.at("ci_high"))};
      }candidates.push_back(std::move(c));
    }
    check(candidates.size()>=2,"At least two threshold candidates are needed");
    std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b){return a.scale<b.scale;});
    fs::path dir(argv[2]);fs::create_directories(dir);TFile out(argv[3],"RECREATE");check(!out.IsZombie(),"Cannot create comparison ROOT");
    const std::string heading=std::string(argv[4])+" GeV | "+argv[5];std::vector<std::string> cuts;std::istringstream ss(argv[6]);std::string token;while(std::getline(ss,token,','))cuts.push_back(token);
    const std::array<int,8> colors={kRed+1,kBlue+1,kGreen+2,kOrange+7,kMagenta+1,kCyan+2,kViolet+1,kGray+2};
    auto legend=[](Page& p){p.legend.SetNColumns(3);p.legend.SetX1(.03);p.legend.SetX2(.97);p.legend.SetTextSize(.13);};
    auto title=[](const Candidate& c){std::ostringstream s;s<<c.name<<" (threshold x"<<c.scale<<")";return s.str();};
    {
      Page page("threshold_scan_hit_counts",heading+" | Threshold candidates: strip and SSA counts",
                "Same generation input and data reference | original axis-valid selection | 95% event-bootstrap intervals",2,1);legend(page);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;std::vector<std::unique_ptr<TH1D>> frames;
      for(int i=0;i<2;++i) {
        const std::string prefix=i?"post/count":"pre/uncut/count";
        auto data=curve(*candidates[0].data,prefix,kBlack,20);double high=graphHigh(*data);
        std::vector<std::unique_ptr<TGraphAsymmErrors>> mc;
        for(size_t j=0;j<candidates.size();++j){auto g=curve(*candidates[j].mc,prefix,colors[j],24+j%4);high=std::max(high,graphHigh(*g));mc.push_back(std::move(g));}
        auto h=std::make_unique<TH1D>(("scan_count_frame_"+std::to_string(i)).c_str(),"",30,-.5,29.5);page.body.cd(i+1);
        frame(*h,i?"After SSA;Layer;Unique positions / event":"Before SSA;Layer;Physical strips / event",0,high*1.2+1);
        data->Draw("LP same");if(!i)page.legend.AddEntry(data.get(),"Experimental data","lp");graphs.push_back(std::move(data));
        for(size_t j=0;j<mc.size();++j){mc[j]->Draw("LP same");if(!i)page.legend.AddEntry(mc[j].get(),title(candidates[j]).c_str(),"lp");graphs.push_back(std::move(mc[j]));}frames.push_back(std::move(h));
      }page.save(out,dir);
    }
    {
      Page page("threshold_scan_signal_spectrum",heading+" | Threshold candidates: strip signal",
                "Signed merged physical-strip signal | normalized distributions | same data reference",2,1);legend(page);
      std::vector<std::unique_ptr<TH1>> hist;
      for(int i=0;i<2;++i) {
        auto d=normalizedSignal(*candidates[0].data,"signal/all");std::vector<std::unique_ptr<TH1>> mc;double high=d->GetMaximum();
        for(size_t j=0;j<candidates.size();++j){auto h=normalizedSignal(*candidates[j].mc,"signal/all");high=std::max(high,h->GetMaximum());h->SetLineColor(colors[j]);h->SetLineWidth(2);h->SetLineStyle(2+j%3);mc.push_back(std::move(h));}
        page.body.cd(i+1);pad();gPad->SetLogy();d->SetTitle(i?"Low-signal region;Strip signal [MIP];Fraction / bin":"Signal spectrum;Strip signal [MIP];Fraction / bin");d->SetLineColor(kBlack);d->SetLineWidth(2);d->SetMaximum(high*2+1e-12);d->SetMinimum(std::max(1e-8,high*1e-5));d->GetXaxis()->SetRangeUser(i?-1.:-2.,i?3.:20.);d->Draw("hist");
        if(!i)page.legend.AddEntry(d.get(),"Experimental data","l");
        hist.push_back(std::move(d));
        for(size_t j=0;j<mc.size();++j){mc[j]->Draw("hist same");if(!i)page.legend.AddEntry(mc[j].get(),title(candidates[j]).c_str(),"l");hist.push_back(std::move(mc[j]));}
      }page.save(out,dir);
    }
    {
      Page page("threshold_scan_mip_cut_difference",heading+" | Threshold candidates: count difference after common MIP cuts",
                "MC / data - 1 | offline cut before SSA | fixed events within each sample | 95% event-bootstrap intervals",2,(cuts.size()+1)/2);legend(page);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;std::vector<std::unique_ptr<TH1D>> frames;std::vector<std::unique_ptr<TLine>> lines;
      for(size_t i=0;i<cuts.size();++i) {
        std::vector<std::unique_ptr<TGraphAsymmErrors>> mc;double low=-.3,high=.3;
        for(size_t j=0;j<candidates.size();++j) {
          auto g=curve(*candidates[j].mc,"comparison/"+cuts[i]+"/count",colors[j],24+j%4,"_difference");
          for(int k=0;k<g->GetN();++k){low=std::min(low,g->GetPointY(k)-g->GetErrorYlow(k)-.1);high=std::max(high,g->GetPointY(k)+g->GetErrorYhigh(k)+.1);}mc.push_back(std::move(g));
        }
        auto h=std::make_unique<TH1D>(("scan_difference_frame_"+std::to_string(i)).c_str(),"",30,-.5,29.5);page.body.cd(i+1);frame(*h,label(cuts[i])+";Layer;MC / data - 1",low,high);
        auto zero=std::make_unique<TLine>(-.5,0,29.5,0);zero->SetLineStyle(2);zero->Draw();lines.push_back(std::move(zero));
        for(size_t j=0;j<mc.size();++j){mc[j]->Draw("LP same");if(!i)page.legend.AddEntry(mc[j].get(),title(candidates[j]).c_str(),"lp");graphs.push_back(std::move(mc[j]));}frames.push_back(std::move(h));
      }page.save(out,dir);
    }
    {
      Page page("threshold_scan_response",heading+" | Total strip-count response to threshold",
                "All layers summed | MC / data - 1 | x: annotated threshold-mean factor | each curve uses a common offline MIP cut",1,1);legend(page);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;double low=-.3,high=.3;
      for(size_t j=0;j<cuts.size();++j) {
        auto g=std::make_unique<TGraphAsymmErrors>();g->SetLineColor(colors[j]);g->SetMarkerColor(colors[j]);g->SetMarkerStyle(20+j%4);g->SetLineWidth(2);
        for(const auto& c:candidates) {
          const auto& v=c.total.at(cuts[j]);if(!std::isfinite(v[0])||!std::isfinite(v[1])||!std::isfinite(v[2]))continue;
          int n=g->GetN();g->SetPoint(n,c.scale,v[0]);g->SetPointError(n,0,0,std::max(0.,v[0]-v[1]),std::max(0.,v[2]-v[0]));low=std::min(low,v[1]-.1);high=std::max(high,v[2]+.1);
        }page.legend.AddEntry(g.get(),label(cuts[j]).c_str(),"lp");graphs.push_back(std::move(g));
      }
      double xmin=candidates.front().scale,xmax=candidates.back().scale;double margin=std::max(.05,(xmax-xmin)*.1);
      TH1D h("scan_response_frame","",100,xmin-margin,xmax+margin);page.body.cd(1);frame(h,";Threshold-mean scale;Total strip count: MC / data - 1",low,high);
      TLine zero(xmin-margin,0,xmax+margin,0);zero.SetLineStyle(2);zero.Draw();for(auto& g:graphs)g->Draw("LP same");page.save(out,dir);
    }
    out.Close();return 0;
  }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
