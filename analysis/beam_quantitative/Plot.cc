#include "PlotCommon.hh"

int main(int argc,char** argv) {
  try {
    check(argc==8,"Plot data_hist mc_hist metrics.tsv plots comparison.root energy condition");
    gROOT->SetBatch();TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetTitleSize(.05,"XYZ");gStyle->SetLabelSize(.045,"XYZ");
    std::unique_ptr<TFile>d(TFile::Open(argv[1])),m(TFile::Open(argv[2]));
    check(d&&!d->IsZombie()&&m&&!m->IsZombie(),"Cannot open evaluation histograms");
    const fs::path dir(argv[4]);fs::create_directories(dir);
    TFile out(argv[5],"RECREATE");check(!out.IsZombie(),"Cannot create comparison ROOT");
    const std::string heading=std::string(argv[6])+" GeV | "+argv[7];
    std::vector<std::unique_ptr<TH1>> keep;
    for(bool log:{false,true}) {
      Page page(log?"entrance_positions_log":"entrance_positions",heading+" | Entrance centroid, one point per event",
                "Before SSA, layers 0-3 by default; common axis-valid selection; normalized distributions",2,1);
      for(int axis=0;axis<2;++axis) {
        const std::string coordinate=axis?"y":"x";
        auto a=clone(*d,"entrance/"+coordinate,true),b=clone(*m,"entrance/"+coordinate,true);
        page.body.cd(axis+1);pairDraw(*a,*b,coordinate+" direction;Entrance "+coordinate+" [mm];Fraction / bin",log);
        if(!axis){page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");}
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }
      page.save(out,dir);
    }
    {
      Page page("entrance_xy",heading+" | Entrance centroid, one point per event","Data and simulation share binning and color scale; flow is retained in ROOT and TV metric",2,1);
      auto a=clone(*d,"entrance/xy",false),b=clone(*m,"entrance/xy",false);
      auto* aa=dynamic_cast<TH2*>(a.get());auto* bb=dynamic_cast<TH2*>(b.get());check(aa&&bb,"Not a 2D histogram");
      if(aa->Integral()>0) aa->Scale(1./aa->Integral());
      if(bb->Integral()>0) bb->Scale(1./bb->Integral());
      const double maximum=std::max(aa->GetMaximum(),bb->GetMaximum());
      int i=0;for(auto* h:{aa,bb}) {page.body.cd(++i);pad();gPad->SetRightMargin(.18);h->SetTitle(i==1?"Experimental data;X [mm];Y [mm]":"Simulation;X [mm];Y [mm]");h->SetMinimum(0);h->SetMaximum(maximum);h->Draw("colz");}
      page.save(out,dir);keep.push_back(std::move(a));keep.push_back(std::move(b));
    }
    for(const std::string frame:{"lab","aligned"})for(const std::string weight:{"count","positive_energy"})for(bool log:{false,true}) {
      const std::string name="profiles_"+frame+"_"+weight+(log?"_log":"");
      Page page(name,heading+" | SSA hit profiles",
                (frame=="lab"?"Detector coordinates":"Relative to each event's entrance centroid")+std::string(" | ")+
                (weight=="count"?"Unique hit-position counts":"Positive energy after merging identical positions")+" | 5.3 mm bins | normalized by layer",5,6);
      for(int l=0;l<30;++l) {
        const std::string key=frame+"/"+weight+"/layer_"+std::to_string(l)+"/u";
        auto a=clone(*d,key,true),b=clone(*m,key,true);
        page.body.cd(l+1);const std::string axis=l%2?"X":"Y";
        pairDraw(*a,*b,"Layer "+std::to_string(l)+" ("+axis+");"+axis+" [mm];Fraction / 5.3 mm bin",log);
        if(frame=="aligned")a->GetXaxis()->SetRangeUser(-115,115);
        if(!l){page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");}
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }
      page.save(out,dir);
    }
    for(const std::string weight:{"count","positive_energy"})for(int parity=0;parity<2;++parity)for(int sample=0;sample<2;++sample) {
      const std::string name="hit_map_"+std::string(sample?"mc":"data")+"_"+(parity?"odd":"even")+"_"+weight;
      Page page(name,heading+" | "+(sample?"Simulation":"Experimental data")+" SSA hit maps",
                "Detector coordinates | "+std::string(weight=="count"?"Unique positions / event":"Positive signal energy [MeV] / event")+" | common data/MC color scale per layer",4,4);
      for(int i=0;i<15;++i) {
        int l=2*i+parity;const std::string key="lab/"+weight+"/layer_"+std::to_string(l)+"/xy";
        auto a=clone(*d,key,false),b=clone(*m,key,false);
        auto* aa=dynamic_cast<TH2*>(a.get());auto* bb=dynamic_cast<TH2*>(b.get());
        check(aa&&bb,"Not a 2D map");page.body.cd(i+1);pad();gPad->SetRightMargin(.19);
        auto* h=sample?bb:aa;h->SetTitle(("Layer "+std::to_string(l)+";X [mm];Y [mm]").c_str());
        h->SetMinimum(0);h->SetMaximum(std::max(aa->GetMaximum(),bb->GetMaximum()));h->Draw("colz");
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }
      page.save(out,dir);
    }
    std::ifstream table(argv[3]);check(bool(table),"Cannot open metrics TSV");
    std::string line;std::getline(table,line);std::vector<std::string>columns;
    {std::istringstream ss(line);while(std::getline(ss,line,'\t'))columns.push_back(line);}
    std::vector<std::map<std::string,std::string>> rows;
    while(std::getline(table,line)) {
      std::istringstream ss(line);std::map<std::string,std::string> r;
      for(const auto& c:columns){std::string v;std::getline(ss,v,'\t');r[c]=v;}rows.push_back(r);
    }
    {
      Page page("assessment",heading+" | Quantitative agreement", "Score = 95% CI bound / tolerance | <= 1 needs sufficient statistics | verdicts in TSV/JSON",2,1);
      std::array<double,30> lab{},aligned{};std::array<bool,30> labKnown{},alignedKnown{};
      std::vector<std::pair<std::string,double>> axis;
      for(const auto& r:rows) {
        double v=std::stod(r.at("score"));
        if(r.at("group")=="entrance"){
          const std::map<std::string,std::string> names={{"mean","Center"},{"rms_relative","Width"},{"q90_relative","90% width"},{"tail_fraction","Tail"},{"correlation","Correlation"},{"tv","2D shape"}};
          axis.push_back({names.at(r.at("metric"))+" "+r.at("axis"),v});
        }
        else {
          int l=std::stoi(r.at("layer"));if(!std::isfinite(v))continue;
          bool relative=r.at("frame")=="aligned";auto& a=relative?aligned:lab;auto& valid=relative?alignedKnown:labKnown;
          a[l]=std::max(a[l],v);valid[l]=true;
        }
      }
      page.body.cd(1);pad();TH1D entrance("entrance_scores","Entrance-position metrics;;Score",axis.size(),0,axis.size());
      double maximum=1.5;
      for(size_t i=0;i<axis.size();++i){entrance.GetXaxis()->SetBinLabel(i+1,axis[i].first.c_str());if(std::isfinite(axis[i].second)){entrance.SetBinContent(i+1,axis[i].second);maximum=std::max(maximum,axis[i].second*1.2);}}
      entrance.SetMaximum(maximum);entrance.SetMinimum(0);entrance.GetYaxis()->SetTitleOffset(1.5);entrance.GetXaxis()->LabelsOption("v");gPad->SetBottomMargin(.29);entrance.Draw("hist");
      TLine limit1(0,1,axis.size(),1);limit1.SetLineColor(kGreen+2);limit1.SetLineStyle(2);limit1.Draw();
      page.body.cd(2);pad();TH1D frame("layer_score_frame","Worst score per layer;Layer;Score",30,-.5,29.5);
      double mx=1.5;for(int l=0;l<30;++l)mx=std::max({mx,lab[l]*1.2,aligned[l]*1.2});frame.SetMinimum(0);frame.SetMaximum(mx);frame.GetXaxis()->SetTitleOffset(1.05);frame.GetYaxis()->SetTitleOffset(1.5);frame.Draw();
      TGraph gl,ga;for(int l=0;l<30;++l){if(labKnown[l])gl.SetPoint(gl.GetN(),l,lab[l]);if(alignedKnown[l])ga.SetPoint(ga.GetN(),l,aligned[l]);}
      gl.SetLineColor(kBlue+1);gl.SetMarkerColor(kBlue+1);gl.SetMarkerStyle(20);ga.SetLineColor(kOrange+7);ga.SetMarkerColor(kOrange+7);ga.SetMarkerStyle(24);
      gl.Draw("LP same");ga.Draw("LP same");TLine limit2(-.5,1,29.5,1);limit2.SetLineColor(kGreen+2);limit2.SetLineStyle(2);limit2.Draw();
      page.legend.AddEntry(&gl,"Detector coordinates","lp");page.legend.AddEntry(&ga,"Entrance-centered (diagnostic)","lp");
      page.save(out,dir);
    }
    out.Close();return 0;
  } catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
