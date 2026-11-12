#include "../beam_quantitative/PlotCommon.hh"
#include <TGraphAsymmErrors.h>
#include <TObjString.h>
#include <numeric>
#include <set>
#include <TTree.h>
#include "../calibration_parameters/EBUdecode.h"

struct StagePage:Page {
  StagePage(const std::string& name,const std::string& title,const std::string& subtitle,int nx,int ny):
    Page(name,title,subtitle,nx,ny) {
    legend.SetX1(.025);legend.SetX2(.975);
    legend.SetY1(.015);legend.SetY2(.34);legend.SetTextSize(.13);
  }
};

void panelLabel(const std::string& title) {
  const std::string caption=title.substr(0,title.find(';'));
  if(caption.empty())return;
  TLatex text;text.SetNDC();text.SetTextFont(42);
  text.SetTextSize(std::min(.05,1.35/caption.size()));text.DrawLatex(.15,.94,caption.c_str());
}

struct Curve { std::string key,label;int color,marker,style; };
const std::vector<Curve> stagesPlot={
 {"truth","MC truth: deposited energy > 0",kGray+2,25,3},
 {"digitized","MC after digitization (threshold applied)",kBlue+1,24,2},
 {"calibrated","MC after calibration",kRed+1,20,1},
 {"data","Experimental data after calibration",kBlack,21,1}};

std::unique_ptr<TGraphAsymmErrors> graph(TFile& f,const std::string& key,int color,int marker,int style,const std::string& prefix="curves/",bool survival=false) {
  auto a=clone(f,prefix+key+"_mean",false),lo=clone(f,prefix+key+"_low",false),hi=clone(f,prefix+key+"_high",false);
  std::unique_ptr<TH1> support;if(survival)support=clone(f,prefix+"support",false);
  auto g=std::make_unique<TGraphAsymmErrors>();
  for(int i=1;i<=a->GetNbinsX();++i){
    if(support && support->GetBinContent(i)==0)continue;
    const double v=a->GetBinContent(i);if(!std::isfinite(v))continue;
    const int j=g->GetN();g->SetPoint(j,a->GetBinCenter(i),v);
    g->SetPointError(j,0,0,std::max(0.,v-lo->GetBinContent(i)),std::max(0.,hi->GetBinContent(i)-v));
  }
  g->SetLineColor(color);g->SetMarkerColor(color);g->SetMarkerStyle(marker);g->SetMarkerSize(.7);g->SetLineStyle(style);g->SetLineWidth(2);return g;
}

int main(int argc,char** argv) {
 try {
  check(argc==6 || argc==7,"Plot histogram.root plots_dir plots.root energy condition [threshold_cut]");
  gROOT->SetBatch();TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetOptTitle(1);gStyle->SetTitleSize(.05,"XYZ");gStyle->SetLabelSize(.045,"XYZ");
  TFile f(argv[1]);check(!f.IsZombie(),"Cannot read stage histograms");fs::path dir(argv[2]);fs::create_directories(dir);
  TFile out(argv[3],"RECREATE");check(!out.IsZombie(),"Cannot create canvas ROOT");
  const std::string cut=argc==7?argv[6]:"allChannel";
  const std::string heading=std::string(argv[4])+" GeV | "+argv[5]+(cut=="allChannel"?"":" | channel cut: "+cut);
  std::vector<std::unique_ptr<TH1>> hs;std::vector<std::unique_ptr<TGraphAsymmErrors>> gs;std::vector<std::unique_ptr<TLine>> lines;
  int serial=0;
  if(cut!="allChannel") {
    {
      StagePage p("threshold_channel_coverage",heading+" | Channels retained by threshold fit quality",
        "Same physical channel mask for data and all MC stages | original events and entrance axes stay fixed",1,1);
      p.body.cd(1);auto all=clone(f,"quality/available_counts",false),selected=clone(f,"quality/selected_counts",false);
      pairDraw(*all,*selected,";Layer;Selected / connected channels",false);
      selected->SetLineColor(kBlue+1);p.legend.AddEntry(all.get(),"Connected channels (210 per layer)","l");
      p.legend.AddEntry(selected.get(),"Channels passing this quality cut","l");p.save(out,dir);
      hs.push_back(std::move(all));hs.push_back(std::move(selected));
    }
    auto* tree=dynamic_cast<TTree*>(f.Get("quality/selected_channels"));check(tree,"Missing selected-channel table");
    int id=0;check(tree->SetBranchAddress("CellID",&id)>=0,"Invalid channel table");std::set<int> selected;
    for(Long64_t i=0;i<tree->GetEntries();++i){check(tree->GetEntry(i)>0,"Cannot read channel mask");selected.insert(id);}
    tree->ResetBranchAddresses();
    for(int parity=0;parity<2;++parity) {
      StagePage p("threshold_channel_mask_"+std::string(parity?"odd":"even"),heading+" | Selected strip locations",
        "0 = excluded, 1 = retained | one bin per strip pitch (5.3 x 45.4 mm) | no event occupancy weighting",4,4);
      for(int j=0;j<15;++j){const int l=2*j+parity;p.body.cd(j+1);pad();gPad->SetRightMargin(.17);
        const int nx=parity?42:5,ny=parity?5:42;const double xmax=parity?111.3:113.5,ymax=parity?113.5:111.3;
        auto h=std::make_unique<TH2D>(("quality_mask_"+std::to_string(l)).c_str(),("Layer "+std::to_string(l)+";X [mm];Y [mm]").c_str(),nx,-xmax,xmax,ny,-ymax,ymax);
        for(int c=0;c<6;++c)for(int ch=0;ch<(c==5?30:36);++ch){const int cell=l*100000+c*10000+ch;
          if(selected.count(cell)){auto* pos=EBUdecode(l,c,ch);h->Fill(pos[0],pos[1]);}}
        h->SetMinimum(0);h->SetMaximum(1);h->Draw("colz");hs.push_back(std::move(h));
      }p.save(out,dir);
    }
  }
  auto drawCurves=[&](Page& page,int panel,const std::vector<Curve>& curves,const std::string& title,double ymin,double ymax,bool addLegend) {
    page.body.cd(panel);pad();std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;
    for(const auto& c:curves) {
      auto g=graph(f,c.key,c.color,c.marker,c.style);
      for(int j=0;j<g->GetN();++j){double x,y;g->GetPoint(j,x,y);if(ymax<0)ymin=std::min(ymin,y-g->GetErrorYlow(j));}
      graphs.push_back(std::move(g));
    }
    if(ymax<0){ymax=0;for(const auto& g:graphs)for(int j=0;j<g->GetN();++j){double x,y;g->GetPoint(j,x,y);ymax=std::max(ymax,y+g->GetErrorYhigh(j));}ymax=std::max(1.,ymax*1.18);}
    auto h=std::make_unique<TH1D>(("frame_"+std::to_string(serial++)).c_str(),title.c_str(),30,-.5,29.5);
    h->SetMinimum(ymin);h->SetMaximum(ymax);h->GetYaxis()->SetTitleOffset(1.5);h->Draw("axis");panelLabel(title);
    for(size_t i=0;i<graphs.size();++i){graphs[i]->Draw("LP same");if(addLegend)page.legend.AddEntry(graphs[i].get(),curves[i].label.c_str(),"lp");gs.push_back(std::move(graphs[i]));}
    hs.push_back(std::move(h));
  };
  {
    StagePage p("stage_strip_counts",heading+" | Physical strip counts at each MC stage",
           "Same accepted events in all stages | one count per physical CellID | error bars: 95% event-bootstrap intervals",2,1);
    drawCurves(p,1,stagesPlot,"All layers;Layer;Physical strips / event",0,-1,true);
    drawCurves(p,2,stagesPlot,"Deep layers 20-28;Layer;Physical strips / event",0,-1,false);
    hs.back()->GetXaxis()->SetRangeUser(19.5,28.5);p.save(out,dir);
  }
  {
    StagePage p("calibrated_strip_excess",heading+" | MC / data strip-count difference",
      "Same calibrated-event cohort for every offline MIP cut | zero line: MC and data agree | 95% event-bootstrap intervals",1,1);
    std::vector<Curve> cuts={{"mc_data_difference","No additional signal cut",kRed+1,20,1},
      {"mc_data_difference_mip_0_5","Signal >= 0.5 MIP",kBlue+1,24,2},
      {"mc_data_difference_mip_1_0","Signal >= 1 MIP",kGreen+2,25,3},
      {"mc_data_difference_mip_2_0","Signal >= 2 MIP",kMagenta+1,26,4}};
    drawCurves(p,1,cuts,";Layer;MC / data - 1",-.5,cut=="allChannel"?1.:-1.,true);
    auto line=std::make_unique<TLine>(-.5,0,29.5,0);line->SetLineColor(kGray+1);line->SetLineStyle(2);line->Draw();lines.push_back(std::move(line));p.save(out,dir);
  }
  {
    StagePage p("stage_strip_retention",heading+" | Which strips survive each stage?",
      "Matched physical CellIDs | fractions use stage-specific input counts | overlapping blue/red curves indicate equal values",2,1);
    drawCurves(p,1,{{"retention_digitization","Truth -> digitization",kBlue+1,24,2},
      {"retention_calibration","Digitization -> calibration",kGreen+2,25,3},
      {"retention_total","Truth -> calibration",kRed+1,20,1}},"Strip survival fraction;Layer;Matched output strips / input strips",0,1.08,true);
    drawCurves(p,2,{{"lost_digitization","Removed by digitization",kBlue+1,24,2},
      {"lost_calibration","Removed by calibration",kRed+1,20,1},
      {"added_digitization","Digitized strips absent in truth",kMagenta+1,26,4},
      {"added_calibration","Calibrated strips absent in digitization",kOrange+7,27,3}},
      "Removed / unexpected added strips;Layer;Physical strips / event",0,-1,false);
    // Panel-specific legend entries fit in the header, away from the curves.
    for(size_t i=gs.size()-4;i<gs.size();++i)p.legend.AddEntry(gs[i].get(),i==gs.size()-4?"Removed by digitization":i==gs.size()-3?"Removed by calibration":i==gs.size()-2?"Added at digitization":"Added at calibration","lp");
    p.legend.SetNColumns(3);p.legend.SetTextSize(.10);p.save(out,dir);
  }
  {
    StagePage p("stage_energy_response",heading+" | Energy selection and calibration response",
      "Truth energy kept: deposited energy on surviving IDs | right: same matched truth-positive strips | ratios of energy sums",2,1);
    drawCurves(p,1,{{"truth_energy","MC truth: all deposited energy",kGray+2,25,3},
      {"truth_energy_digitized","MC truth energy on digitized strips",kBlue+1,24,2},
      {"calib_energy","MC reconstructed positive energy",kRed+1,20,1},
      {"data_energy","Data reconstructed positive energy",kBlack,21,1}},"Layer energy;Layer;Energy / event [MeV]",0,-1,true);
    drawCurves(p,2,{{"calibration_response","MC reconstructed / deposited energy",kRed+1,20,1}},
      "Response on matched strips;Layer;Reconstructed energy / deposited energy",0,-1,false);
    auto line=std::make_unique<TLine>(-.5,1,29.5,1);line->SetLineStyle(2);line->Draw();lines.push_back(std::move(line));p.save(out,dir);
  }
  {
    StagePage p("stage_strip_width",heading+" | Strip distribution width by layer",
      "X in odd layers; Y in even layers | entrance-centroid position subtraction, same axis for every MC stage | count-weighted RMS",1,1);
    auto curves=stagesPlot;for(auto& c:curves)c.key="width_"+c.key;
    drawCurves(p,1,curves,";Layer;Fine-axis RMS [mm]",0,-1,true);p.save(out,dir);
  }
  // Pooled strip diagnostics are separated by SiPM cell size. The saved
  // sample remains layers 0..29 (15 um layers 30/31 are outside this analysis).
  const std::array<std::string,2> sipmKeys={"10um","15um"};
  const std::array<std::string,2> sipmLabels={"10 #mum (layers 4-27)","15 #mum (layers 0-3, 28-29)"};
  {
    StagePage p("truth_signal_survival",heading+" | Strip survival by SiPM size and truth energy",
      "Top: 10 um | bottom: 15 um | same truth denominator for both stages | 95% event-bootstrap intervals",4,2);
    const std::vector<std::pair<std::string,std::string>> groups={
      {"10um/all","10 #mum: all layers 4-27"}, {"10um/layers_4_9","10 #mum: layers 4-9"},
      {"10um/layers_10_19","10 #mum: layers 10-19"}, {"10um/layers_20_27","10 #mum: layers 20-27"},
      {"15um/all","15 #mum: all layers 0-3, 28-29"}, {"15um/layers_0_3","15 #mum: layers 0-3"},
      {"15um/layers_28_29","15 #mum: layers 28-29"}, {"15um/layer_29","15 #mum: layer 29"}};
    for(size_t i=0;i<groups.size();++i){p.body.cd(i+1);pad();gPad->SetLogx();
      auto h=std::make_unique<TH1D>(("survival_frame_"+std::to_string(i)).c_str(),(groups[i].second+";Deposited energy / nominal MIP;Fraction of truth strips surviving").c_str(),100,.001,100);
      h->SetMinimum(0);h->SetMaximum(1.08);h->GetYaxis()->SetTitleOffset(1.5);h->Draw("axis");panelLabel(groups[i].second);
      const auto split=groups[i].first.find('/');
      const std::string prefix="sipm/"+groups[i].first.substr(0,split)+"/survival/"+groups[i].first.substr(split+1)+"/";
      for(const auto& stage:{stagesPlot[1],stagesPlot[2]}){auto g=graph(f,stage.key,stage.color,stage.marker,stage.style,prefix,true);g->Draw("LP same");
        if(!i){p.legend.AddEntry(g.get(),stage.label.c_str(),"lp");}gs.push_back(std::move(g));}hs.push_back(std::move(h));
    }p.save(out,dir);
  }
  {
    StagePage p("stage_signal_spectra",heading+" | Strip signals by SiPM size",
      "Top: 10 um | bottom: 15 um | left: truth deposits, right: calibrated signals | each curve normalized within its SiPM type",2,2);
    for(int panel=0;panel<2;++panel){
      const std::vector<Curve> curves=panel?std::vector<Curve>{{"calibrated","MC reconstructed signal",kRed+1,20,1},{"data","Data reconstructed signal",kBlack,21,1}}:
        std::vector<Curve>{{"truth","All MC truth strips",kGray+2,25,3},{"truth_digitized","Truth strips surviving digitization",kBlue+1,24,2},{"truth_calibrated","Truth strips surviving calibration",kRed+1,20,1}};
      std::array<std::vector<std::unique_ptr<TH1>>,2> v;double maximum=1e-6;
      for(int row=0;row<2;++row)for(const auto& c:curves){auto h=clone(f,"sipm/"+sipmKeys[row]+"/signal/"+c.key,true);
        h->SetLineColor(c.color);h->SetLineStyle(c.style);h->SetLineWidth(2);maximum=std::max(maximum,h->GetMaximum());v[row].push_back(std::move(h));}
      for(int row=0;row<2;++row){p.body.cd(2*row+panel+1);pad();gPad->SetLogx();gPad->SetLogy();
        for(size_t i=0;i<v[row].size();++i){auto& h=v[row][i];
          const std::string title=sipmLabels[row]+(panel?": calibrated signals;Reconstructed signal [MIP];Fraction / variable-width bin":": truth deposits;Deposited energy / nominal MIP;Fraction / variable-width bin");
          h->SetTitle(title.c_str());h->GetXaxis()->SetRangeUser(.001,200.);h->GetYaxis()->SetTitleOffset(1.5);
          h->SetMinimum(1e-6);h->SetMaximum(maximum*2);h->Draw(i?"hist same":"hist");
          if(!row){p.legend.AddEntry(h.get(),curves[i].label.c_str(),"l");}hs.push_back(std::move(h));}
      }
    }p.legend.SetNColumns(3);p.legend.SetTextSize(.10);p.save(out,dir);
  }
  {
    StagePage p("digitized_adc",heading+" | Saved ADC by SiPM size",
      "Top: 10 um | bottom: 15 um | left: HG, right: LG | threshold-selected, pedestal included | single-row strips only",2,2);
    for(int col=0;col<2;++col){std::array<std::unique_ptr<TH1>,2> v;double maximum=1e-6;
      for(int row=0;row<2;++row){v[row]=clone(f,"sipm/"+sipmKeys[row]+(col?"/signal/lg_adc":"/signal/hg_adc"),true);maximum=std::max(maximum,v[row]->GetMaximum());}
      for(int row=0;row<2;++row){p.body.cd(2*row+col+1);pad();gPad->SetLogy();auto& h=v[row];
        h->SetTitle((sipmLabels[row]+(col?": low gain;Raw LG ADC;Fraction / 20 ADC bin":": high gain;Raw HG ADC;Fraction / 20 ADC bin")).c_str());
        h->SetLineColor(kBlue+1);h->SetLineWidth(2);h->SetMinimum(1e-6);h->SetMaximum(maximum*2);h->GetYaxis()->SetTitleOffset(1.5);h->Draw("hist");
        if(!row && !col){p.legend.AddEntry(h.get(),"MC saved ADC (each SiPM type normalized separately)","l");}hs.push_back(std::move(h));}
    }p.save(out,dir);
  }
  {
    StagePage p("truth_calibration_correlation",heading+" | Matched strip response by SiPM size",
      "Left: 10 um | right: 15 um | same physical strip, truth > 0 | shared color scale: strip entries | red line: equal energy",2,1);
    std::array<std::unique_ptr<TH1>,2> v;double maximum=1;
    for(int col=0;col<2;++col){v[col]=clone(f,"sipm/"+sipmKeys[col]+"/response/truth_calib",false);maximum=std::max(maximum,v[col]->GetMaximum());}
    for(int col=0;col<2;++col){p.body.cd(col+1);pad();gPad->SetRightMargin(.17);gPad->SetLogz();auto& h=v[col];
      h->SetTitle((sipmLabels[col]+";Deposited energy / nominal MIP;Reconstructed signal [MIP]").c_str());
      h->SetMinimum(1);h->SetMaximum(std::max(2.,maximum));h->Draw("colz");
      auto line=std::make_unique<TLine>(0,0,50,50);line->SetLineColor(kRed+1);line->SetLineWidth(2);line->Draw();
      if(!col)p.legend.AddEntry(line.get(),"Equal deposited and reconstructed energy","l");
      hs.push_back(std::move(h));lines.push_back(std::move(line));}
    p.save(out,dir);
  }
  for(const std::string frame:{"lab","aligned"}){
    StagePage p("stage_strip_profiles_"+frame,heading+" | Strip position profiles",
      frame=="lab"?"Detector coordinates | X in odd layers; Y in even layers | contiguous 5.3 mm strip-pitch bins | each layer/stage normalized separately":
      "Entrance-centroid subtracted | same axis in all MC stages | X in odd layers; Y in even layers | each layer/stage normalized separately",5,6);
    for(int l=0;l<30;++l){p.body.cd(l+1);pad();std::vector<std::unique_ptr<TH1>> v;double maximum=0;
      for(const auto& c:stagesPlot){auto h=clone(f,"profiles/"+frame+"/"+c.key+"/layer_"+std::to_string(l),true);
        h->SetLineColor(c.color);h->SetLineStyle(c.style);h->SetLineWidth(2);maximum=std::max(maximum,h->GetMaximum());v.push_back(std::move(h));}
      for(size_t i=0;i<v.size();++i){auto& h=v[i];const std::string axis=l%2?"X":"Y";
        h->SetTitle(("Layer "+std::to_string(l)+" ("+axis+");"+axis+(frame=="lab"?" [mm]":" - entrance centroid [mm]")+";Fraction / 5.3 mm bin").c_str());
        h->GetXaxis()->SetRangeUser(-111.3,111.3);h->SetMinimum(0);h->SetMaximum(maximum*1.2+1e-6);h->GetYaxis()->SetTitleOffset(1.5);h->Draw(i?"hist same":"hist");
        if(!l){p.legend.AddEntry(h.get(),stagesPlot[i].label.c_str(),"l");}hs.push_back(std::move(h));}
    }p.save(out,dir);
  }
  std::cout<<"Wrote stage plots to "<<dir<<'\n';
 }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
