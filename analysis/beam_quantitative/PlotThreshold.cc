#include "ThresholdPlotCommon.hh"
#include <TObjString.h>

int main(int argc,char** argv) {
  try {
    check(argc==8,"PlotThreshold data_hist mc_hist plots comparison.root energy condition cut_keys");
    gROOT->SetBatch();TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetOptTitle(1);
    gStyle->SetTitleSize(.05,"XYZ");gStyle->SetLabelSize(.045,"XYZ");
    std::unique_ptr<TFile>d(TFile::Open(argv[1])),m(TFile::Open(argv[2]));
    check(d&&!d->IsZombie()&&m&&!m->IsZombie(),"Cannot open pre-SSA histograms");
    fs::path dir(argv[3]);fs::create_directories(dir);TFile out(argv[4],"RECREATE");check(!out.IsZombie(),"Cannot create threshold canvas ROOT");
    std::string heading=std::string(argv[5])+" GeV | "+argv[6];
    auto* status=dynamic_cast<TObjString*>(m->Get("audit/threshold_mode"));
    const bool fallback=status&&std::string(status->GetString().Data())=="fallback_detected";
    if(fallback)heading+=" | threshold fallback detected";
    {
      Page page("threshold_input_audit",heading+" | Digitization threshold input",
                "Production-log evidence | offline diagnostics cannot restore hits removed during acquisition",1,1);
      page.body.cd(1);TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.045);
      text.SetTextColor(fallback?kRed+1:kBlack);
      text.DrawLatex(.07,.80,fallback?"Threshold parameter loading FAILED in matching MC logs":"Threshold status: see matching production logs and saved JSON");
      text.SetTextColor(kBlack);text.SetTextSize(.035);
      text.DrawLatex(.07,.62,fallback?"MC continued with fallback thresholds; inspect fallback ADC values in JSON":"A supplied parameter ROOT is a provenance annotation, not proof of producer use");
      text.DrawLatex(.07,.46,"Pedestal and calibration selections still apply");
      text.DrawLatex(.07,.30,"Compare physical strips before SSA, then repeat with common offline MIP cuts");
      text.DrawLatex(.07,.14,"Parameter scans require the same original generation ROOT and unchanged data reference");
      page.save(out,dir);
    }
    std::vector<std::string> cuts;std::istringstream ss(argv[7]);std::string token;
    while(std::getline(ss,token,','))cuts.push_back(token);
    std::vector<std::unique_ptr<TH1>> keep;
    {
      Page page("strip_signal_spectrum",heading+" | Strip signal",
                "Same axis-valid events | signed signal after merging memory-cell duplicates | normalized strip distributions",2,1);
      for(int i=0;i<2;++i) {
        auto a=normalizedSignal(*d,"signal/all"),b=normalizedSignal(*m,"signal/all");
        page.body.cd(i+1);pairDraw(*a,*b,i?"Low-signal region;Strip signal [MIP];Fraction / bin":"Signal spectrum;Strip signal [MIP];Fraction / bin",true);
        a->GetXaxis()->SetRangeUser(i?-1.:-2.,i?3.:20.);
        if(!i){page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");}
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }
      page.save(out,dir);
    }
    {
      Page page("strip_signal_by_layer",heading+" | Strip signal by layer",
                "Signed physical-strip signal | same events in every panel | normalized separately by layer",5,6);
      for(int l=0;l<30;++l){
        auto a=normalizedSignal(*d,"signal/layer_"+std::to_string(l),l),b=normalizedSignal(*m,"signal/layer_"+std::to_string(l),l);
        page.body.cd(l+1);pairDraw(*a,*b,"Layer "+std::to_string(l)+";Strip signal [MIP];Fraction / bin",true);a->GetXaxis()->SetRangeUser(-1,5);
        if(!l){page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");}
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }page.save(out,dir);
    }
    {
      Page page("strip_hit_counts",heading+" | Counts before and after SSA",
                "Left: physical strips | Right: original SSA unique positions | no additional MIP cut | 95% event-bootstrap intervals",2,1);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;std::vector<std::unique_ptr<TH1D>> frames;
      for(int i=0;i<2;++i){
        auto a=curve(*d,i?"post/count":"pre/uncut/count",kBlack,20),b=curve(*m,i?"post/count":"pre/uncut/count",kRed+1,24);
        auto da=clone(*d,i?"post/count_mean":"pre/uncut/count_mean",false),ma=clone(*m,i?"post/count_mean":"pre/uncut/count_mean",false);
        auto h=std::make_unique<TH1D>(("counts_frame_"+std::to_string(i)).c_str(),"",30,-.5,29.5);page.body.cd(i+1);
        frame(*h,i?"After SSA;Layer;Unique positions / event":"Before SSA;Layer;Physical strips / event",0,std::max(upperBound(*a),upperBound(*b))*1.2+1);
        a->Draw("LP same");b->Draw("LP same");
        if(!i){page.legend.AddEntry(a.get(),"Experimental data","lp");page.legend.AddEntry(b.get(),"Simulation","lp");}
        frames.push_back(std::move(h));graphs.push_back(std::move(a));graphs.push_back(std::move(b));
      }page.save(out,dir);
    }
    {
      Page page("strip_ssa_multiplicity",heading+" | SSA positions per input strip",
                "Ratio of mean original SSA-position count to mean pre-SSA physical-strip count | no additional MIP cut",1,1);
      auto a=clone(*d,"ssa_multiplicity/mean",false),b=clone(*m,"ssa_multiplicity/mean",false);
      page.body.cd(1);pairDraw(*a,*b,";Layer;Mean SSA positions / mean input strips",false);
      page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");page.save(out,dir);
    }
    const int nx=2,ny=(cuts.size()+1)/2;
    for(const std::string quantity:{"count","rms"}) {
      Page page("strip_mip_cut_"+quantity,heading+" | Strip "+(quantity=="count"?"strip counts":"position widths")+" versus MIP cut",
                "Offline cut on merged physical-strip signal | event selection fixed | 95% event-bootstrap intervals | SSA is unchanged",nx,ny);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;std::vector<std::unique_ptr<TH1D>> frames;
      for(size_t i=0;i<cuts.size();++i) {
        const auto prefix="pre/"+cuts[i]+"/"+quantity;
        auto a=curve(*d,prefix,kBlack,20),b=curve(*m,prefix,kRed+1,24);
        auto da=clone(*d,prefix+"_mean",false),ma=clone(*m,prefix+"_mean",false);
        auto h=std::make_unique<TH1D>((quantity+"_frame_"+std::to_string(i)).c_str(),"",30,-.5,29.5);
        page.body.cd(i+1);frame(*h,label(cuts[i])+";Layer;"+(quantity=="count"?"Physical strips / event":"Fine-axis RMS [mm]"),0,std::max(upperBound(*a),upperBound(*b))*1.2+1);
        a->Draw("LP same");b->Draw("LP same");
        if(!i){page.legend.AddEntry(a.get(),"Experimental data","lp");page.legend.AddEntry(b.get(),"Simulation","lp");}
        frames.push_back(std::move(h));graphs.push_back(std::move(a));graphs.push_back(std::move(b));
      }page.save(out,dir);
    }
    {
      Page page("strip_mip_cut_difference",heading+" | Strip count difference versus MIP cut",
                "MC / data - 1 | same events for all cuts | 95% event-bootstrap intervals | green lines: configured count tolerance",nx,ny);
      std::vector<std::unique_ptr<TGraphAsymmErrors>> graphs;std::vector<std::unique_ptr<TH1D>> frames;
      std::vector<std::unique_ptr<TLine>> lines;
      const double tolerance=clone(*m,"comparison/count_tolerance",false)->GetBinContent(1);
      for(size_t i=0;i<cuts.size();++i) {
        const auto prefix="comparison/"+cuts[i]+"/count";
        auto a=clone(*m,prefix+"_difference",false),lo=clone(*m,prefix+"_low",false),hi=clone(*m,prefix+"_high",false);
        auto g=std::make_unique<TGraphAsymmErrors>();double low=-.3,high=.3;
        for(int l=1;l<=30;++l) {
          const double v=a->GetBinContent(l),dl=lo->GetBinContent(l),du=hi->GetBinContent(l);
          if(!std::isfinite(v)||!std::isfinite(dl)||!std::isfinite(du))continue;
          int j=g->GetN();g->SetPoint(j,l-1,v);g->SetPointError(j,0,0,std::max(0.,v-dl),std::max(0.,du-v));low=std::min(low,dl-.1);high=std::max(high,du+.1);
        }
        auto h=std::make_unique<TH1D>(("difference_frame_"+std::to_string(i)).c_str(),"",30,-.5,29.5);
        page.body.cd(i+1);frame(*h,label(cuts[i])+";Layer;MC / data - 1",low,high);
        g->SetLineColor(kRed+1);g->SetMarkerColor(kRed+1);g->SetMarkerStyle(20);g->Draw("LP same");
        for(double v:{-tolerance,0.,tolerance}) {auto line=std::make_unique<TLine>(-.5,v,29.5,v);line->SetLineColor(v==0?kGray+1:kGreen+2);line->SetLineStyle(2);line->Draw();lines.push_back(std::move(line));}
        if(!i)page.legend.AddEntry(g.get(),"Count difference, with 95% interval","lp");
        frames.push_back(std::move(h));graphs.push_back(std::move(g));
      }page.save(out,dir);
    }
    for(const auto& cut:cuts) {
      Page page("strip_profiles_"+cut,heading+" | Strip position profiles",label(cut)+" | 5.3 mm bins (one strip pitch) | X in odd layers, Y in even layers | normalized by layer",5,6);
      for(int l=0;l<30;++l) {
        auto key="pre/"+cut+"/layer_"+std::to_string(l)+"/u";
        auto a=clone(*d,key,true),b=clone(*m,key,true);page.body.cd(l+1);const std::string axis=l%2?"X":"Y";
        pairDraw(*a,*b,"Layer "+std::to_string(l)+" ("+axis+");"+axis+" [mm];Fraction / strip bin",false);
        if(!l){page.legend.AddEntry(a.get(),"Experimental data","l");page.legend.AddEntry(b.get(),"Simulation","l");}
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }page.save(out,dir);
    }
    for(int parity=0;parity<2;++parity)for(int sample=0;sample<2;++sample) {
      const std::string name="strip_occupancy_"+std::string(sample?"mc":"data")+"_"+(parity?"odd":"even");
      Page page(name,heading+" | "+(sample?"Simulation":"Experimental data")+" strip occupancy",
                "One bin per strip: 5.3 x 45.4 mm pitch | occupancy / event | same data/MC color scale per layer",4,4);
      for(int i=0;i<15;++i) {
        int l=2*i+parity;std::string key="occupancy/layer_"+std::to_string(l);
        auto a=clone(*d,key,false),b=clone(*m,key,false);auto* aa=dynamic_cast<TH2*>(a.get());auto* bb=dynamic_cast<TH2*>(b.get());check(aa&&bb,"Not a 2D occupancy map");
        page.body.cd(i+1);pad();gPad->SetRightMargin(.19);auto* h=sample?bb:aa;
        h->SetTitle(("Layer "+std::to_string(l)+";X [mm];Y [mm]").c_str());h->GetZaxis()->SetTitle("Occupancy / event");h->SetMinimum(0);h->SetMaximum(std::max(aa->GetMaximum(),bb->GetMaximum()));h->Draw("colz");
        keep.push_back(std::move(a));keep.push_back(std::move(b));
      }page.save(out,dir);
    }
    out.Close();return 0;
  }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
