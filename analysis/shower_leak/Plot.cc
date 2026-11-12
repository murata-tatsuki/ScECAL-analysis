#include "Common.hh"
#include "IntervalPlots.hh"
#include "TrackSelection.hh"
#include "PlotStyle.hh"
#include <TKey.h>
using Events=std::vector<Event>;
void trackAudit(const fs::path& dest,const Events& all,const std::vector<TrackTag>& tags,const TrackConfig& cfg,int cog,const fs::path& source) {
 TFile f((dest/"track_selection.root").c_str(),"RECREATE");require(!f.IsZombie(),"Cannot write track audit");
 TTree tree("selection","Track diagnostics for all sampled events; identity refers to source events.root");
 Event id;TrackTag t;bool before=false,after=false;
 tree.Branch("file_index",&id.file_index);tree.Branch("entry",&id.entry);tree.Branch("run",&id.run);tree.Branch("time",&id.time);tree.Branch("event",&id.event);tree.Branch("energy",&id.energy);tree.Branch("hits",&id.hits);tree.Branch("tail_rich",&id.tail_rich);tree.Branch("cog_valid",&id.cog_valid);
 tree.Branch("active_layers",&t.active);tree.Branch("active_front5",&t.front);tree.Branch("active_back5",&t.back);tree.Branch("nx",&t.nx);tree.Branch("ny",&t.ny);tree.Branch("hits_per_active_layer",&t.occupancy);tree.Branch("rms_x_mm",&t.rmsX);tree.Branch("rms_y_mm",&t.rmsY);
 tree.Branch("low",&t.low);tree.Branch("sparse",&t.sparse);tree.Branch("spanning",&t.spanning);tree.Branch("straight",&t.straight);tree.Branch("old_track_tag",&t.oldTrack);tree.Branch("low_signal_tag",&t.lowSignal);tree.Branch("track_tag",&t.tag);tree.Branch("before_track",&before);tree.Branch("selected",&after);
 long long tail=0,positive=0,low=0,sparse=0,span=0,tagged=0,lowSignal=0,combined=0,afterTrack=0,undefined=0,selected=0,base=0;
 for(size_t i=0;i<all.size();++i){id=all[i];t=tags[i];before=!id.tail_rich&&id.positive&&(cog==200||id.pass20);after=before&&(!cfg.enabled||!t.tag);base+=before;selected+=after;tree.Fill();
  if(id.tail_rich)continue;
  ++tail;if(!id.positive)continue;
  ++positive;low+=t.low;sparse+=t.low&&t.sparse;span+=t.low&&t.sparse&&t.spanning;tagged+=t.oldTrack;lowSignal+=t.lowSignal;combined+=t.tag;
  if(cfg.enabled&&t.tag)continue;
  ++afterTrack;undefined+=!id.cog_valid;
 }
 tree.Write();TNamed("definition",cfg.definition().c_str()).Write();TNamed("event_source",source.c_str()).Write();f.Close();
 std::ofstream out(dest/"track_cutflow.tsv");out<<"stage\tn\nall_sampled\t"<<all.size()<<"\nlow_LG_tail_rejected\t"<<all.size()-tail<<"\nafter_low_LG_tail\t"<<tail<<"\nnonpositive_rejected\t"<<tail-positive<<"\nafter_positive\t"<<positive<<"\ncandidate_low_energy\t"<<low<<"\ncandidate_also_sparse\t"<<sparse<<"\ncandidate_also_spanning\t"<<span<<"\ncandidate_also_straight\t"<<tagged<<"\ncandidate_low_signal\t"<<lowSignal<<"\ncandidate_combined_veto\t"<<combined<<"\ntrack_rejected_before_cog\t"<<positive-afterTrack<<"\nafter_track_before_cog\t"<<afterTrack<<"\nundefined_cog_after_track\t"<<undefined<<"\ncog_rejected_after_track\t"<<afterTrack-selected<<"\nselected\t"<<selected<<"\nbaseline_same_cog\t"<<base<<"\ntrack_rejected_same_cog\t"<<base-selected<<'\n';
}
std::vector<double> energies(const Events&v){std::vector<double>x;for(auto&e:v)x.push_back(e.energy);return x;}
void save(TCanvas&c,const fs::path&p){readableAxes(c);TIter pads(c.GetListOfPrimitives());while(auto*o=pads())if(auto*p=dynamic_cast<TPad*>(o))p->Modified();c.Modified();c.Update();c.Write(c.GetName(),TObject::kOverwrite);c.SaveAs(p.c_str());}
void label(const std::string&s){TLatex t;t.SetNDC();t.SetTextSize(.027);t.DrawLatex(.13,.94,s.c_str());}
void overlay(TH1&d,TH1&m,bool normalize=false){d.SetLineColor(kBlue+1);m.SetLineColor(kRed+1);d.SetLineWidth(2);m.SetLineWidth(2);double nd=normalize&&d.Integral(0,d.GetNbinsX()+1)>0?d.Integral(0,d.GetNbinsX()+1):1;double nm=normalize&&m.Integral(0,m.GetNbinsX()+1)>0?m.Integral(0,m.GetNbinsX()+1):1;d.SetMaximum(1.25*std::max(d.GetMaximum()/nd,m.GetMaximum()/nm)*nd);if(normalize){d.DrawNormalized("hist",1);m.DrawNormalized("hist same",1);}else{d.Draw("hist E");m.Draw("hist E same");}auto*l=new TLegend(.65,.73,.89,.88);l->AddEntry(&d,"data","l");l->AddEntry(&m,"MC","l");l->Draw();}
// The strip map and its legend share the same color palette.
const std::array<int,5> stripMapColors={kRed-7,kBlue-7,kMagenta-4,kGray+1,kWhite};
void stripMapLegend(double width) {
  gPad->Range(0,0,1,1);
  TLatex text; text.SetNDC(); text.SetTextSize(.05);
  text.DrawLatex(.08,.88,"Physical strip boundary map");
  text.DrawLatex(.08,.77,Form("Outer band: %g mm",width));
  const char* labels[]={"x edge", "y edge", "x and y edges", "Masked strips", "Interior strips"};
  for (int i=0;i<5;++i) {
    const double y=.65-.10*i;
    auto* box=new TBox(.08,y-.025,.16,y+.025);
    box->SetFillColor(stripMapColors[i]); box->SetFillStyle(1001);
    box->SetLineColor(kGray+2); box->SetBit(TObject::kCanDelete); box->Draw("l");
    text.SetTextAlign(12); text.DrawLatex(.20,y,labels[i]);
  }
  text.SetTextAlign(11); text.SetTextSize(.045);
  text.DrawLatex(.08,.12,"Whole-strip signal; no splitting");
  text.DrawLatex(.08,.04,"NOT escaped energy");
}

void stripMaps(const fs::path&dest,double width,const std::set<int>&mask){auto geom=geometry(width,mask);TFile f((dest/"strip_maps.root").c_str(),"RECREATE");for(int page=0;page<2;++page){TCanvas c(Form("strip_map_%d",page),"",1600,1600);c.Divide(4,4);std::vector<std::unique_ptr<TH2D>> frames;for(int j=0;j<15;++j){int l=page*15+j;c.cd(j+1);gPad->SetLeftMargin(.22);gPad->SetRightMargin(.06);gPad->SetBottomMargin(.20);gPad->SetTopMargin(.12);auto h=std::make_unique<TH2D>(Form("layer%d",l),Form("Layer %d;x [mm];y [mm]",l),1,-120,120,1,-120,120);h->SetDirectory(nullptr);h->GetYaxis()->SetTitleOffset(1.4);h->Draw();for(auto&[id,g]:geom){if(id/100000!=l)continue;auto*b=new TBox(g.x-g.hx,g.y-g.hy,g.x+g.hx,g.y+g.hy);b->SetFillColor(stripMapColors[g.masked?3:g.edgeX&&g.edgeY?2:g.edgeX?0:g.edgeY?1:4]);b->SetLineColor(kGray);b->Draw();}frames.push_back(std::move(h));}c.cd(16);stripMapLegend(width);c.Write();save(c,dest/("strip_map_"+std::to_string(page)+".png"));}}
#include "DisplayPlots.hh"
int main(int argc,char**argv){try{
 gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPadLeftMargin(.18);gStyle->SetPadRightMargin(.06);gStyle->SetPadBottomMargin(.19);gStyle->SetPadTopMargin(.12);TH1::AddDirectory(false);
 auto o=options(argc,argv,{{"input-base",""},{"low-signal-max-energy-gev","0.1"},{"scratch-dir",".scratch/shower_leak_style"},{"style-only","0"},{"base",""},{"energy","5"},{"mask-source",""},{"intervals-only","0"},{"maps-only","0"},{"track-cut","off"},{"track-max-energy-gev","0.03"},{"track-max-hits","60"},{"track-min-layers","15"},{"track-min-end-layers","2"},{"track-max-hits-per-layer","2"},{"track-max-rms-mm","5.3"}});fs::path base=o["base"],inputBase=o["input-base"].empty()?base:fs::path(o["input-base"]);double en=beamEnergy(o["energy"]);Events all[2];double edgeWidth=-1;
 const auto track=TrackConfig::parse(o);const std::string condition=track.v2?"exclude_track_v2":track.enabled?"exclude_track":"exclude";
 auto resultPath=[&](const std::string& sample,int cog){return base/sample/condition/(std::to_string(cog)+"mm")/(energyLabel(en)+"GeV");};
 if(o["maps-only"]=="1" || o["style-only"]=="1") {
  if(o["style-only"]=="1")for(int cog:{200,20})redrawSaved(base,condition,en,cog,fs::path(o["scratch-dir"])/base.filename()/condition/energyLabel(en)/std::to_string(cog));
  for(int s=0;s<2;++s) {
   auto f=openFile((samplePath(inputBase,s,200,en)/"events.root").string());
   double w=object<TParameter<double>>(*f,"edge_mm")->GetVal();
   require(edgeWidth<0||w==edgeWidth,"Different edge definitions"); edgeWidth=w;
  }
   auto dest=resultPath("comparison",200);
  require(fs::is_directory(dest),"Missing saved map directory");
  stripMaps(dest,edgeWidth,readMask(o["mask-source"]));
  return 0;
 }

 if(o["intervals-only"]=="1") {
  for(int cog:{200,20}) {
   fs::path dest=resultPath("comparison",cog);
   TFile saved((dest/"comparison.root").c_str(),"UPDATE");
   require(!saved.IsZombie(),"Cannot open saved comparison canvas");
   if(track.enabled)require(std::string(object<TNamed>(saved,"track_definition")->GetTitle())==track.definition(),"Saved track thresholds differ; pass original thresholds for display-only redraw, or use --plots-only to reselect events");
   for(const auto& name:indicators) {
    auto* previous=object<TCanvas>(saved,name+"_interval_metrics");
    auto values=readIntervalCanvas(*previous);
    std::string caption=name+", "+energyLabel(en)+" GeV, "+(cog==200?std::string("nocut"):std::string("CoG #pm20 mm"))+track.caption();
    delete previous;
    TCanvas updated((name+"_interval_metrics").c_str(),"",1600,1000);
    drawIntervalCanvas(updated,values,caption);
    saved.cd();updated.Write(updated.GetName(),TObject::kOverwrite);
    save(updated,dest/(name+"_interval_metrics.png"));
   }
  }
  return 0;
 }
 for(int s=0;s<2;++s){auto f=openFile((samplePath(inputBase,s,200,en)/"events.root").string());double w=object<TParameter<double>>(*f,"edge_mm")->GetVal();require(edgeWidth<0||w==edgeWidth,"Different edge definitions");edgeWidth=w;auto*t=object<TTree>(*f,"events");Event e;e.branches(*t,false);for(Long64_t j=0;j<t->GetEntries();++j){require(t->GetEntry(j)>0,"Read events failed");all[s].push_back(e);}t->ResetBranchAddresses();}
 const auto layerZ=trackLayerZ();std::vector<TrackTag> tags[2];for(int s=0;s<2;++s)for(const auto&e:all[s])tags[s].push_back(trackTag(e,track,layerZ));
 for(int cog:{200,20}){
  fs::path dest=resultPath("comparison",cog);fs::create_directories(dest);
  std::string caption=energyLabel(en)+" GeV, "+(cog==200?std::string("nocut"):std::string("CoG #pm20 mm"))+", low-LG tail excluded; leakage PROXIES";
  caption+=track.caption();
  Events v[2],beforeTrack[2],removed[2];for(int s=0;s<2;++s)for(size_t i=0;i<all[s].size();++i){auto&e=all[s][i];if(!e.tail_rich&&e.positive&&(cog==200||e.pass20)){beforeTrack[s].push_back(e);if(!track.enabled||!tags[s][i].tag)v[s].push_back(e);else removed[s].push_back(e);}}
  require(!v[0].empty()&&!v[1].empty(),"No events in a comparison sample");Metrics baseline[2]={metrics(energies(v[0])),metrics(energies(v[1]))};
  double emax=0;for(auto&vv:v)for(auto&e:vv)emax=std::max(emax,e.energy);emax=std::max(1e-3,emax*1.01);
  std::unique_ptr<TFile> output[2];for(int s=0;s<2;++s){auto dir=resultPath(s?"mc":"data",cog);fs::create_directories(dir);trackAudit(dir,all[s],tags[s],track,cog,samplePath(inputBase,s,200,en)/"events.root");output[s]=std::make_unique<TFile>((dir/"histograms.root").c_str(),"RECREATE");TNamed("event_source",(samplePath(inputBase,s,200,en)/"events.root").c_str()).Write();TNamed("track_definition",track.definition().c_str()).Write();TParameter<Long64_t>("selected_events",v[s].size()).Write();TNamed("interpretation","Signals in instrumented boundary strips/layers are leakage proxies, NOT escaped energy fractions").Write();}
  TFile canvases((dest/"comparison.root").c_str(),"RECREATE");
  TNamed("track_definition",track.definition().c_str()).Write();
  std::ofstream trackConfig(dest/"track_configuration.txt");trackConfig<<track.definition()<<'\n';
  std::ofstream trackImpact(dest/"track_cut_impact.tsv");trackImpact<<"sample\tstage\tretention\t";metricHeader(trackImpact);trackImpact<<'\n';
  for(int s=0;s<2;++s)for(int stage=0;stage<3;++stage){const auto& ev=stage==0?beforeTrack[s]:stage==1?v[s]:removed[s];trackImpact<<(s?"mc":"data")<<'\t'<<(stage==0?"before":stage==1?"after":"removed")<<'\t'<<double(ev.size())/beforeTrack[s].size()<<'\t';metricRow(trackImpact,metrics(energies(ev)));trackImpact<<'\n';}
  std::ofstream summary(dest/"summary.tsv");summary<<"sample\t";metricHeader(summary);summary<<'\n';for(int s=0;s<2;++s){summary<<(s?"mc":"data")<<'\t';metricRow(summary,baseline[s]);summary<<'\n';}
  std::ofstream flows(dest/"histogram_flows.tsv");flows<<"sample\thistogram\tentries\tzero\tunderflow\toverflow\n";
  auto writeHist=[&](TH1&h,int s,long long zero){output[s]->cd();h.Write();flows<<(s?"mc":"data")<<'\t'<<h.GetName()<<'\t'<<h.GetEntries()<<'\t'<<zero<<'\t'<<h.GetBinContent(0)<<'\t'<<h.GetBinContent(h.GetNbinsX()+1)<<'\n';};
  TH1D he0("energy",";Total calibrated signal [GeV];Probability / bin",120,0,emax),he1(he0);for(auto&e:v[0])he0.Fill(e.energy);for(auto&e:v[1])he1.Fill(e.energy);writeHist(he0,0,0);writeHist(he1,1,0);
  TCanvas overview("overview","",1600,1100);overview.Divide(2,2);overview.cd(1);overlay(he0,he1,true);label(caption);overview.cd(2);gPad->SetLogy();overlay(he0,he1,true);
  TProfile le[2]={{"longitudinal_energy",";Layer (upstream 0 #rightarrow 29 downstream);Mean layer signal [GeV]",30,-.5,29.5},{"longitudinal_energy_mc","",30,-.5,29.5}};
  TProfile lf[2]={{"longitudinal_fraction",";Layer;Mean layer / total signal (proxy)",30,-.5,29.5},{"longitudinal_fraction_mc","",30,-.5,29.5}};
  TProfile lh[2]={{"longitudinal_hits",";Layer;Mean accepted hits",30,-.5,29.5},{"longitudinal_hits_mc","",30,-.5,29.5}};
  std::ofstream longitudinal(dest/"longitudinal.tsv");longitudinal<<"sample\tlayer\tmean_energy\tsem_energy\tmean_fraction\tsem_fraction\tmean_hits\n";
  for(int s=0;s<2;++s){for(auto&e:v[s])for(int l=0;l<30;++l){le[s].Fill(l,e.layer_energy[l]);lf[s].Fill(l,e.layer_fraction[l]);lh[s].Fill(l,e.layer_hits[l]);}output[s]->cd();le[s].Write("longitudinal_energy");lf[s].Write("longitudinal_fraction");lh[s].Write("longitudinal_hits");for(int l=0;l<30;++l)longitudinal<<(s?"mc":"data")<<'\t'<<l<<'\t'<<le[s].GetBinContent(l+1)<<'\t'<<le[s].GetBinError(l+1)<<'\t'<<lf[s].GetBinContent(l+1)<<'\t'<<lf[s].GetBinError(l+1)<<'\t'<<lh[s].GetBinContent(l+1)<<'\n';}
  overview.cd(3);overlay(le[0],le[1]);overview.cd(4);overlay(lf[0],lf[1]);canvases.cd();overview.Write();save(overview,dest/"energy_longitudinal.png");
  std::ofstream bins(dest/"indicator_intervals.tsv"),impact(dest/"impact.tsv"),delta(dest/"relative_variance.tsv"),dependence(dest/"run_dependence.tsv"),tails(dest/"tail_association.tsv");
  bins<<"indicator\tsample\tinterval\tlo\thi\tlow_energy_fraction\t";metricHeader(bins);bins<<'\n';
  impact<<"indicator\tthreshold\tsample\tn_before\tn_above\tabove_fraction\tabove_binomial_se\tabove_wilson_lo\tabove_wilson_hi\tretention\tdelta_mean\tdelta_gaussian\tdelta_central68\tdelta_rms\tlow_energy_boundary\tlow_fraction_before\tlow_fraction_after\t";metricHeader(impact);impact<<'\n';
  delta<<"indicator\tthreshold\twidth_estimator\tdata_retention\tmc_retention\tdifference_before\tdifference_after\tdifference_change\n";
  dependence<<"indicator\tsample\trun\tn\tmean_fraction\tsem_fraction\tmean_absolute_GeV\tsem_absolute\tmean_energy\tstatistics\n";
  tails<<"indicator\tthreshold\tsample\tenergy_region\tn\tn_tagged\ttagged_fraction\tbinomial_se\n";
  const std::vector<double> cuts={0,.001,.003,.01,.03,.05,.1,.2};
  const std::vector<double> boundaries={0,.001,.003,.01,.03,.05,.1,.2,1.0000001};
  for(int k=0;k<5;++k){std::string name=indicators[k];
   TH1D hf0((name+"_fraction").c_str(),(";"+name+" / total signal (leakage proxy);Probability / bin").c_str(),100,0,1.0000001),hf1(hf0);long long zeros[2]{};
   for(int s=0;s<2;++s){auto&h=s?hf1:hf0;for(auto&e:v[s]){double f=fraction(e,k);h.Fill(f);if(f==0)++zeros[s];}writeHist(h,s,zeros[s]);}
   TCanvas dist((name+"_distribution").c_str(),"",1300,650);dist.Divide(2,1);dist.cd(1);overlay(hf0,hf1,true);label(caption);dist.cd(2);gPad->SetLogy();overlay(hf0,hf1,true);label("Exact zero: data="+std::to_string(zeros[0])+", MC="+std::to_string(zeros[1]));canvases.cd();dist.Write();save(dist,dest/(name+"_fraction.png"));
   double amax=0;for(auto&vv:v)for(auto&e:vv)amax=std::max(amax,absolute(e,k));amax=std::max(1e-6,amax*1.01);
   TCanvas corr((name+"_correlations").c_str(),"",1600,1000);corr.Divide(2,2);std::vector<std::unique_ptr<TH2D>> hs;std::vector<std::unique_ptr<TProfile>> ps;
   for(int kind=0;kind<2;++kind){double xmax=kind?amax:1.0000001;double zmax=0;
    for(int s=0;s<2;++s){auto h=std::make_unique<TH2D>(Form("%s_%s_vs_energy_%d",name.c_str(),kind?"absolute":"fraction",s),Form("%s;%s %s;Total calibrated signal [GeV]",s?"MC":"data",name.c_str(),kind?"signal [GeV]":"/ total (proxy)"),60,0,xmax,80,0,emax);auto p=std::make_unique<TProfile>(Form("%s_energy_profile_%d_%d",name.c_str(),kind,s),"",60,0,xmax);for(auto&e:v[s]){double x=kind?absolute(e,k):fraction(e,k);h->Fill(x,e.energy);p->Fill(x,e.energy);}output[s]->cd();h->Write();p->Write();h->Scale(1./v[s].size());h->GetZaxis()->SetTitle("Probability / bin");zmax=std::max(zmax,h->GetMaximum());hs.push_back(std::move(h));ps.push_back(std::move(p));}
    for(int s=0;s<2;++s){int i=2*kind+s;corr.cd(i+1);gPad->SetRightMargin(.19);gPad->SetLeftMargin(.13);hs[i]->SetTitle("");hs[i]->GetZaxis()->SetLabelSize(.028);hs[i]->GetZaxis()->SetTitleSize(.03);hs[i]->GetZaxis()->SetTitleOffset(1.3);hs[i]->SetMaximum(zmax);hs[i]->Draw("colz");ps[i]->SetLineColor(kBlack);ps[i]->SetMarkerStyle(20);ps[i]->SetMarkerSize(.35);ps[i]->Draw("same");label(std::string(s?"MC | ":"data | ")+caption);}
   }canvases.cd();corr.Write();save(corr,dest/(name+"_correlations.png"));
   // Explicit singleton zero interval, then (lo,hi], common for data and MC.
   IntervalValues intervalValues;
   for(int b=0;b<int(boundaries.size());++b)for(int s=0;s<2;++s){std::vector<double> ev;int low=0;double lo=b?boundaries[b-1]:0,hi=boundaries[b];for(auto&e:v[s]){double f=fraction(e,k);if((b==0&&f==0)||(b>0&&f>lo&&f<=hi)){ev.push_back(e.energy);if(e.energy<.8*baseline[s].median)++low;}}auto m=metrics(ev);bins<<name<<'\t'<<(s?"mc":"data")<<'\t'<<b<<'\t'<<lo<<'\t'<<hi<<'\t'<<(ev.empty()?NaN:double(low)/ev.size())<<'\t';metricRow(bins,m);bins<<'\n';double values[]={m.mean,m.gaus,m.w68,m.rms};for(int q=0;q<4;++q)if(std::isfinite(values[q])){intervalValues[q].push_back({b,s,m.n<100,values[q],q==0?m.sem:q==1?m.gaus_err:0});}}
   TCanvas widths((name+"_interval_metrics").c_str(),"",1600,1000);
   drawIntervalCanvas(widths,intervalValues,name+", "+energyLabel(en)+" GeV, "+(cog==200?std::string("nocut"):std::string("CoG #pm20 mm"))+track.caption());
   canvases.cd();widths.Write();save(widths,dest/(name+"_interval_metrics.png"));
   for(double threshold:cuts){Metrics after[2];double retained[2];for(int s=0;s<2;++s){std::vector<double> kept;long long lowBefore=0,lowAfter=0;long long nr[3]{},nt[3]{};for(auto&e:v[s]){bool tag=fraction(e,k)>threshold;bool low=e.energy<.8*baseline[s].median;int region=low?0:e.energy>1.2*baseline[s].median?2:1;++nr[region];nt[region]+=tag;if(low)++lowBefore;if(!tag){kept.push_back(e.energy);if(low)++lowAfter;}}after[s]=metrics(kept);auto&m=after[s];double n=v[s].size(),above=n-kept.size(),p=above/n;retained[s]=kept.size()/n;double den=1+1/n,center=(p+.5/n)/den,half=std::sqrt(p*(1-p)/n+.25/(n*n))/den;
    impact<<name<<'\t'<<threshold<<'\t'<<(s?"mc":"data")<<'\t'<<n<<'\t'<<above<<'\t'<<p<<'\t'<<std::sqrt(p*(1-p)/n)<<'\t'<<center-half<<'\t'<<center+half<<'\t'<<retained[s]<<'\t'<<m.mean-baseline[s].mean<<'\t'<<m.gaus-baseline[s].gaus<<'\t'<<m.w68-baseline[s].w68<<'\t'<<m.rms-baseline[s].rms<<'\t'<<.8*baseline[s].median<<'\t'<<lowBefore/n<<'\t'<<(kept.empty()?NaN:double(lowAfter)/kept.size())<<'\t';metricRow(impact,m);impact<<'\n';for(int r=0;r<3;++r){double p=nr[r]?double(nt[r])/nr[r]:NaN;tails<<name<<'\t'<<threshold<<'\t'<<(s?"mc":"data")<<'\t'<<(r==0?"E<0.8baseline_median":r==1?"central":"E>1.2baseline_median")<<'\t'<<nr[r]<<'\t'<<nt[r]<<'\t'<<p<<'\t'<<(nr[r]?std::sqrt(p*(1-p)/nr[r]):NaN)<<'\n';}}
    for(int q=0;q<3;++q){auto val=[&](const Metrics&m){return q==0?m.gaus:q==1?m.w68:m.rms;};double before=std::pow(val(baseline[0]),2)-std::pow(val(baseline[1]),2),aft=std::pow(val(after[0]),2)-std::pow(val(after[1]),2);delta<<name<<'\t'<<threshold<<'\t'<<(q==0?"gaussian":q==1?"central68":"rms")<<'\t'<<retained[0]<<'\t'<<retained[1]<<'\t'<<before<<'\t'<<aft<<'\t'<<aft-before<<'\n';}
   }
   TCanvas dep((name+"_dependence").c_str(),"",1600,1100);dep.Divide(2,3);std::vector<std::unique_ptr<TProfile>> profiles;
   for(int kind=0;kind<2;++kind)for(int axis=0;axis<2;++axis){for(int s=0;s<2;++s){auto p=std::make_unique<TProfile>(Form("%s_cog_%d_%d_%d",name.c_str(),kind,axis,s),Form(";%s [mm];Mean %s %s",axis?"layer 10 CoG y":"layer 9 CoG x",name.c_str(),kind?"[GeV]":"/ total (proxy)"),24,-120,120);for(auto&e:v[s])if(e.cog_valid)p->Fill(axis?e.cog_y[10]:e.cog_x[9],kind?absolute(e,k):fraction(e,k));output[s]->cd();p->Write();profiles.push_back(std::move(p));}int i=(kind*2+axis)*2;dep.cd(kind*2+axis+1);overlay(*profiles[i],*profiles[i+1]);}
   std::vector<std::unique_ptr<TH1D>> runhist;for(int s=0;s<2;++s){std::map<int,std::array<Stats,3>> runs;for(auto&e:v[s]){auto&r=runs[e.run];r[0].add(fraction(e,k));r[1].add(absolute(e,k));r[2].add(e.energy);}auto h=std::make_unique<TH1D>(Form("%s_run_%d",name.c_str(),s),Form("%s;Run (categorical);Mean %s / total (proxy)",s?"MC":"data",name.c_str()),runs.size(),0,runs.size());int b=0;for(auto&[run,r]:runs){++b;h->GetXaxis()->SetBinLabel(b,std::to_string(run).c_str());h->SetBinContent(b,r[0].mean);h->SetBinError(b,r[0].sem());dependence<<name<<'\t'<<(s?"mc":"data")<<'\t'<<run<<'\t'<<r[0].n<<'\t'<<r[0].mean<<'\t'<<r[0].sem()<<'\t'<<r[1].mean<<'\t'<<r[1].sem()<<'\t'<<r[2].mean<<'\t'<<(r[0].n<100?"LOW_STAT":"OK")<<'\n';}output[s]->cd();h->Write();dep.cd(5+s);h->SetLineColor(s?kRed+1:kBlue+1);h->Draw("E");runhist.push_back(std::move(h));}canvases.cd();dep.Write();save(dep,dest/(name+"_dependence.png"));
  }
  // Per-layer center/spread and shower depth are saved as additional diagnostic profiles.
  for(int s=0;s<2;++s){output[s]->cd();for(int q=0;q<4;++q){const char*names[]={"cog_x","cog_y","width_x","width_y"};TProfile p(names[q],";Layer;Strip-center moment [mm]",30,-.5,29.5);for(auto&e:v[s])for(int l=0;l<30;++l){double a=q==0?e.cog_x[l]:q==1?e.cog_y[l]:q==2?e.width_x[l]:e.width_y[l];if(std::isfinite(a))p.Fill(l,a);}p.Write();}TH1D depth("depth",";Energy-weighted depth [mm];Events",60,0,300),maximum("max_layer",";Maximum-response layer;Events",30,-.5,29.5);for(auto&e:v[s]){depth.Fill(e.depth);maximum.Fill(e.max_layer);}depth.Write();maximum.Write();output[s]->Close();}
  canvases.Close();if(cog==200)stripMaps(dest,edgeWidth,readMask(o["mask-source"]));std::cout<<"Wrote "<<dest<<'\n';
 }
 return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
