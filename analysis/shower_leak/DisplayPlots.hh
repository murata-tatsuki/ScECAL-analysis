#pragma once
// Rebuild presentation from persisted histogram bins/profile sums and saved
// interval graph points. No calibrated inputs, event loop, selection or fit.
inline void redrawSaved(const fs::path&base,const std::string&condition,double en,int cog,
                        const fs::path&scratch) {
 auto dir=[&](const std::string&s){return base/s/condition/(std::to_string(cog)+"mm")/(energyLabel(en)+"GeV");};
 const auto dest=dir("comparison");
 std::array<IntervalValues,5> intervals;std::string definition;
 {auto saved=openFile((dest/"comparison.root").string());for(int k=0;k<5;++k)intervals[k]=readIntervalCanvas(*object<TCanvas>(*saved,indicators[k]+"_interval_metrics"));if(auto*d=dynamic_cast<TNamed*>(saved->Get("track_definition")))definition=d->GetTitle();}
 std::unique_ptr<TFile> input[2]={openFile((dir("data")/"histograms.root").string()),openFile((dir("mc")/"histograms.root").string())};
 std::vector<std::unique_ptr<TH1>> owned;
 auto get=[&](int s,const std::string&name){auto*h=dynamic_cast<TH1*>(object<TH1>(*input[s],name)->Clone());h->SetDirectory(nullptr);owned.emplace_back(h);return h;};
 const std::string selectionCaption=condition=="exclude_track_v2"?"; v2: low-signal + track veto":condition=="exclude_track"?"; track veto":"";
 const std::string caption=energyLabel(en)+" GeV, "+(cog==200?"nocut":"CoG #pm20 mm")+", low-LG tail excluded; leakage PROXIES"+selectionCaption;
 fs::create_directories(scratch);TFile output((scratch/"comparison.root").c_str(),"RECREATE");require(!output.IsZombie(),"Cannot write styled canvas");
 if(!definition.empty())TNamed("track_definition",definition.c_str()).Write();
 auto*h0=get(0,"energy"),*h1=get(1,"energy");const double n[]={h0->GetEntries(),h1->GetEntries()};
 require(n[0]>0&&n[1]>0,"Empty saved histogram");
 {TCanvas c("overview","",1600,1100);c.Divide(2,2);c.cd(1);overlay(*h0,*h1,true);label(caption);c.cd(2);gPad->SetLogy();overlay(*h0,*h1,true);
  c.cd(3);auto*d=get(0,"longitudinal_energy"),*m=get(1,"longitudinal_energy");m->SetName("longitudinal_energy_mc");overlay(*d,*m);
  c.cd(4);d=get(0,"longitudinal_fraction");m=get(1,"longitudinal_fraction");m->SetName("longitudinal_fraction_mc");overlay(*d,*m);output.cd();save(c,scratch/"energy_longitudinal.png");}
 std::map<std::string,std::array<long long,2>> zeros;
 {std::ifstream in(dest/"histogram_flows.tsv");require(bool(in),"Missing histogram flow table");std::string line;std::getline(in,line);while(std::getline(in,line)){std::istringstream row(line);std::string sample,name;long long entries,z,u,o;row>>sample>>name>>entries>>z>>u>>o;require(bool(row),"Invalid histogram flow row");zeros[name][sample=="mc"]=z;}}
 for(int k=0;k<5;++k){const auto&name=indicators[k];
  {auto*d=get(0,name+"_fraction"),*m=get(1,name+"_fraction");TCanvas c((name+"_distribution").c_str(),"",1300,650);c.Divide(2,1);c.cd(1);overlay(*d,*m,true);label(caption);c.cd(2);gPad->SetLogy();overlay(*d,*m,true);label("Exact zero: data="+std::to_string(zeros[name+"_fraction"][0])+", MC="+std::to_string(zeros[name+"_fraction"][1]));output.cd();save(c,scratch/(name+"_fraction.png"));}
  {TCanvas c((name+"_correlations").c_str(),"",1600,1000);c.Divide(2,2);
   for(int kind=0;kind<2;++kind){TH1*h[2];for(int s=0;s<2;++s){h[s]=get(s,name+(kind?"_absolute":"_fraction")+"_vs_energy_"+std::to_string(s));h[s]->Scale(1./n[s]);h[s]->GetZaxis()->SetTitle("Probability / bin");}double zmax=std::max(h[0]->GetMaximum(),h[1]->GetMaximum());
    for(int s=0;s<2;++s){c.cd(2*kind+s+1);gPad->SetRightMargin(.24);h[s]->SetMaximum(zmax);h[s]->SetTitle("");h[s]->Draw("colz");auto*p=get(s,name+"_energy_profile_"+std::to_string(kind)+"_"+std::to_string(s));p->SetLineColor(kBlack);p->SetMarkerStyle(20);p->SetMarkerSize(.35);p->Draw("same");label(std::string(s?"MC | ":"data | ")+caption);}}
   output.cd();save(c,scratch/(name+"_correlations.png"));}
  {TCanvas c((name+"_interval_metrics").c_str(),"",1600,1000);drawIntervalCanvas(c,intervals[k],name+", "+energyLabel(en)+" GeV, "+(cog==200?"nocut":"CoG #pm20 mm")+selectionCaption);output.cd();save(c,scratch/(name+"_interval_metrics.png"));}
  {TCanvas c((name+"_dependence").c_str(),"",1600,1100);c.Divide(2,3);
   for(int kind=0;kind<2;++kind)for(int axis=0;axis<2;++axis){c.cd(kind*2+axis+1);auto*d=get(0,name+"_cog_"+std::to_string(kind)+"_"+std::to_string(axis)+"_0");auto*m=get(1,name+"_cog_"+std::to_string(kind)+"_"+std::to_string(axis)+"_1");overlay(*d,*m);}
   for(int s=0;s<2;++s){c.cd(5+s);auto*h=get(s,name+"_run_"+std::to_string(s));h->SetLineColor(s?kRed+1:kBlue+1);h->Draw("E");}output.cd();save(c,scratch/(name+"_dependence.png"));}
 }
 output.Close();
 // Publish only after every canvas and PNG was generated successfully.
 for(const auto&file:fs::directory_iterator(scratch))fs::copy_file(file.path(),dest/file.path().filename(),fs::copy_options::overwrite_existing);
 fs::remove_all(scratch);
}
