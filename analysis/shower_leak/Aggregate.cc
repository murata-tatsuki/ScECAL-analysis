#include "Common.hh"
struct Input{std::string path;Long64_t n=0,offset=0;};
int main(int argc,char**argv){try{
 auto o=options(argc,argv,{{"manifest",""},{"sample","data"},{"tail-results",""},{"mask-source",""},{"out",""},{"energy","5"},{"edge-mm","10"},{"max-events","6000"},{"max-files","4"},{"seed","20261005"}});
 require(o["sample"]=="data"||o["sample"]=="mc","Invalid sample");double en=beamEnergy(o["energy"]);
 fs::path dest=o["out"];fs::create_directories(dest/"inputs");require(!fs::exists(dest/"events.root"),"events.root exists; choose another label");
 auto mask=readMask(o["mask-source"]);require(mask==readMask((fs::path(o["tail-results"])/"inputs/comparison_source.cc").string()),"Tail and analysis masks differ");
 auto geom=geometry(std::stod(o["edge-mm"]),mask);
 // Reject results older than ANY calibrated file contributing to the tail medians/common set.
 auto stamp=fs::last_write_time(fs::path(o["tail-results"])/"tail_study.root");
 for(auto name:{"data.tsv","mc.tsv"}){std::ifstream in{fs::path(o["tail-results"])/"inputs"/name};require(bool(in),"Missing tail manifest");std::string line;while(std::getline(in,line)){auto t=line.find('\t');require(t!=std::string::npos,"Malformed tail manifest");auto p=fs::canonical(line.substr(t+1));require(fs::last_write_time(p)<=stamp,"Stale tail result: "+p.string());}}
 hl_tail::Selection tail(o["tail-results"]);
 std::set<std::string> expectedPaths;{std::ifstream in{fs::path(o["tail-results"])/"inputs"/(o["sample"]+".tsv")};std::string line;while(std::getline(in,line)){auto tab=line.find('\t');require(tab!=std::string::npos,"Malformed tail manifest");expectedPaths.insert(fs::canonical(line.substr(tab+1)).string());}}
 std::vector<std::string> paths;{std::ifstream in(o["manifest"]);require(bool(in),"Missing manifest");std::string p;std::set<std::string>unique;while(std::getline(in,p)){p=fs::canonical(p).string();require(expectedPaths.count(p),"Input absent from requested sample tail manifest: "+p);require(unique.insert(p).second,"Duplicate input");paths.push_back(p);}}
 require(!paths.empty(),"Empty inputs");long long maxf=std::stoll(o["max-files"]),maxe=std::stoll(o["max-events"]);require(maxf>=0&&maxe>=0,"Negative limits");
 std::vector<Input> inputs;size_t nf=maxf?std::min<size_t>(maxf,paths.size()):paths.size();Long64_t total=0;
 for(size_t j=0;j<nf;++j){size_t k=std::min(paths.size()-1,size_t((j+.5)*paths.size()/nf));auto f=openFile(paths[k]);auto*t=object<TTree>(*f,"Calib_Hit");require(t->GetEntries()>0,"Empty tree");inputs.push_back({paths[k],t->GetEntries(),total});total+=t->GetEntries();}
 Long64_t count=maxe?std::min<Long64_t>(maxe,total):total;
 std::mt19937_64 rng(std::stoull(o["seed"]));std::vector<std::vector<Long64_t>> entries(nf);size_t file=0;
 for(Long64_t j=0;j<count;++j){Long64_t lo=j*total/count,hi=(j+1)*total/count-1;Long64_t idx=std::uniform_int_distribution<Long64_t>(lo,hi)(rng);while(file+1<nf&&idx>=inputs[file+1].offset)++file;entries[file].push_back(idx-inputs[file].offset);}
 std::ofstream inventory(dest/"inputs/files.tsv");inventory<<"file_index\tpath\tbytes\tmtime_ticks\tentries\tsampled\n";
 std::ofstream schema(dest/"inputs/branches.tsv");schema<<"file_index\tbranch\tclass\ttitle\n";
 std::ofstream map(dest/"strip_map.tsv");map<<"physical_id\tlayer\tx_mm\ty_mm\tz_mm\thalf_x\thalf_y\tboundary_x\tboundary_y\tedge_x\tedge_y\tmasked\n";
 for(auto&[id,g]:geom)map<<id<<'\t'<<id/100000<<'\t'<<g.x<<'\t'<<g.y<<'\t'<<g.z<<'\t'<<g.hx<<'\t'<<g.hy<<'\t'<<g.bx<<'\t'<<g.by<<'\t'<<g.edgeX<<'\t'<<g.edgeY<<'\t'<<g.masked<<'\n';
 TFile output((dest/"events.partial.root").c_str(),"RECREATE");require(!output.IsZombie(),"Cannot create output");TTree tree("events","All sampled calibrated events, including rejected events");Event e;e.branches(tree,true);
 long long rich=0,nonpositive=0,valid=0,pass=0,undefined=0,negativeHits=0,maskedHits=0,disconnectedHits=0,memoryHits=0;double dxmax=0,dymax=0,dzmax=0;
 for(size_t fi=0;fi<nf;++fi){auto&in=inputs[fi];inventory<<fi<<'\t'<<in.path<<'\t'<<fs::file_size(in.path)<<'\t'<<fs::last_write_time(in.path).time_since_epoch().count()<<'\t'<<in.n<<'\t'<<entries[fi].size()<<'\n';
  auto flags=tail.mask(in.path);auto f=openFile(in.path);auto*t=object<TTree>(*f,"Calib_Hit");t->SetBranchStatus("*",0);
  int run=0,time=0,event=0;std::vector<int>*ids=nullptr;std::vector<double>*energy=nullptr,*x=nullptr,*y=nullptr,*z=nullptr;
  bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,"Event_Num",&event);bind(t,"CellID",&ids);bind(t,"Hit_Energy",&energy);bind(t,"Hit_X",&x);bind(t,"Hit_Y",&y);bind(t,"Hit_Z",&z);
  for(auto name:{"Run_Num","Event_Time","Event_Num","CellID","Hit_Energy","Hit_X","Hit_Y","Hit_Z"}){auto*b=t->GetBranch(name);schema<<fi<<'\t'<<name<<'\t'<<b->GetClassName()<<'\t'<<b->GetTitle()<<'\n';}
  for(auto entry:entries[fi]){require(t->GetEntry(entry)>0,"Read failure");e=Event{};e.file_index=fi;e.entry=entry;e.run=run;e.time=time;e.event=event;e.tail_rich=flags.at(entry);require(ids->size()==energy->size()&&ids->size()==x->size()&&ids->size()==y->size()&&ids->size()==z->size(),"Hit vector length mismatch");
   double sx[30]{},sy[30]{},sxx[30]{},syy[30]{},sz=0,sl=0;
   for(size_t j=0;j<ids->size();++j){int full=ids->at(j),id=physical(full);if(!connected(full)){++e.disconnected_hits;continue;}if(id!=full)++memoryHits;const auto&g=geom.at(id);
    require(std::isfinite(energy->at(j))&&std::isfinite(x->at(j))&&std::isfinite(y->at(j))&&std::isfinite(z->at(j)),"Nonfinite calibrated hit");
    dxmax=std::max(dxmax,std::abs(x->at(j)-g.x));dymax=std::max(dymax,std::abs(y->at(j)-g.y));dzmax=std::max(dzmax,std::abs(z->at(j)-g.z));
    require(std::abs(x->at(j)-g.x)<.01&&std::abs(y->at(j)-g.y)<.01&&std::abs(z->at(j)-g.z)<.01,"Hit coordinates disagree with EBUdecode");
    if(g.masked){++e.masked_hits;continue;}double w=energy->at(j)/1000.;if(w<0){++e.negative_hits;e.negative_energy+=w;continue;}
    int l=id/100000;++e.hits;++e.layer_hits[l];e.energy+=w;e.layer_energy[l]+=w;sx[l]+=w*g.x;sy[l]+=w*g.y;sxx[l]+=w*g.x*g.x;syy[l]+=w*g.y*g.y;sz+=w*g.z;sl+=w*l;
    if(l>=27)e.back3+=w;
    if(l>=25)e.back5+=w;
    if(g.edgeX)e.edge_x+=w;
    if(g.edgeY)e.edge_y+=w;
    if(g.edgeX||g.edgeY)e.edge+=w;
   }
   e.positive=e.energy>0;
   double layerSum=0;for(double value:e.layer_energy)layerSum+=value;require(std::abs(layerSum-e.energy)<1e-9*std::max(1.,e.energy),"Layer energy conservation failed");
   require(e.back3<=e.back5+1e-10&&e.back5<=e.energy+1e-10&&e.edge<=e.energy+1e-10&&e.edge_x<=e.edge+1e-10&&e.edge_y<=e.edge+1e-10,"Indicator invariant failed");
   for(int k=0;k<5;++k)e.indicator_fraction[k]=fraction(e,k);
   for(int l=0;l<30;++l){double w=e.layer_energy[l];e.layer_fraction[l]=e.positive?w/e.energy:NaN;e.cog_x[l]=w>0?sx[l]/w:NaN;e.cog_y[l]=w>0?sy[l]/w:NaN;e.width_x[l]=w>0?std::sqrt(std::max(0.,sxx[l]/w-e.cog_x[l]*e.cog_x[l])):NaN;e.width_y[l]=w>0?std::sqrt(std::max(0.,syy[l]/w-e.cog_y[l]*e.cog_y[l])):NaN;}
   if(e.positive){e.depth=sz/e.energy;e.depth_layer=sl/e.energy;e.max_layer=std::max_element(e.layer_energy,e.layer_energy+30)-e.layer_energy;}
   e.cog_valid=e.layer_energy[9]>0&&e.layer_energy[10]>0;e.pass20=e.cog_valid&&std::abs(e.cog_x[9])<20&&std::abs(e.cog_y[9])<20&&std::abs(e.cog_x[10])<20&&std::abs(e.cog_y[10])<20;
   if(e.tail_rich)++rich;else if(!e.positive)++nonpositive;else{++valid;if(!e.cog_valid)++undefined;if(e.pass20)++pass;}
   negativeHits+=e.negative_hits;maskedHits+=e.masked_hits;disconnectedHits+=e.disconnected_hits;tree.Fill();
  }
  t->ResetBranchAddresses();delete ids;delete energy;delete x;delete y;delete z;
  std::cout<<o["sample"]<<" file "<<fi+1<<'/'<<nf<<": "<<entries[fi].size()<<" sampled, complete tail identity verified\n";
 }
 output.cd();tree.Write();TNamed("tail_definition",tail.definition().c_str()).Write();TNamed("sample",o["sample"].c_str()).Write();TNamed("selection","tail_rich==false && energy>0; CoG20 additionally requires |x9|,|y9|,|x10|,|y10|<20 mm; no energy-window cut").Write();TParameter<double>("edge_mm",std::stod(o["edge-mm"])).Write();TParameter<double>("beam_energy_GeV",en).Write();output.Close();fs::rename(dest/"events.partial.root",dest/"events.root");
 std::ofstream flow(dest/"cutflow.tsv");flow<<"stage\tn\ninput_files_available\t"<<paths.size()<<"\ninput_files_selected\t"<<nf<<"\nentries_available_selected_files\t"<<total<<"\nsampled\t"<<count<<"\ntail_rejected\t"<<rich<<"\nafter_tail\t"<<count-rich<<"\nnonpositive_rejected_after_tail\t"<<nonpositive<<"\nnocut_selected\t"<<valid<<"\nundefined_cog_kept_in_nocut\t"<<undefined<<"\ncog20_rejected_after_positive\t"<<valid-pass<<"\ncog20_selected\t"<<pass<<'\n';
 std::ofstream check(dest/"validation.tsv");check<<"check\tvalue\nmax_abs_dx_mm\t"<<dxmax<<"\nmax_abs_dy_mm\t"<<dymax<<"\nmax_abs_dz_mm\t"<<dzmax<<"\nnegative_good_hits\t"<<negativeHits<<"\nmasked_hits\t"<<maskedHits<<"\ndisconnected_hits\t"<<disconnectedHits<<"\nhits_with_memory_digits\t"<<memoryHits<<"\ntail_input_mtime_and_full_identity\tPASS\nenergy_sum_invariant\tPASS\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
