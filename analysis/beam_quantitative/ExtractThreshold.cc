#include "../ssa_shower_shape/Common.hh"

// Preserve the sampled event/axis selection of the completed beam evaluation.
int main(int argc, char** argv) {
  try {
    require(argc == 5, "ExtractThreshold manifest beam_events output axis_layers");
    const int axisLayers = std::stoi(argv[4]);
    auto cached = openRoot(argv[2]); auto* original = getTree(*cached,"events");
    int savedFile=0,savedRun=0,savedEvent=0,savedTime=0;
    Long64_t savedEntry=0,cursor=0; bool savedValid=false; double savedX=0,savedY=0;
    std::vector<int>* savedLayer=nullptr;
    original->SetBranchStatus("*",0);
    bind(original,"fileIndex",&savedFile); bind(original,"entry",&savedEntry);
    bind(original,"run",&savedRun); bind(original,"event",&savedEvent);
    bind(original,"eventTime",&savedTime); bind(original,"axisValid",&savedValid);
    bind(original,"axisX",&savedX); bind(original,"axisY",&savedY); bind(original,"layer",&savedLayer);
    auto out=openRoot(argv[3],"RECREATE");
    TTree events("events","Physical strips before SSA; fixed original sampled events");
    int fileIndex=0,run=0,event=0,eventTime=0,invalidPre=0,mergedEntries=0;
    Long64_t entry=0; bool axisValid=false;
    std::vector<int> cellID,layer,postCounts; std::vector<double> x,y,energy;
#define BR(n,t) events.Branch(#n,&n,#n "/" t)
    BR(fileIndex,"I"); BR(entry,"L"); BR(run,"I"); BR(event,"I"); BR(eventTime,"I");
    BR(axisValid,"O"); BR(invalidPre,"I"); BR(mergedEntries,"I");
#undef BR
    events.Branch("cellID",&cellID); events.Branch("layer",&layer);
    events.Branch("x",&x); events.Branch("y",&y); events.Branch("energy",&energy); events.Branch("postCounts",&postCounts);
    std::ifstream manifest(argv[1]); require(bool(manifest),"Cannot open manifest"); std::string line;
    while(std::getline(manifest,line)) {
      if(line.empty()) continue;
      std::istringstream ss(line); std::string prePath,postPath,truthPath,indices;
      require(bool(std::getline(ss,prePath,'\t')) && bool(std::getline(ss,postPath,'\t')) &&
              bool(std::getline(ss,truthPath,'\t')) && bool(std::getline(ss,indices)),"Invalid manifest");
      auto file=openRoot(prePath); Reader pre(getTree(*file,"Calib_Hit"),false);
      std::replace(indices.begin(),indices.end(),',',' '); std::istringstream selected(indices);
      while(selected>>entry) {
        require(cursor<original->GetEntries() && original->GetEntry(cursor++)>0,"Cached event count mismatch");
        require(savedFile==fileIndex && savedEntry==entry,"Cached manifest identity mismatch"); pre.read(entry);
        require(std::tie(pre.run,pre.event,pre.time)==std::tie(savedRun,savedEvent,savedTime),"Cached event identity mismatch");
        auto hits=pre.hits(invalidPre); const auto axis=entranceAxis(hits,axisLayers);
        require(axis.valid==savedValid,"Cached entrance selection changed");
        if(axis.valid) require(std::abs(axis.x-savedX)<1e-9 && std::abs(axis.y-savedY)<1e-9,"Cached entrance axis changed");
        std::map<int,Hit> strips;
        for(auto hit:hits) {
          const int id=physical(hit.id); hit.id=id; auto inserted=strips.emplace(id,hit);
          if(!inserted.second) inserted.first->second.e+=hit.e;
        }
        mergedEntries=static_cast<int>(hits.size()-strips.size());
        run=pre.run; event=pre.event; eventTime=pre.time; axisValid=savedValid;
        cellID.clear(); layer.clear(); x.clear(); y.clear(); energy.clear(); postCounts.assign(NL,0);
        for(const auto& kv:strips) {const auto& h=kv.second; cellID.push_back(h.id); layer.push_back(h.l); x.push_back(h.x); y.push_back(h.y); energy.push_back(h.e);}
        for(int l:*savedLayer) {require(l>=0 && l<NL,"Invalid cached SSA layer"); ++postCounts[l];}
        events.Fill();
      }
      ++fileIndex;
    }
    require(cursor==original->GetEntries() && cursor>0,"Manifest does not cover cached events");
    out->cd(); events.Write();
    TNamed("strip_rule","One physical CellID per event; memory-cell duplicates summed with signed energies before any MIP cut").Write();
    TNamed("selection_rule","Original axis-valid events, unchanged across MIP cuts; uncut includes nonpositive signals").Write();
    std::cout<<"Extracted pre-SSA strips for "<<cursor<<" sampled events\n"; return 0;
  } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n'; return 1;}
}
