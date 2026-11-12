#include "../ssa_shower_shape/Common.hh"

// Reuse the established strip geometry, branch checks and entrance-axis definition.
int main(int argc, char** argv) {
  try {
    require(argc == 7, "Extract manifest output energy axis_layers split|all data|mc");
    const int axisLayers = std::stoi(argv[4]);
    const bool split = std::string(argv[5]) == "split";
    const bool mc = std::string(argv[6]) == "mc";
    require(axisLayers >= 2 && axisLayers <= 30, "Invalid axis layers");
    require(split || std::string(argv[5]) == "all", "Invalid SSA selection");
    std::ifstream manifest(argv[1]);
    require(bool(manifest), "Cannot open manifest");
    auto out = openRoot(argv[2], "RECREATE");
    TTree events("events", "Paired sampled events; invalid axes retained for cutflow");
    int fileIndex = 0, run = 0, event = 0, eventTime = 0, invalidPre = 0, invalidPost = 0;
    Long64_t entry = 0;
    bool axisValid = false, hasTruth = false;
    double axisX = NAN, axisY = NAN, truthX = NAN, truthY = NAN;
    double primaryX = NAN, primaryY = NAN, primaryZ = NAN;
    int missingParents = 0;
    double missingEnergy = 0;
    std::vector<int> layer;
    std::vector<double> x, y, energy;
#define BR(n,t) events.Branch(#n,&n,#n "/" t)
    BR(fileIndex,"I"); BR(entry,"L"); BR(run,"I"); BR(event,"I"); BR(eventTime,"I");
    BR(axisValid,"O"); BR(hasTruth,"O"); BR(axisX,"D"); BR(axisY,"D");
    BR(truthX,"D"); BR(truthY,"D"); BR(primaryX,"D"); BR(primaryY,"D"); BR(primaryZ,"D");
    BR(invalidPre,"I"); BR(invalidPost,"I"); BR(missingParents,"I"); BR(missingEnergy,"D");
#undef BR
    events.Branch("layer", &layer); events.Branch("x", &x); events.Branch("y", &y); events.Branch("energy", &energy);
    std::string line;
    while (std::getline(manifest, line)) {
      if (line.empty()) continue;
      std::istringstream ss(line);
      std::string a, b, c, indices;
      require(bool(std::getline(ss,a,'\t')) && bool(std::getline(ss,b,'\t')) && bool(std::getline(ss,c,'\t')) && bool(std::getline(ss,indices)), "Invalid manifest");
      auto f = openRoot(a), g = openRoot(b);
      Reader pre(getTree(*f,"Calib_Hit"),false), post(getTree(*g,"SSA_Hit"),true);
      require(pre.t->GetEntries() == post.t->GetEntries(), "Calib/SSA entry count mismatch");
      std::unique_ptr<TFile> tf;
      std::unique_ptr<Truth> truth;
      if (c != "-") {
        require(mc, "Data must not use MC truth");
        tf = openRoot(c); truth = std::make_unique<Truth>(getTree(*tf,"MC_Truth"));
        require(truth->t->GetEntries() == pre.t->GetEntries(), "Truth entry count mismatch");
      }
      std::replace(indices.begin(),indices.end(),',',' ');
      std::istringstream selected(indices);
      while (selected >> entry) {
        require(entry >= 0 && entry < pre.t->GetEntries(), "Entry out of range");
        pre.read(entry); post.read(entry);
        require(std::tie(pre.run,pre.time,pre.event,pre.detector) == std::tie(post.run,post.time,post.event,post.detector), "Calib/SSA event identity mismatch");
        run = pre.run; event = pre.event; eventTime = pre.time;
        auto ah = pre.hits(invalidPre), bh = post.hits(invalidPost);
        std::map<int,double> pa, pb, absPre;
        for (const auto& h : ah) { pa[h.id] += h.e; absPre[h.id] += std::abs(h.e); }
        for (const auto& h : bh) pb[h.id] += h.e;
        for (const auto& kv : pb) {
          require(pa.count(kv.first), "SSA contains a parent absent from calib");
          require(std::abs(pa.at(kv.first)-kv.second) <= 1e-7 + 1e-7*absPre.at(kv.first), "SSA parent energy is not conserved");
        }
        missingParents = 0; missingEnergy = 0;
        for (const auto& kv : pa) if (!pb.count(kv.first)) { ++missingParents; missingEnergy += kv.second; }
        const auto axis = entranceAxis(ah,axisLayers);
        axisValid = axis.valid; axisX = axis.valid ? axis.x : NAN; axisY = axis.valid ? axis.y : NAN;
        hasTruth = bool(truth); truthX = truthY = primaryX = primaryY = primaryZ = NAN;
        if (truth) {
          require(mc && pre.run == 0 && pre.time == entry && pre.event == entry, "MC truth entry identity mismatch");
          truth->read(entry);
          require(std::abs(truth->energy-std::stod(argv[3])*1000) < 1e-4, "MC primary energy mismatch");
          primaryX = truth->x; primaryY = truth->y; primaryZ = truth->z;
          truthX = truth->atX(zLayer(0)); truthY = truth->atY(zLayer(0));
        }
        layer.clear(); x.clear(); y.clear(); energy.clear();
        for (const auto& h : positions(bh,split)) {
          layer.push_back(h.l); x.push_back(h.x); y.push_back(h.y); energy.push_back(h.e);
        }
        events.Fill();
      }
      ++fileIndex;
    }
    require(events.GetEntries() > 0, "No selected events");
    out->cd(); events.Write();
    TNamed("position_rule","Exact (layer,x,y) merged per event; signed energy summed before max(0,sum)").Write();
    TNamed("truth_rule","Initial primary straight-line extrapolation, not transported ECAL-entry truth").Write();
    std::cout << "Extracted " << events.GetEntries() << " paired events\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << "ERROR: " << e.what() << '\n'; return 1; }
}
