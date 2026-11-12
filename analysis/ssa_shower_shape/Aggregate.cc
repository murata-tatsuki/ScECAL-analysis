#include "Common.hh"
struct Pair {
  std::string pre, post, truth;
  Long64_t entries = 0;
};
struct Options {
  fs::path manifest, out, beam;
  std::string sample, energy;
  Long64_t limit = 6000;
  unsigned seed = 20261005;
  int axisLayers = 4;
};
Options options(int argc, char** argv) {
  Options o;
  for (int i = 1; i < argc; i += 2) {
    require(i + 1 < argc, "Option requires value");
    std::string k = argv[i], v = argv[i + 1];
    if (k == "--manifest")o.manifest = v;
    else if (k == "--out")o.out = v;
    else if (k == "--sample")o.sample = v;
    else if (k == "--energy")o.energy = v;
    else if (k == "--max-events")o.limit = std::stoll(v);
    else if (k == "--seed")o.seed = std::stoul(v);
    else if (k == "--axis-layers")o.axisLayers = std::stoi(v);
    else if (k == "--beam-map")o.beam = v;
    else throw std::runtime_error("Unknown option " + k);
  }
  require(!o.manifest.empty() && !o.out.empty() && (o.sample == "data" || o.sample == "mc") && !o.energy.empty(), "Missing required argument");
  require(o.limit >= 0 && o.axisLayers >= 2 && o.axisLayers <= 30, "Invalid limit/axis layers");
  return o;
}
struct CellAudit {
  Long64_t seen = 0, missing = 0, split = 0;
  double pre = 0, post = 0, maxDiff = 0;
};
int main(int argc, char** argv) {
  try {
    gROOT->SetBatch();
    TH1::AddDirectory(false);
    auto o = options(argc, argv);
    std::ifstream in(o.manifest);
    require(bool(in), "Cannot read manifest");
    std::vector<Pair>pairs;
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty())continue;
      std::istringstream s(line);
      Pair p;
      require(bool(std::getline(s, p.pre, '\t')) && bool(std::getline(s, p.post, '\t')) && bool(std::getline(s, p.truth, '\t')), "Manifest needs pre,post,truth TSV columns");
      pairs.push_back(p);
    }
    require(!pairs.empty(), "Empty manifest");
    bool hasTruth = pairs[0].truth != "-";
    require(!hasTruth || o.sample == "mc", "Truth is MC only");
    Long64_t total = 0;
    for (auto& p : pairs) {
      require((p.truth != "-") == hasTruth, "Mixed truth availability is forbidden");
      auto f = openRoot(p.pre), g = openRoot(p.post);
      auto* t = getTree(*f, "Calib_Hit"), *u = getTree(*g, "SSA_Hit");
      p.entries = t->GetEntries();
      require(p.entries == u->GetEntries(), "Pre/post entry counts differ: " + p.pre);
      if (hasTruth) {
        auto h = openRoot(p.truth);
        require(getTree(*h, "MC_Truth")->GetEntries() == p.entries, "Truth entry count mismatch");
      }
      total += p.entries;
    }
    require(total > 0, "Empty input trees");
    Long64_t target = o.limit ? std::min(total, o.limit) : total;
    std::mt19937_64 rng(o.seed);
    std::vector<Long64_t> selected;
    selected.reserve(target);
    for (Long64_t i = 0; i < target; ++i) {
      Long64_t lo = i * total / target, hi = (i + 1) * total / target - 1;
      selected.push_back(std::uniform_int_distribution<Long64_t>(lo, hi)(rng));
    }
    auto diagnostic = o.out / o.sample / "diagnostics" / (o.energy + "GeV");
    fs::create_directories(diagnostic);
    auto df = openRoot(diagnostic / "diagnostics.root", "RECREATE");
    std::map<std::string, std::unique_ptr<Profile>> profiles;
    for (auto& s : stages)for (int f = 0; f < (hasTruth ? 3 : 2); ++f)for (auto& w : weights) {
          std::string k = s + "/" + frames[f] + "/" + w;
          profiles[k] = std::make_unique<Profile>(s + frames[f] + w, s == "post_split");
        }
    TH2D preAxis("entrance_axis", "Entrance axis;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D generated("primary_xy", "Primary vertex, one per sampled event;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D projected("primary_at_layer0", "Truth axis at layer 0;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D primarySelected("primary_at_layer0_selected", "Truth axis at layer 0, common cohort;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D original("post_layer0_all", "Split positions, all sampled events;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D smooth("post_layer0_smooth", "Split positions + uniform +/-2.7 mm, all sampled events;X [mm];Y [mm]", 92, -115, 115, 92, -115, 115);
    TH2D axisResidual("axis_residual", "Entrance centroid minus truth at layer 0;dX [mm];dY [mm]", 100, -100, 100, 100, -100, 100);
    TH2D splitByU("split_by_absu", "SSA tag 1 output points;layer;|u-axis| [mm]", 30, -.5, 29.5, 100, 0, 250);
    TH2D unsplitByU("unsplit_by_absu", "SSA tag 0 output points;layer;|u-axis| [mm]", 30, -.5, 29.5, 100, 0, 250);
    TH2D stripByU("strips_by_absu", "Pre-SSA strip positions;layer;|u-axis| [mm]", 30, -.5, 29.5, 100, 0, 250);
    TH2D energyRadiusPre("pre_energy_vs_absu", "Pre-SSA signed strip energy;|u-axis| [mm];energy [MeV]", 100, 0, 250, 220, -10, 100);
    TH2D energyRadiusPost("post_energy_vs_absu", "SSA signed point energy;|u-axis| [mm];energy [MeV]", 100, 0, 250, 220, -10, 100);
    for (auto* h : {
           &preAxis, &generated, &projected, &primarySelected, &original, &smooth, &axisResidual, &splitByU, &unsplitByU, &stripByU, &energyRadiusPre, &energyRadiusPost
         })h->Sumw2();
    int fileIndex = 0, run = 0, event = 0, eventTime = 0;
    Long64_t entry = 0;
    bool axisValid = false;
    double axisX = 0, axisY = 0, truthX = NAN, truthY = NAN;
    int invalidPre = 0, invalidPost = 0;
    int npre[NL] {}, npost[NL] {}, nsplit[NL] {}, nunsplit[NL] {}, missing[NL] {}, negativePre[NL] {}, negativePost[NL] {};
    double energyPre[NL] {}, energyPost[NL] {}, missingEnergy[NL] {}, positivePre[NL] {}, positivePost[NL] {};
    TTree events("events", "Sampled event audit; all comparisons use axis_valid events");
    events.SetDirectory(nullptr);
#define BR(x,typ) events.Branch(#x,&x,#x "/" typ)
    BR(fileIndex, "I");
    BR(entry, "L");
    BR(run, "I");
    BR(event, "I");
    BR(eventTime, "I");
    BR(axisValid, "O");
    BR(axisX, "D");
    BR(axisY, "D");
    BR(truthX, "D");
    BR(truthY, "D");
    BR(invalidPre, "I");
    BR(invalidPost, "I");
#undef BR
#define AR(x,typ) events.Branch(#x,x,#x "[30]/" typ)
    AR(npre, "I");
    AR(npost, "I");
    AR(nsplit, "I");
    AR(nunsplit, "I");
    AR(missing, "I");
    AR(negativePre, "I");
    AR(negativePost, "I");
    AR(energyPre, "D");
    AR(energyPost, "D");
    AR(missingEnergy, "D");
    AR(positivePre, "D");
    AR(positivePost, "D");
#undef AR
    std::map<int, CellAudit>cellAudit;
    Long64_t accepted = 0, cursor = 0, offset = 0, missingEvents = 0;
    std::uniform_real_distribution<double> jitter(-2.7, 2.7);
    std::mt19937_64 jitterRng(o.seed + 1);
    std::ofstream validation(diagnostic / "validation.tsv");
    validation << "file_index\tentries\tsampled\tidentity\tcommon_parent_energy\ttruth_mapping\n";
    for (size_t fi = 0; fi < pairs.size(); ++fi) {
      auto& p = pairs[fi];
      fileIndex = fi;
      auto f = openRoot(p.pre), g = openRoot(p.post);
      Reader pre(getTree(*f, "Calib_Hit"), false), post(getTree(*g, "SSA_Hit"), true);
      std::unique_ptr<TFile>tf;
      std::unique_ptr<Truth>truth;
      if (hasTruth) {
        tf = openRoot(p.truth);
        truth = std::make_unique<Truth>(getTree(*tf, "MC_Truth"));
      }
      Long64_t sampled = 0;
      while (cursor < target && selected[cursor] < offset + p.entries) {
        entry = selected[cursor++] - offset;
        ++sampled;
        pre.read(entry);
        post.read(entry);
        require(std::tie(pre.run, pre.time, pre.event, pre.detector) == std::tie(post.run, post.time, post.event, post.detector), "Pre/post event identity mismatch: " + p.pre + " entry " + std::to_string(entry));
        run = pre.run;
        event = pre.event;
        eventTime = pre.time;
        auto a = pre.hits(invalidPre), b = post.hits(invalidPost);
        auto axis = entranceAxis(a, o.axisLayers);
        axisValid = axis.valid;
        axisX = axis.valid ? axis.x : NAN;
        axisY = axis.valid ? axis.y : NAN;
        if (truth) {
          require(pre.run == 0 && pre.time == entry && pre.event == entry, "Truth mapping requires Run_Num=0 and Event_Time=Event_Num=entry");
          truth->read(entry);
          require(std::abs(truth->energy - std::stod(o.energy) * 1000) < 1e-4, "Truth primary energy is not requested GeV value");
          truthX = truth->atX(zLayer(0));
          truthY = truth->atY(zLayer(0));
          generated.Fill(truth->x, truth->y);
          projected.Fill(truthX, truthY);
        }
        for (auto* v : {
               npre, npost, nsplit, nunsplit, missing, negativePre, negativePost
             })std::fill(v, v + NL, 0);
        for (auto* v : {
               energyPre, energyPost, missingEnergy, positivePre, positivePost
             })std::fill(v, v + NL, 0.);
        std::map<int, double>pa, pb, paAbs;
        std::set<int>splitParents;
        for (auto& h : a) {
          pa[h.id] += h.e;
          paAbs[h.id] += std::abs(h.e);
          energyPre[h.l] += h.e;
          positivePre[h.l] += std::max(0., h.e);
          negativePre[h.l] += h.e < 0;
        }
        for (auto& h : b) {
          pb[h.id] += h.e;
          if (h.tag == 1)splitParents.insert(h.id);
          energyPost[h.l] += h.e;
          positivePost[h.l] += std::max(0., h.e);
          negativePost[h.l] += h.e < 0;
        }
        for (auto& kv : pb) {
          require(pa.count(kv.first), "SSA parent absent from calibration");
          double diff = std::abs(pa.at(kv.first) - kv.second);
          require(diff <= 1e-7 + 1e-7 * paAbs.at(kv.first), "Stored SSA/calib parent energy mismatch: " + p.pre + " entry " + std::to_string(entry) + " CellID " + std::to_string(kv.first));
        }
        bool lost = false;
        for (auto& kv : pa) {
          auto& ca = cellAudit[physical(kv.first)];
          ++ca.seen;
          ca.pre += kv.second;
          if (!pb.count(kv.first)) {
            ++ca.missing;
            ++missing[kv.first / 100000];
            missingEnergy[kv.first / 100000] += kv.second;
            lost = true;
          } else {
            ca.post += pb.at(kv.first);
            ca.maxDiff = std::max(ca.maxDiff, std::abs(kv.second - pb.at(kv.first)));
          }
          ca.split += splitParents.count(kv.first);
        }
        missingEvents += lost;
        std::vector<Hit>present, splitPre;
        for (auto& h : a) {
          if (pb.count(h.id))present.push_back(h);
          if (splitParents.count(h.id))splitPre.push_back(h);
        }
        std::vector<std::vector<Hit>>hs = {positions(a, false), positions(present, false), positions(splitPre, false), positions(b, false), positions(b, true)};
        for (auto& h : hs[0])++npre[h.l];
        for (auto& h : hs[3])++npost[h.l];
        for (auto& h : hs[4]) {
          ++nsplit[h.l];
          if (h.l == 0) {
            original.Fill(h.x, h.y);
            smooth.Fill(h.x + jitter(jitterRng), h.y + jitter(jitterRng));
          }
        }
        // Unsplit counts are unique positions among tag 0, independent of tag 1 overlap.
        std::vector<Hit>zero;
        for (auto& h : b)if (h.tag == 0)zero.push_back(h);
        for (auto& h : positions(zero, false))++nunsplit[h.l];
        if (axis.valid) {
          ++accepted;
          preAxis.Fill(axis.x, axis.y);
          if (truth) {
            primarySelected.Fill(truthX, truthY);
            axisResidual.Fill(axis.x - truthX, axis.y - truthY);
          }
          for (size_t s = 0; s < stages.size(); ++s)for (int fr = 0; fr < (hasTruth ? 3 : 2); ++fr)for (int w = 0; w < 2; ++w)profiles.at(stages[s] + "/" + frames[fr] + "/" + weights[w])->fill(hs[s], fr, w, axis, truth.get());
          for (auto& h : hs[0]) {
            double u = std::abs(h.l % 2 ? h.x - axis.x : h.y - axis.y);
            stripByU.Fill(h.l, u);
            energyRadiusPre.Fill(u, h.e);
          }
          for (auto& h : hs[3])energyRadiusPost.Fill(std::abs(h.l % 2 ? h.x - axis.x : h.y - axis.y), h.e);
          for (auto& h : hs[4])splitByU.Fill(h.l, std::abs(h.l % 2 ? h.x - axis.x : h.y - axis.y));
          for (auto& h : positions(zero, false))unsplitByU.Fill(h.l, std::abs(h.l % 2 ? h.x - axis.x : h.y - axis.y));
        }
        events.Fill();
      }
      offset += p.entries;
      validation << fi << '\t' << p.entries << '\t' << sampled << "\tOK\tOK\t" << (hasTruth ? "entry_id_checked" : "disabled") << '\n';
      std::cout << o.sample << " file " << fi + 1 << '/' << pairs.size() << " sampled " << sampled << " accepted total " << accepted << std::endl;
    }
    require(cursor == target, "Sampling incomplete");
    require(accepted > 0, "No events with a valid entrance axis");
    df->cd();
    events.Write();
    for (auto* h : {
           &preAxis, &generated, &projected, &primarySelected, &original, &smooth, &axisResidual, &splitByU, &unsplitByU, &stripByU, &energyRadiusPre, &energyRadiusPost
         })h->Write();
    number(*df, "available_events", total);
    number(*df, "sampled_events", target);
    number(*df, "accepted_events", accepted);
    number(*df, "axis_invalid_events", target - accepted);
    number(*df, "events_with_missing_parents", missingEvents);
    number(*df, "has_truth", hasTruth);
    TNamed("truth_matching", "same filename and entry count; Run_Num=0, Event_Time=Event_Num=entry required; MC_Truth has no independent event ID").Write();
    TNamed("energy_weight", "max(0, sum of signed energies at an exact position); occupancy includes negative/zero energy").Write();
    if (!o.beam.empty()) {
      auto beam = openRoot(o.beam);
      auto* h = dynamic_cast<TH2*>(beam->Get("h_hitmap_smooth"));
      require(h, "Missing beam template h_hitmap_smooth");
      df->cd();
      h->Write("beam_input_template");
    }
    df->Close();
    std::ofstream cells(diagnostic / "parent_audit.tsv");
    cells << std::setprecision(12) << "CellID\tlayer\tchip\tchannel\tparent_event_occurrences\tmissing_occurrences\tsplit_occurrences\tpre_signed_MeV\tpost_signed_MeV\tmax_common_difference_MeV\n";
    for (auto& kv : cellAudit) {
      int id = kv.first;
      auto& c = kv.second;
      cells << id << '\t' << id / 100000 << '\t' << (id % 100000) / 10000 << '\t' << id % 100 << '\t' << c.seen << '\t' << c.missing << '\t' << c.split << '\t' << c.pre << '\t' << c.post << '\t' << c.maxDiff << '\n';
    }
    std::ofstream cut(diagnostic / "cutflow.tsv");
    cut << "stage\tevents\navailable\t" << total << "\nsampled_and_paired\t" << target << "\naxis_valid_common_cohort\t" << accepted << "\naxis_invalid\t" << target - accepted << "\nwith_missing_parents_diagnostic_only\t" << missingEvents << '\n';
    for (auto& kv : profiles)kv.second->save(o.out / o.sample / kv.first / (o.energy + "GeV"));
    std::cout << "Complete: " << diagnostic << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
}
