#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using EventKey = std::tuple<int, int, int>;
using FilePtr = std::unique_ptr<TFile>;

FilePtr openFile(const std::string& path, const char* mode = "READ") {
    FilePtr f(TFile::Open(path.c_str(), mode));
    if (!f || f->IsZombie()) throw std::runtime_error("Cannot open ROOT file: " + path);
    return f;
}
TTree* tree(TFile& f, const char* name) {
    auto* t = dynamic_cast<TTree*>(f.Get(name));
    if (!t) throw std::runtime_error(std::string("Missing tree ") + name + " in " + f.GetName());
    t->SetBranchStatus("*", false);
    return t;
}
template<class T> void bind(TTree* t, const char* name, T* value) {
    if (!t->GetBranch(name)) throw std::runtime_error(std::string("Missing branch: ") + name);
    t->SetBranchStatus(name, true);
    if (t->SetBranchAddress(name, value) < 0)
        throw std::runtime_error(std::string("Incompatible branch type: ") + name);
}
int physicalID(int c) {
    int layer = c / 100000, chip = c % 100000 / 10000, channel = c % 100;
    if (c < 0 || layer >= 32 || chip >= 6 || channel >= 36)
        throw std::runtime_error("Invalid CellID: " + std::to_string(c));
    return layer * 100000 + chip * 10000 + channel;
}
bool connected(int c) { return c % 100000 / 10000 != 5 || c % 100 < 30; }
int channelSlot(int physical) {
    return (physical / 100000 * 6 + physical % 100000 / 10000) * 36 + physical % 100;
}
size_t hitSlot(int full) {
    const int p = physicalID(full);
    return ((p / 100000 * 6 + p % 100000 / 10000) * 100 + full % 10000 / 100) * 36 + p % 100;
}
struct HitIndex {
    uint64_t calibrated = 0, raw = 0, signal = 0;
    size_t index = 0;
};
void cacheActiveBranches(TTree* t) {
    t->SetCacheSize(0);
    t->SetCacheSize(32 * 1024 * 1024);
    TIter next(t->GetListOfBranches());
    while (auto* branch = next()) if (t->GetBranchStatus(branch->GetName()))
        t->AddBranchToCache(branch->GetName(), true);
    t->StopCacheLearningPhase();
}

struct Options {
    std::string manifest, pedestal, output, sample, channels = "all", layers = "0-29";
    long long maxEvents = 0;
    int bins = 128;
    double hgMax = 4200, lgMax = 3200, hgEnergyMax = 40, energyMax = 400;
};
Options options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 == argc) throw std::runtime_error("Missing argument value");
        std::string k = argv[i], v = argv[i + 1];
        if (k == "--manifest") o.manifest = v;
        else if (k == "--pedestal") o.pedestal = v;
        else if (k == "--output") o.output = v;
        else if (k == "--sample") o.sample = v;
        else if (k == "--channels") o.channels = v;
        else if (k == "--layers") o.layers = v;
        else if (k == "--max-events") o.maxEvents = std::stoll(v);
        else if (k == "--bins") o.bins = std::stoi(v);
        else if (k == "--hg-max") o.hgMax = std::stod(v);
        else if (k == "--lg-max") o.lgMax = std::stod(v);
        else if (k == "--hg-energy-max") o.hgEnergyMax = std::stod(v);
        else if (k == "--energy-max") o.energyMax = std::stod(v);
        else throw std::runtime_error("Unknown argument: " + k);
    }
    if (o.manifest.empty() || o.pedestal.empty() || o.output.empty() || o.sample.empty() ||
        o.bins < 16 || o.bins > 1024 || o.maxEvents < 0 ||
        !(o.hgMax > 3000) || !(o.lgMax > 130) || !(o.hgEnergyMax > 0) || !(o.energyMax > 0))
        throw std::runtime_error("Invalid options; use run.sh --help");
    return o;
}
std::set<int> integerList(const std::string& s, bool ranges) {
    std::set<int> result;
    std::istringstream input(s);
    for (std::string item; std::getline(input, item, ',');) {
        auto dash = item.find('-');
        int lo = std::stoi(item), hi = lo;
        if (dash != std::string::npos) {
            if (!ranges) throw std::runtime_error("Channel ranges are not supported");
            hi = std::stoi(item.substr(dash + 1));
        }
        if (lo < 0 || hi < lo || (ranges && hi >= 32)) throw std::runtime_error("Invalid selection: " + s);
        for (int i = lo; i <= hi; ++i) result.insert(i);
    }
    if (result.empty()) throw std::runtime_error("Empty selection");
    return result;
}
struct Pedestal { double hg, lg; };
template<class Real> std::map<int, Pedestal> readPedestalTree(TTree* t) {
    std::vector<int>* ids = nullptr;
    std::vector<Real>* h = nullptr;
    std::vector<Real>* l = nullptr;
    bind(t, "CellID", &ids); bind(t, "PedHighMean", &h); bind(t, "PedLowMean", &l);
    std::map<int, Pedestal> result;
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
        if (t->GetEntry(e) <= 0) throw std::runtime_error("Cannot read pedestal entry");
        if (ids->size() != h->size() || ids->size() != l->size()) throw std::runtime_error("Pedestal vector length mismatch");
        for (size_t j = 0; j < ids->size(); ++j) {
            int id = physicalID(ids->at(j));
            Pedestal value{double(h->at(j)), double(l->at(j))};
            if (!std::isfinite(value.hg) || !std::isfinite(value.lg)) throw std::runtime_error("Nonfinite pedestal");
            auto [it, fresh] = result.emplace(id, value);
            if (!fresh && (std::abs(it->second.hg - value.hg) > 1e-6 || std::abs(it->second.lg - value.lg) > 1e-6))
                throw std::runtime_error("Conflicting pedestal values for channel " + std::to_string(id));
        }
    }
    t->ResetBranchAddresses();
    delete ids; delete h; delete l;
    return result;
}
std::map<int, Pedestal> readPedestal(const std::string& name) {
    auto f = openFile(name); auto* t = tree(*f, "ChnLevel");
    auto* branch = t->GetBranch("PedHighMean");
    if (!branch) throw std::runtime_error("Missing PedHighMean");
    std::string type = branch->GetClassName();
    if (type.find("double") != std::string::npos) return readPedestalTree<double>(t);
    if (type.find("float") != std::string::npos) return readPedestalTree<float>(t);
    throw std::runtime_error("Unsupported pedestal branch type: " + type);
}

struct Axes { std::string name, x, y; double xmin, xmax, ymin, ymax; };
// Sparse bin counts preserve every hit, including flow bins. Only the channel
// being written is expanded into a dense TH2D; no event/point sampling is used.
struct Plot {
    Axes a;
    int n;
    long long entries = 0, outside = 0;
    std::unordered_map<unsigned int, uint64_t> counts;
    std::unique_ptr<TProfile> profile;
    std::array<double, 7> stats{};
    Plot(Axes axes, int bins) : a(std::move(axes)), n(bins),
        profile(std::make_unique<TProfile>((a.name + "_profile").c_str(), "", n, a.xmin, a.xmax)) {
        profile->SetDirectory(nullptr);
        profile->GetXaxis()->SetTitle(a.x.c_str()); profile->GetYaxis()->SetTitle(a.y.c_str());
    }
    int bin(double x, double low, double high) const {
        if (x < low) return 0;
        if (!(x < high)) return n + 1;
        return 1 + std::min(n - 1, int((x - low) / (high - low) * n));
    }
    void fill(double x, double y) {
        if (!std::isfinite(x) || !std::isfinite(y)) return;
        int bx = bin(x, a.xmin, a.xmax), by = bin(y, a.ymin, a.ymax);
        ++counts[bx + (n + 2) * by]; ++entries;
        if (bx == 0 || by == 0 || bx == n + 1 || by == n + 1) ++outside;
        else {
            stats[0] += 1; stats[1] += 1; stats[2] += x; stats[3] += x*x;
            stats[4] += y; stats[5] += y*y; stats[6] += x*y;
        }
        // Profiles use exact y values (including y beyond the display range).
        profile->Fill(x, y);
    }
    void write() {
        TH2D hist(a.name.c_str(), (";" + a.x + ";" + a.y + ";Hits").c_str(), n, a.xmin, a.xmax, n, a.ymin, a.ymax);
        hist.SetDirectory(nullptr);
        for (const auto& [binID, count] : counts) hist.SetBinContent(binID, double(count));
        hist.PutStats(stats.data()); hist.SetEntries(entries); hist.Write();
        profile->Write();
        TParameter<Long64_t>((a.name + "_outside").c_str(), outside).Write();
    }
};
struct Channel {
    int id;
    Pedestal ped;
    std::vector<std::unique_ptr<Plot>> plots;
    Channel(int c, Pedestal p, const Options& o) : id(c), ped(p) {
        std::vector<Axes> definitions = {
            {"hg_lg", "HG - pedestal [ADC]", "LG - pedestal [ADC]", -100, o.hgMax, -50, o.lgMax},
            {"hg_lg_zoom", "HG - pedestal [ADC]", "LG - pedestal [ADC]", -50, 3000, -5, 130},
            {"hg_energy", "HG - pedestal [ADC]", "Hit_HG_Energy [MeV]", -100, o.hgMax, -5, o.hgEnergyMax},
            {"lg_energy", "LG - pedestal [ADC]", "Hit_LG_Energy [MeV]", -50, o.lgMax, -5, o.energyMax},
            {"hg_selected_energy", "HG - pedestal [ADC]", "Hit_Energy [MeV] (all matched hits)", -100, o.hgMax, -5, o.energyMax},
            {"lg_selected_energy", "LG - pedestal [ADC]", "Hit_Energy [MeV] (all matched hits)", -50, o.lgMax, -5, o.energyMax}
        };
        for (auto& d : definitions) plots.push_back(std::make_unique<Plot>(d, o.bins));
    }
    void write(TDirectory* out) {
        auto* d = out->mkdir(("channel_" + std::to_string(id)).c_str()); d->cd();
        TParameter<int>("physical_cellid", id).Write();
        TParameter<double>("pedestal_hg", ped.hg).Write();
        TParameter<double>("pedestal_lg", ped.lg).Write();
        for (auto& plot : plots) plot->write();
    }
};

void run(const Options& o) {
    TH1::AddDirectory(false);
    const auto ped = readPedestal(o.pedestal);
    const auto layers = integerList(o.layers, true);
    const auto cells = o.channels == "all" ? std::set<int>{} : integerList(o.channels, false);
    for (int c : cells) if (physicalID(c) != c || !connected(c) || !layers.count(c / 100000))
        throw std::runtime_error("Specify a connected physical CellID in the selected layers: " + std::to_string(c));
    std::array<bool, 32 * 6 * 36> selectedSlots{};
    for (int l = 0; l < 32; ++l) for (int chip = 0; chip < 6; ++chip) for (int ch = 0; ch < 36; ++ch) {
        int c = l * 100000 + chip * 10000 + ch;
        selectedSlots[channelSlot(c)] = connected(c) && layers.count(l) && (cells.empty() || cells.count(c));
    }
    auto selected = [&](int c) { return selectedSlots[channelSlot(c)]; };
    std::map<int, std::unique_ptr<Channel>> channels;
    std::array<Channel*, 32 * 6 * 36> channelPointers{};
    // Reuse dense indices across events. A generation counter makes clearing
    // unnecessary; all 100 representable memory-cell values remain distinct.
    std::vector<HitIndex> hitIndex(32 * 6 * 100 * 36);
    auto getChannel = [&](int c) -> Channel& {
        auto*& cached = channelPointers[channelSlot(c)];
        if (cached) return *cached;
        if (!ped.count(c)) throw std::runtime_error("Missing pedestal for channel " + std::to_string(c));
        auto& ptr = channels[c];
        if (!ptr) ptr = std::make_unique<Channel>(c, ped.at(c), o);
        cached = ptr.get();
        return *ptr;
    };
    for (int c : cells) getChannel(c);
    std::ifstream input(o.manifest);
    if (!input) throw std::runtime_error("Cannot open manifest: " + o.manifest);
    std::string manifestText, line;
    long long processed = 0, nonfinite = 0;
    for (; std::getline(input, line);) {
        if (line.empty()) continue;
        manifestText += line + "\n";
        auto sep = line.find('\t');
        if (sep == std::string::npos) throw std::runtime_error("Manifest requires decode TAB calibration");
        if (o.maxEvents && processed >= o.maxEvents) continue;
        auto rawFile = openFile(line.substr(0, sep)); auto calFile = openFile(line.substr(sep + 1));
        auto* raw = tree(*rawFile, "Raw_Hit"); auto* cal = tree(*calFile, "Calib_Hit");
        int rr = 0, rt = 0, ri = 0, cr = 0, ct = 0, ci = 0;
        bind(raw, "Run_Num", &rr); bind(raw, "Event_Time", &rt); bind(raw, "TriggerID", &ri);
        bind(cal, "Run_Num", &cr); bind(cal, "Event_Time", &ct); bind(cal, "Event_Num", &ci);
        std::map<EventKey, Long64_t> lookup;
        std::vector<EventKey> calibrationKeys;
        for (Long64_t e = 0; e < cal->GetEntries(); ++e) {
            if (cal->GetEntry(e) <= 0) throw std::runtime_error("Cannot read calibration event key");
            EventKey key{cr, ct, ci};
            calibrationKeys.push_back(key);
        }
        // Calibration writes one entry per decode entry, but drops CycleID.
        // Real data can repeat (run,time,trigger). For these repeated keys only,
        // use preserved entry order after verifying ALL event keys in both files.
        bool ordered = raw->GetEntries() == cal->GetEntries();
        if (ordered) for (Long64_t e = 0; e < raw->GetEntries(); ++e) {
            if (raw->GetEntry(e) <= 0) throw std::runtime_error("Cannot read raw event key");
            if (EventKey{rr, rt, ri} != calibrationKeys.at(e)) { ordered = false; break; }
        }
        if (!ordered) for (size_t e = 0; e < calibrationKeys.size(); ++e) {
            if (!lookup.emplace(calibrationKeys[e], e).second)
                throw std::runtime_error("Ambiguous duplicate calibration event key; files do not preserve entry order: " + std::string(calFile->GetName()));
        }
        std::vector<int>* rawIDs = nullptr; std::vector<int>* tags = nullptr;
        std::vector<double>* hg = nullptr; std::vector<double>* lg = nullptr;
        std::vector<int>* calIDs = nullptr;
        std::vector<double>* eh = nullptr; std::vector<double>* el = nullptr; std::vector<double>* energy = nullptr;
        bind(raw, "CellID", &rawIDs); bind(raw, "HitTag", &tags);
        bind(raw, "HG_Charge", &hg); bind(raw, "LG_Charge", &lg);
        bind(cal, "CellID", &calIDs); bind(cal, "Hit_HG_Energy", &eh);
        bind(cal, "Hit_LG_Energy", &el); bind(cal, "Hit_Energy", &energy);
        std::vector<bool> usedCalibrationEntries(ordered ? 0 : cal->GetEntries(), false);
        cacheActiveBranches(raw); cacheActiveBranches(cal);
        std::cout << o.sample << ": " << rawFile->GetName() << " (" << raw->GetEntries()
                  << " events; " << (ordered ? "entry order verified" : "unique event-key join") << ")" << std::endl;
        for (Long64_t e = 0; e < raw->GetEntries() && (!o.maxEvents || processed < o.maxEvents); ++e) {
            if (raw->GetEntry(e) <= 0) throw std::runtime_error("Cannot read raw event");
            Long64_t calEntry = e;
            if (!ordered) {
                auto found = lookup.find(EventKey{rr, rt, ri});
                if (found == lookup.end()) throw std::runtime_error("No matching calibration event for run/time/trigger " + std::to_string(rr) + "/" + std::to_string(rt) + "/" + std::to_string(ri));
                calEntry = found->second;
                if (usedCalibrationEntries[calEntry]) throw std::runtime_error("Ambiguous duplicate raw event key");
                usedCalibrationEntries[calEntry] = true;
            }
            if (cal->GetEntry(calEntry) <= 0) throw std::runtime_error("Cannot read matched calibration event");
            if (rawIDs->size() != hg->size() || rawIDs->size() != lg->size() || rawIDs->size() != tags->size() ||
                calIDs->size() != eh->size() || calIDs->size() != el->size() || calIDs->size() != energy->size())
                throw std::runtime_error("Hit vector length mismatch");
            const uint64_t generation = static_cast<uint64_t>(processed) + 1;
            for (size_t j = 0; j < calIDs->size(); ++j) {
                auto& slot = hitIndex[hitSlot(calIDs->at(j))];
                if (slot.calibrated == generation) throw std::runtime_error("Duplicate full CellID in calibration event");
                slot.calibrated = generation; slot.index = j;
            }
            for (size_t j = 0; j < rawIDs->size(); ++j) {
                int full = rawIDs->at(j), c = physicalID(full);
                if (!selected(c)) continue;
                auto& slot = hitIndex[hitSlot(full)];
                if (slot.raw == generation) throw std::runtime_error("Duplicate full CellID in raw event");
                slot.raw = generation;
                if (tags->at(j) != 1) continue;
                slot.signal = generation;
                auto& ch = getChannel(c);
                double qh = hg->at(j) - ch.ped.hg, ql = lg->at(j) - ch.ped.lg;
                if (!std::isfinite(qh) || !std::isfinite(ql)) { ++nonfinite; continue; }
                ch.plots[0]->fill(qh, ql); ch.plots[1]->fill(qh, ql);
                if (slot.calibrated != generation) continue;
                size_t k = slot.index;
                ch.plots[2]->fill(qh, eh->at(k)); ch.plots[3]->fill(ql, el->at(k));
                ch.plots[4]->fill(qh, energy->at(k)); ch.plots[5]->fill(ql, energy->at(k));
            }
            for (int full : *calIDs)
                if (selected(physicalID(full)) && hitIndex[hitSlot(full)].signal != generation)
                    throw std::runtime_error("Calibration hit has no corresponding raw HitTag=1 hit: " + std::to_string(full));
            ++processed;
            if (processed % 10000 == 0) std::cout << "  processed " << processed << " events" << std::endl;
        }
        raw->ResetBranchAddresses(); cal->ResetBranchAddresses();
        delete rawIDs; delete tags; delete hg; delete lg;
        delete calIDs; delete eh; delete el; delete energy;
    }
    if (!processed) throw std::runtime_error("No events processed");
    fs::path outPath(o.output); fs::create_directories(outPath.parent_path());
    const std::string temporary = o.output + ".tmp";
    {
        auto f = openFile(temporary, "RECREATE"); f->cd();
        TNamed("sample", o.sample.c_str()).Write();
        TNamed("pedestal_file", o.pedestal.c_str()).Write();
        TNamed("input_pairs", manifestText.c_str()).Write();
        TNamed("selection", ("HitTag == 1; no energy, CoG, position, or bad-channel cut; layers=" + o.layers + "; channels=" + o.channels).c_str()).Write();
        TNamed("matching", "Within each file pair: unique (Run_Num, Event_Time, TriggerID/Event_Num); repeated keys use preserved entry order only if entry counts and the complete event-key sequence agree. Then join full CellID including memory cell and validate calibration hits against raw HitTag=1. Histograms aggregate physical CellID.").Write();
        TParameter<Long64_t>("processed_events", processed).Write();
        TParameter<Long64_t>("nonfinite_adc_pairs", nonfinite).Write();
        for (auto& [id, ch] : channels) { (void)id; ch->write(f.get()); }
        f->Close();
    }
    fs::rename(temporary, o.output);
    std::cout << "Saved " << channels.size() << " channels from " << processed << " events to " << o.output << std::endl;
}
int main(int argc, char** argv) {
    try { run(options(argc, argv)); return 0; }
    catch (const std::exception& e) { std::cerr << "ERROR: " << e.what() << std::endl; return 1; }
}
