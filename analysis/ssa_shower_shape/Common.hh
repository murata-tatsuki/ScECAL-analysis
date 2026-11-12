#pragma once
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include "../calibration_parameters/EBUdecode.h"
namespace fs = std::filesystem;
inline void require(bool b, const std::string& s) {
  if (!b)throw std::runtime_error(s);
}
inline constexpr int NL = 30;
inline const std::vector<std::string> stages = {"pre", "pre_present", "pre_split_parents", "post_all", "post_split"};
inline const std::vector<std::string> frames = {"lab", "aligned", "truth"};
inline const std::vector<std::string> weights = {"count", "positive_energy"};
inline double zLayer(int l) {
  return l % 2 ? 12.2 + (l - 1) / 2 * 19.9 : 1 + l / 2 * 19.9;
}
inline int physical(int id) {
  return id / 100000 * 100000 + (id % 100000) / 10000 * 10000 + id % 100;
}
inline bool connected(int id) {
  if (id < 0)return false;
  int l = id / 100000, c = (id % 100000) / 10000, h = id % 100;
  return l < 30 && c < 6 && h < (c == 5 ? 30 : 36);
}
inline std::unique_ptr<TFile> openRoot(const fs::path& p, const char* mode = "READ") {
  auto f = std::unique_ptr<TFile>(TFile::Open(p.c_str(), mode));
  require(f && !f->IsZombie(), "Cannot open " + p.string());
  return f;
}
inline TTree* getTree(TFile& f, const char* n) {
  auto* t = dynamic_cast<TTree*>(f.Get(n));
  require(t, std::string("Missing tree ") + n + " in " + f.GetName());
  return t;
}
template<class T> inline void bind(TTree* t, const char* n, T* p) {
  require(t->GetBranch(n), std::string("Missing branch ") + n);
  t->SetBranchStatus(n, 1);
  require(t->SetBranchAddress(n, p) >= 0, std::string("Branch type mismatch ") + n);
}
struct Hit {
  int id, l, tag;
  double x, y, z, e;
};
struct Reader {
  TTree* t;
  int run = 0, time = 0, event = 0, detector = 0;
  std::vector<int>* ids = nullptr, *tags = nullptr;
  std::vector<double>* x = nullptr, *y = nullptr, *z = nullptr, *e = nullptr;
  bool ssa;
  Reader(TTree* tree, bool s): t(tree), ssa(s) {
    t->SetBranchStatus("*", 0);
    bind(t, "Run_Num", &run);
    bind(t, "Event_Time", &time);
    bind(t, "Event_Num", &event);
    bind(t, "DetectorID", &detector);
    bind(t, "CellID", &ids);
    bind(t, "Hit_X", &x);
    bind(t, "Hit_Y", &y);
    bind(t, "Hit_Z", &z);
    bind(t, "Hit_Energy", &e);
    if (ssa)bind(t, "ssaTag", &tags);
  }
  void read(Long64_t n) {
    require(t->GetEntry(n) > 0, "Cannot read entry");
    size_t k = ids->size();
    require(x->size() == k && y->size() == k && z->size() == k && e->size() == k && (!ssa || tags->size() == k), "Vector length mismatch");
  }
  std::vector<Hit> hits(int& invalid)const {
    std::vector<Hit> out;
    invalid = 0;
    for (size_t i = 0; i < ids->size(); ++i) {
      int id = ids->at(i);
      require(std::isfinite(x->at(i)) && std::isfinite(y->at(i)) && std::isfinite(z->at(i)) && std::isfinite(e->at(i)), "Nonfinite hit");
      if (!connected(id)) {
        ++invalid;
        continue;
      }
      int l = id / 100000;
      require(std::abs(z->at(i) - zLayer(l)) < 1e-4, "Unexpected layer Z");
      if (!ssa) {
        auto* p = EBUdecode(l, (id % 100000) / 10000, id % 100);
        require(std::abs(x->at(i) - p[0]) < 1e-4 && std::abs(y->at(i) - p[1]) < 1e-4, "Calib strip geometry mismatch");
      }
      if (ssa)require(tags->at(i) == 0 || tags->at(i) == 1, "Invalid ssaTag");
      out.push_back({id, l, ssa ? tags->at(i) : -1, x->at(i), y->at(i), z->at(i), e->at(i)});
    }
    return out;
  }
};
struct Truth {
  TTree* t;
  int pdg = 0;
  double energy = 0, x = 0, y = 0, z = 0;
  std::vector<double>* dir = nullptr;
  Truth(TTree* tree): t(tree) {
    t->SetBranchStatus("*", 0);
    bind(t, "PrimaryPDG", &pdg);
    bind(t, "PrimaryEnergy", &energy);
    bind(t, "PrimaryPosX", &x);
    bind(t, "PrimaryPosY", &y);
    bind(t, "PrimaryPosZ", &z);
    bind(t, "PrimaryDirection", &dir);
  }
  void read(Long64_t n) {
    require(t->GetEntry(n) > 0, "Cannot read truth");
    require(pdg == 11 && dir && dir->size() == 3, "Invalid electron primary");
    require(std::isfinite(x) && std::isfinite(y) && std::isfinite(z) && std::isfinite(energy) && std::all_of(dir->begin(), dir->end(), [](double v) {
      return std::isfinite(v);
    }) && dir->at(2) > 0, "Invalid primary direction/position");
  }
  double atX(double zz)const {
    return x + (zz - z) * dir->at(0) / dir->at(2);
  }
  double atY(double zz)const {
    return y + (zz - z) * dir->at(1) / dir->at(2);
  }
};
struct Axis {
  bool valid = false;
  double x = 0, y = 0;
};
inline Axis entranceAxis(const std::vector<Hit>& h, int layers) {
  double sx = 0, sy = 0, ex = 0, ey = 0;
  for (auto& p : h)if (p.l < layers && p.e > 0) {
      if (p.l % 2) {
        sx += p.x * p.e;
        ex += p.e;
      } else {
        sy += p.y * p.e;
        ey += p.e;
      }
    }
  return {ex > 0 && ey > 0, ex > 0 ? sx / ex : 0, ey > 0 ? sy / ey : 0};
}
// One occupancy count per exact (layer,x,y), with energies summed at duplicate positions.
inline std::vector<Hit> positions(const std::vector<Hit>& hits, bool split) {
  std::map<std::tuple<int, double, double>, Hit> out;
  for (auto h : hits) {
    if (split && h.tag != 1)continue;
    auto key = std::make_tuple(h.l, h.x, h.y);
    auto it = out.find(key);
    if (it == out.end())out[key] = h;
    else it->second.e += h.e;
  }
  std::vector<Hit> v;
  for (auto& kv : out)v.push_back(kv.second);
  return v;
}
struct Moments {
  double w = 0, s = 0, s2 = 0;
  void add(double x, double ww) {
    w += ww;
    s += ww * x;
    s2 += ww * x * x;
  } double mean()const {
    return w > 0 ? s / w : NAN;
  } double rms()const {
    return w > 0 ? std::sqrt(std::max(0., s2 / w - s / w * s / w)) : NAN;
  }
};
inline double containment(const TH1D& h, double q) {
  double total = h.Integral(0, h.GetNbinsX() + 1);
  if (total <= 0)return NAN;
  double sum = h.GetBinContent(0);
  for (int i = 1; i <= h.GetNbinsX(); ++i) {
    double n = h.GetBinContent(i);
    if (sum + n >= q * total && n > 0)return h.GetBinLowEdge(i) + h.GetBinWidth(i) * (q * total - sum) / n;
    sum += n;
  }
  return NAN;
}
inline void number(TDirectory& d, const char* n, Long64_t v) {
  d.cd();
  TParameter<Long64_t>(n, v).Write();
}
struct Profile {
  std::array<std::unique_ptr<TH1D>, NL> u, absu, radial;
  std::array<std::unique_ptr<TH2D>, NL> maps;
  std::array<Moments, NL> moments;
  Long64_t events = 0;
  bool map;
  Profile(const std::string& key, bool m): map(m) {
    for (int l = 0; l < NL; ++l) {
      u[l] = std::make_unique<TH1D>((key + "_u" + std::to_string(l)).c_str(), "", 200, -250, 250);
      absu[l] = std::make_unique<TH1D>((key + "_abs" + std::to_string(l)).c_str(), "", 250, 0, 250);
      for (auto* h : {
             u[l].get(), absu[l].get()
           }) {
        h->SetDirectory(nullptr);
        h->Sumw2();
      }
      if (map) {
        maps[l] = std::make_unique<TH2D>((key + "_xy" + std::to_string(l)).c_str(), "", 84, -226.6, 226.6, 84, -226.6, 226.6);
        maps[l]->SetDirectory(nullptr);
        maps[l]->Sumw2();
        radial[l] = std::make_unique<TH1D>((key + "_r" + std::to_string(l)).c_str(), "", 350, 0, 350);
        radial[l]->SetDirectory(nullptr);
        radial[l]->Sumw2();
      }
    }
  }
  void fill(const std::vector<Hit>& h, int frame, int weight, const Axis& a, const Truth* tr) {
    ++events;
    for (auto& p : h) {
      double x = p.x, y = p.y;
      if (frame == 1) {
        x -= a.x;
        y -= a.y;
      }
      if (frame == 2) {
        x -= tr->atX(p.z);
        y -= tr->atY(p.z);
      }
      double w = weight ? std::max(0., p.e) : 1.;
      if (w == 0)continue;
      double v = p.l % 2 ? x : y;
      u[p.l]->Fill(v, w);
      absu[p.l]->Fill(std::abs(v), w);
      moments[p.l].add(v, w);
      if (map) {
        maps[p.l]->Fill(x, y, w);
        radial[p.l]->Fill(std::hypot(x, y), w);
      }
    }
  }
  void save(const fs::path& p) {
    fs::create_directories(p);
    auto f = openRoot(p / "histograms.root", "RECREATE");
    number(*f, "events", events);
    std::ofstream table(p / "summary.tsv");
    table << std::setprecision(12) << "layer\tfine_axis\tevents\tweight_sum\tmean_mm\trms_mm\tabs_q50_mm\tabs_q90_mm\tabs_tail30_fraction\tu_underflow\tu_overflow\tradial_q50_mm\tradial_q90_mm\n";
    for (int l = 0; l < NL; ++l) {
      u[l]->Write(("u_" + std::to_string(l)).c_str());
      absu[l]->Write(("abs_u_" + std::to_string(l)).c_str());
      if (map) {
        maps[l]->Write(("xy_" + std::to_string(l)).c_str());
        radial[l]->Write(("r_" + std::to_string(l)).c_str());
      }
      auto& m = moments[l];
      table << l << '\t' << (l % 2 ? "x" : "y") << '\t' << events << '\t' << m.w << '\t' << m.mean() << '\t' << m.rms() << '\t' << containment(*absu[l], .5) << '\t' << containment(*absu[l], .9) << '\t' << (m.w ? absu[l]->Integral(31, 251) / m.w : NAN) << '\t' << u[l]->GetBinContent(0) << '\t' << u[l]->GetBinContent(201) << '\t' << (map ? containment(*radial[l], .5) : NAN) << '\t' << (map ? containment(*radial[l], .9) : NAN) << '\n';
    }
  }
};
