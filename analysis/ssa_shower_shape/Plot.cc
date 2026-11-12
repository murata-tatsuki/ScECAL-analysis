#include "Common.hh"
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGraph.h>
struct Series {
  std::unique_ptr<TFile> file;
  Long64_t n = 0;
  std::string name;
  int color;
};
Long64_t param(TFile& f, const char* n) {
  auto* p = dynamic_cast<TParameter<Long64_t>*>(f.Get(n));
  require(p, std::string("Missing parameter ") + n);
  return p->GetVal();
}
std::unique_ptr<TH1D> one(TFile& f, const std::string& name) {
  auto* h = dynamic_cast<TH1D*>(f.Get(name.c_str()));
  require(h, "Missing histogram " + name);
  auto c = std::unique_ptr<TH1D>(static_cast<TH1D*>(h->Clone()));
  c->SetDirectory(nullptr);
  return c;
}
std::unique_ptr<TH2D> two(TFile& f, const std::string& name) {
  auto* h = dynamic_cast<TH2D*>(f.Get(name.c_str()));
  require(h, "Missing map " + name);
  auto c = std::unique_ptr<TH2D>(static_cast<TH2D*>(h->Clone()));
  c->SetDirectory(nullptr);
  return c;
}
void normalize(TH1& h) {
  double sum = h.Integral(0, h.GetNbinsX() + 1);
  if (sum > 0)h.Scale(1 / sum);
}
void style(TH1& h, int color, int line) {
  h.SetLineColor(color);
  h.SetLineWidth(2);
  h.SetLineStyle(line);
  h.SetStats(false);
}
void save(TCanvas& c, TFile& f, const fs::path& p, const std::string& n) {
  c.SaveAs((p / (n + ".png")).c_str());
  f.cd();
  c.Write(n.c_str());
}
// Human-readable labels are presentation only; persisted keys/paths retain their names.
std::string stageLabel(const std::string& stage) {
  const std::map<std::string, std::string> labels = {
    {"pre", "Before SSA: all strips"},
    {"pre_present", "Before SSA: strips retained after SSA"},
    {"pre_split_parents", "Before SSA: strips that were split"},
    {"post_all", "After SSA: all hits (split + unsplit)"},
    {"post_split", "After SSA: successfully split hits"}
  };
  return labels.at(stage);
}
std::string frameLabel(const std::string& frame) {
  if (frame == "lab")return "Detector coordinates";
  if (frame == "aligned")return "Centered on entrance centroid";
  return "Relative to simulated electron track";
}
std::string weightLabel(const std::string& weight) {
  return weight == "count" ? "Hit-position counts" : "Positive signal-energy weights";
}
std::string sampleLabel(const std::string& sample) {
  return sample == "data" ? "Experimental data" : "Simulation";
}
std::string positionAxis(int dimension, const std::string& frame) {
  return std::string(dimension ? "Y" : "X") + (frame == "lab" ? " position [mm]" : " offset [mm]");
}
struct Page {
  TPad* header;
  TPad* body;
};
void headerText(double x, double y, const std::string& text, int size, bool bold = false) {
  TLatex label;
  label.SetNDC();
  label.SetTextFont(bold ? 63 : 43);
  label.SetTextSize(size);
  label.SetTextAlign(12);
  label.DrawLatex(x, y, text.c_str());
}
Page page(TCanvas& canvas, const std::string& title, const std::string& note,
          int nx, int ny, double headerHeight = .14) {
  canvas.cd();
  auto* header = new TPad((std::string(canvas.GetName()) + "_header").c_str(), "", 0, 1 - headerHeight, 1, 1);
  auto* body = new TPad((std::string(canvas.GetName()) + "_body").c_str(), "", 0, 0, 1, 1 - headerHeight);
  header->SetMargin(0, 0, 0, 0);
  header->Draw();
  body->Draw();
  header->cd();
  headerText(.025, .84, title, 30, true);
  headerText(.025, .62, note, 23);
  body->cd();
  body->Divide(nx, ny);
  return {header, body};
}
std::unique_ptr<TLegend> headerLegend(const Page& page, int columns = 2, bool tall = false, double right = .975) {
  page.header->cd();
  auto legend = std::make_unique<TLegend>(.025, .03, right, tall ? .48 : .40);
  legend->SetNColumns(columns);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(43);
  legend->SetTextSize(24);
  legend->SetMargin(.065);
  return legend;
}
void panelStyle(bool map = false) {
  gPad->SetLeftMargin(map ? .13 : .16);
  gPad->SetRightMargin(map ? .20 : .05);
  gPad->SetBottomMargin(.16);
  gPad->SetTopMargin(.13);
}
void profiles(const fs::path& base, const std::string& energy, const std::string& stage, const std::string& frame, const std::string& weight) {
  std::vector<Series>ss;
  for (auto sample : {
         std::string("data"), std::string("mc")
       }) {
    if (frame == "truth" && sample == "data")continue;
    auto p = base / sample / stage / frame / weight / (energy + "GeV") / "histograms.root";
    if (frame == "truth" && !fs::exists(p))return;
    auto f = openRoot(p);
    auto n = param(*f, "events");
    ss.push_back({std::move(f), n, sample, sample == "data" ? kBlack : kRed + 1});
  }
  auto dest = base / "comparison" / stage / frame / weight / (energy + "GeV");
  fs::create_directories(dest);
  auto out = openRoot(dest / "comparison.root", "RECREATE");
  for (bool log : {
         false, true
       }) {
    TCanvas c("fine_profiles", "", 2400, 2250);
    // Objects stay alive until their canvas has been serialized.
    std::vector<std::unique_ptr<TH1D>>keep;
    auto layout = page(c, energy + " GeV | " + stageLabel(stage) + " | " + weightLabel(weight),
                       frameLabel(frame) + "; each layer normalized to unit total weight" + (log ? "; logarithmic vertical scale" : ""), 5, 6, .10);
    auto legend = headerLegend(layout);
    for (int l = 0; l < NL; ++l) {
      layout.body->cd(l + 1);
      panelStyle();
      gPad->SetLogy(log);
      double ymax = 0;
      std::vector<TH1D*>hs;
      for (auto& sample : ss) {
        auto h = one(*sample.file, "u_" + std::to_string(l));
        h->SetName((sample.name + "_u_" + std::to_string(l)).c_str());
        normalize(*h);
        style(*h, sample.color, sample.name == "data" ? 1 : 2);
        ymax = std::max(ymax, h->GetMaximum());
        hs.push_back(h.get());
        keep.push_back(std::move(h));
      }
      for (size_t j = 0; j < hs.size(); ++j) {
        auto* h = hs[j];
        auto title = "Layer " + std::to_string(l) + " - " + (l % 2 ? "X" : "Y") + " direction";
        h->SetTitle((title + ";" + positionAxis(l % 2 ? 0 : 1, frame) + ";" +
                     (weight == "count" ? "Hit fraction / bin" : "Positive signal fraction / bin")).c_str());
        h->SetMaximum(std::max(1e-8, ymax) * (log ? 8 : 1.2));
        h->SetMinimum(log ? 1e-6 : 0);
        h->GetXaxis()->SetRangeUser(-120, 120);
        h->Draw(j ? "hist same" : "hist");
        if (l == 0)legend->AddEntry(h, (sampleLabel(ss[j].name) + " (" + std::to_string(ss[j].n) + " events)").c_str(), "l");
      }
    }
    layout.header->cd();
    legend->Draw();
    save(c, *out, dest, log ? "fine_profiles_log" : "fine_profiles");
  }
  if (stage != "post_split")return;
  // Per-event maps retain the common Z range for data/MC at each layer.
  std::vector<std::unique_ptr<TH2D>>maps;
  double maxima[NL] {};
  for (auto& sample : ss)for (int l = 0; l < NL; ++l) {
      auto h = two(*sample.file, "xy_" + std::to_string(l));
      h->SetName((sample.name + "_xy_" + std::to_string(l)).c_str());
      if (sample.n > 0)h->Scale(1. / sample.n);
      maxima[l] = std::max(maxima[l], h->GetMaximum());
      maps.push_back(std::move(h));
    }
  for (size_t j = 0; j < ss.size(); ++j)for (int parity = 0; parity < 2; ++parity)for (bool log : {
            false, true
          }) {
        TCanvas c("maps", "", 2560, 1600);
        auto layout = page(c, energy + " GeV | " + sampleLabel(ss[j].name) + " | " + stageLabel(stage),
                           frameLabel(frame) + " | " + weightLabel(weight) + " | " + (parity ? "Odd" : "Even") + " layers", 4, 4, .10);
        layout.header->cd();
        headerText(.025, .28, "Per-event signal; common color range for data and simulation at each layer" +
                   std::string(log ? "; logarithmic color scale" : ""), 23);
        for (int l = parity; l < NL; l += 2) {
          layout.body->cd(l / 2 + 1);
          panelStyle(true);
          gPad->SetLogz(log);
          auto* h = maps[j * NL + l].get();
          h->SetTitle(("Layer " + std::to_string(l) + ";" + positionAxis(0, frame) + ";" + positionAxis(1, frame) + ";" +
                       (weight == "count" ? "Positions / event / bin" : "Positive signal [MeV] / event / bin")).c_str());
          h->GetXaxis()->SetRangeUser(-115, 115);
          h->GetYaxis()->SetRangeUser(-115, 115);
          h->GetZaxis()->SetTitleSize(.043);
          h->GetZaxis()->SetTitleOffset(1.25);
          h->SetMaximum(std::max(maxima[l], 1e-9));
          h->SetMinimum(log ? std::max(maxima[l] * 1e-4, 1e-12) : 0);
          h->Draw("colz");
        }
        save(c, *out, dest, ss[j].name + "_map_" + (parity ? "odd" : "even") + (log ? "_log" : ""));
      }
}
std::vector<std::vector<std::string>> readTable(const fs::path& p) {
  std::ifstream f(p);
  require(bool(f), "Missing table " + p.string());
  std::string line;
  std::getline(f, line);
  std::vector<std::vector<std::string>>out;
  while (std::getline(f, line)) {
    std::istringstream s(line);
    std::vector<std::string>r;
    std::string v;
    while (std::getline(s, v, '\t'))r.push_back(v);
    out.push_back(r);
  }
  return out;
}
void widthOverview(const fs::path& base, const std::string& energy, const fs::path& dest, TFile& out) {
  std::ofstream summary(dest / "summary.tsv");
  summary << "sample\tstage\tframe\tweight\tlayer\tfine_axis\tevents\tweight_sum\tmean_mm\trms_mm\tabs_q50_mm\tabs_q90_mm\tabs_tail30_fraction\tu_underflow\tu_overflow\tradial_q50_mm\tradial_q90_mm\n";
  for (auto sample : {
         "data", "mc"
       })for (auto& s : stages)for (auto& fr : frames)for (auto& w : weights) {
          auto p = base / sample / s / fr / w / (energy + "GeV") / "summary.tsv";
          if (!fs::exists(p))continue;
          std::ifstream in(p);
          std::string line;
          std::getline(in, line);
          while (std::getline(in, line))summary << sample << '\t' << s << '\t' << fr << '\t' << w << '\t' << line << '\n';
        }
  for (auto& w : weights) {
    TCanvas c("width_overview", "", 2000, 2400);
    auto layout = page(c, energy + " GeV | Transverse width by layer | " + weightLabel(w),
                       "Left: detector coordinates | Right: centered on entrance centroid", 2, 5, .10);
    std::vector<std::unique_ptr<TGraph>>gs;
    auto legend = headerLegend(layout, 2, false, .53);
    layout.header->cd();
    headerText(.56, .25, "Measured direction: X in odd layers, Y in even layers", 23);
    int panel = 0;
    for (auto& stage : stages)for (int fr = 0; fr < 2; ++fr) {
        layout.body->cd(++panel);
        panelStyle();
        std::vector<TGraph*>pair;
        double max = 1;
        for (auto sample : {
               "data", "mc"
             }) {
          auto rows = readTable(base / sample / stage / frames[fr] / w / (energy + "GeV") / "summary.tsv");
          auto g = std::make_unique<TGraph>();
          for (auto& row : rows) {
            require(row.size() == 13, "Bad summary table");
            double value = std::stod(row[5]);
            if (std::isfinite(value)) {
              g->SetPoint(g->GetN(), std::stod(row[0]), value);
              max = std::max(max, value);
            }
          }
          g->SetName((stage + "_" + frames[fr] + "_" + sample).c_str());
          g->SetLineColor(std::string(sample) == "data" ? kBlack : kRed + 1);
          g->SetMarkerColor(g->GetLineColor());
          g->SetMarkerStyle(std::string(sample) == "data" ? 20 : 24);
          g->SetLineWidth(2);
          pair.push_back(g.get());
          gs.push_back(std::move(g));
        }
        gPad->DrawFrame(-.5, 0, 29.5, max * 1.15,
                        (stageLabel(stage) + ";Layer index (0 = entrance);Transverse std. deviation [mm]").c_str());
        for (size_t i = 0; i < pair.size(); ++i) {
          pair[i]->Draw("LP same");
          if (panel == 1)legend->AddEntry(pair[i], i ? "Simulation" : "Experimental data", "lp");
        }
      }
    layout.header->cd();
    legend->Draw();
    save(c, out, dest, "rms_overview_" + w);
  }
}

std::string closureLabel(const std::string& key) {
  const std::map<std::string, std::string> labels = {
    {"post_layer0_smooth_data", "SSA positions after reconstruction (data)"},
    {"post_layer0_smooth_mc", "SSA positions after reconstruction (simulation)"},
    {"primary_at_layer0", "Simulated electron track at layer 0"},
    {"primary_xy", "Simulated electron at generation vertex"},
    {"beam_input_template", "Input beam-position distribution"}
  };
  return labels.at(key);
}
void closure(const fs::path& base, const std::string& energy, const fs::path& dest, TFile& out) {
  auto d = openRoot(base / "data/diagnostics" / (energy + "GeV") / "diagnostics.root"), m = openRoot(base / "mc/diagnostics" / (energy + "GeV") / "diagnostics.root");
  // Stable TSV labels are kept separately from the descriptive display labels.
  std::vector<std::tuple<TFile*, std::string, std::string, int>>src = {
    {d.get(), "post_layer0_smooth", "data reconstructed", kBlack},
    {m.get(), "post_layer0_smooth", "MC reconstructed", kRed + 1}
  };
  if (param(*m, "has_truth")) {
    src.emplace_back(m.get(), "primary_at_layer0", "MC primary at layer 0", kBlue + 1);
    src.emplace_back(m.get(), "primary_xy", "MC generated vertex", kMagenta + 1);
  }
  if (m->Get("beam_input_template"))src.emplace_back(m.get(), "beam_input_template", "beam template", kGreen + 2);
  TCanvas c("beam_closure", "", 1800, 1100);
  auto layout = page(c, energy + " GeV | Input beam positions and reconstructed positions at layer 0",
                     "All sampled events; distributions normalized separately; SSA positions include +/-2.7 mm smoothing", 2, 1, .28);
  std::vector<std::unique_ptr<TH1D>>keep;
  auto legend = headerLegend(layout, 2, true);
  std::ofstream metrics(dest / "beam_closure.tsv");
  metrics << "series\taxis\tin_range_weight\tmean_mm\trms_mm\tunderflow\toverflow\n";
  for (int dimension = 0; dimension < 2; ++dimension) {
    layout.body->cd(dimension + 1);
    panelStyle();
    gPad->SetLogy();
    bool first = true;
    for (auto& item : src) {
      auto* h = dynamic_cast<TH2*>(std::get<0>(item)->Get(std::get<1>(item).c_str()));
      require(h, "Missing closure map");
      auto p = std::unique_ptr<TH1D>(dimension ? h->ProjectionY((std::get<1>(item) + std::to_string(keep.size())).c_str(), 0, -1) : h->ProjectionX((std::get<1>(item) + std::to_string(keep.size())).c_str(), 0, -1));
      p->SetDirectory(nullptr);
      metrics << std::get<2>(item) << '\t' << (dimension ? "y" : "x") << '\t' << p->Integral() << '\t' << p->GetMean() << '\t' << p->GetStdDev() << '\t' << p->GetBinContent(0) << '\t' << p->GetBinContent(p->GetNbinsX() + 1) << '\n';
      normalize(*p);
      style(*p, std::get<3>(item), first ? 1 : 2);
      p->SetTitle((std::string(dimension ? "Y" : "X") + " direction;" + positionAxis(dimension, "lab") + ";Fraction / bin").c_str());
      p->SetMaximum(1);
      p->SetMinimum(1e-6);
      p->Draw(first ? "hist" : "hist same");
      first = false;
      auto key = std::get<1>(item);
      if (key == "post_layer0_smooth")key += std::get<0>(item) == d.get() ? "_data" : "_mc";
      if (dimension == 0)legend->AddEntry(p.get(), closureLabel(key).c_str(), "l");
      keep.push_back(std::move(p));
    }
  }
  layout.header->cd();
  legend->Draw();
  save(c, out, dest, "beam_closure");
  for (auto sample : {
         "data", "mc"
       }) {
    TFile* f = std::string(sample) == "data" ? d.get() : m.get();
    TCanvas diag("diagnostics", "", 2100, 1250);
    auto diagnosticLayout = page(diag, energy + " GeV | " + sampleLabel(sample) + " | SSA position and signal diagnostics",
                                 "Distance from entrance centroid: X in odd layers, Y in even layers; signal includes negative values", 3, 2, .12);
    std::vector<std::unique_ptr<TH2D>>hs;
    std::vector<std::string>names = {"strips_by_absu", "split_by_absu", "unsplit_by_absu", "pre_energy_vs_absu", "post_energy_vs_absu", "axis_residual"};
    const std::vector<std::string>titles = {
      "Before SSA: strip counts", "After SSA: split positions", "After SSA: unsplit positions",
      "Before SSA: signed signal", "After SSA: signed signal", "Entrance centroid minus electron track"
    };
    for (size_t i = 0; i < names.size(); ++i) {
      diagnosticLayout.body->cd(i + 1);
      panelStyle(true);
      auto h = two(*f, names[i]);
      h->SetName((std::string(sample) + "_" + names[i]).c_str());
      std::string title = titles[i];
      if (i < 3)title += ";Layer index (0 = entrance);Absolute transverse offset [mm];Position counts";
      else if (i < 5)title += ";Absolute transverse offset [mm];Signed signal [MeV];Position counts";
      else {
        if (std::string(sample) == "data")title = "Electron-track truth unavailable for data";
        title += ";Centroid X - track X [mm];Centroid Y - track Y [mm];Events";
      }
      h->SetTitle(title.c_str());
      h->GetZaxis()->SetTitleSize(.04);
      h->GetZaxis()->SetLabelSize(.032);
      h->GetZaxis()->SetTitleOffset(1.7);
      h->Draw("colz");
      hs.push_back(std::move(h));
    }
    save(diag, out, dest, std::string(sample) + "_diagnostics");
  }
}

int main(int argc, char** argv) {
  try {
    require(argc == 3, "Usage: Plot RESULT_BASE ENERGY");
    gROOT->SetBatch();
    TH1::AddDirectory(false);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetTitleFont(43, "t");
    gStyle->SetTitleFontSize(0);
    gStyle->SetTitleSize(23, "t");
    gStyle->SetTitleX(.5);
    gStyle->SetTitleY(.98);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleW(.96);
    gStyle->SetTitleH(.09);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    fs::path base = argv[1];
    std::string energy = argv[2];
    for (auto& s : stages)for (auto& fr : frames)for (auto& w : weights)profiles(base, energy, s, fr, w);
    auto dest = base / "comparison/overview" / (energy + "GeV");
    fs::create_directories(dest);
    auto out = openRoot(dest / "comparison.root", "RECREATE");
    widthOverview(base, energy, dest, *out);
    closure(base, energy, dest, *out);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
  }
}
