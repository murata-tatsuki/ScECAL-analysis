#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include "FastPng.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
struct Options {
    std::string data, mc, output, figures, channels = "all", energy = "100";
    std::string layers = "all", dataOrigin, mcOrigin;
    bool saveCanvases = false;
};
Options options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; i += 2) {
        if (i + 1 == argc) throw std::runtime_error("Missing argument value");
        std::string key = argv[i], value = argv[i + 1];
        if (key == "--data") o.data = value;
        else if (key == "--mc") o.mc = value;
        else if (key == "--output") o.output = value;
        else if (key == "--figure-dir") o.figures = value;
        else if (key == "--channels") o.channels = value;
        else if (key == "--energy") o.energy = value;
        else if (key == "--layers") o.layers = value;
        else if (key == "--data-origin") o.dataOrigin = value;
        else if (key == "--mc-origin") o.mcOrigin = value;
        else if (key == "--save-canvases") {
            if (value != "0" && value != "1") throw std::runtime_error("save-canvases must be 0 or 1");
            o.saveCanvases = value == "1";
        }
        else throw std::runtime_error("Unknown option: " + key);
    }
    if (o.data.empty() || o.mc.empty() || (o.saveCanvases && o.output.empty()) || o.figures.empty() ||
        (o.channels != "all" && !std::regex_match(o.channels, std::regex("[0-9]+(,[0-9]+)*"))))
        throw std::runtime_error("Invalid options; use run.sh --help");
    return o;
}
std::unique_ptr<TFile> openFile(const std::string& name, const char* mode = "READ") {
    std::unique_ptr<TFile> f(TFile::Open(name.c_str(), mode));
    if (!f || f->IsZombie()) throw std::runtime_error("Cannot open ROOT file: " + name);
    return f;
}
template<class H> std::unique_ptr<H> histogram(TFile& file, const std::string& name) {
    auto* h = dynamic_cast<H*>(file.Get(name.c_str()));
    if (h) h->SetDirectory(nullptr);
    return std::unique_ptr<H>(h);
}
void text(double x, double y, const std::string& value, double size) {
    TLatex label;
    label.SetNDC(); label.SetTextFont(42); label.SetTextSize(size);
    label.DrawLatex(x, y, value.c_str());
}
void pad() {
    gPad->SetLeftMargin(.16); gPad->SetRightMargin(.16);
    gPad->SetBottomMargin(.16); gPad->SetTopMargin(.20); gPad->SetTicks(1, 1);
}
bool sameAxes(const TH2D& a, const TH2D& b) {
    return a.GetNbinsX() == b.GetNbinsX() && a.GetNbinsY() == b.GetNbinsY() &&
        a.GetXaxis()->GetXmin() == b.GetXaxis()->GetXmin() && a.GetXaxis()->GetXmax() == b.GetXaxis()->GetXmax() &&
        a.GetYaxis()->GetXmin() == b.GetYaxis()->GetXmin() && a.GetYaxis()->GetXmax() == b.GetYaxis()->GetXmax();
}
TH1F* frame(const TH2D& h, bool profile = false) {
    std::string title = ";" + std::string(h.GetXaxis()->GetTitle()) + ";" +
        (profile ? "Mean " : "") + h.GetYaxis()->GetTitle();
    auto* f = gPad->DrawFrame(h.GetXaxis()->GetXmin(), h.GetYaxis()->GetXmin(),
                             h.GetXaxis()->GetXmax(), h.GetYaxis()->GetXmax(), title.c_str());
    f->GetXaxis()->SetTitleOffset(1.3); f->GetYaxis()->SetTitleOffset(1.75);
    return f;
}
void render(const Options& o) {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetPalette(kViridis);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetEndErrorSize(0); gStyle->SetTitleSize(.043, "XYZ"); gStyle->SetLabelSize(.036, "XYZ");
    TH1::AddDirectory(false);
    std::array<std::unique_ptr<TFile>, 2> files{openFile(o.data), openFile(o.mc)};
    std::set<int> ids;
    const std::regex channelName("channel_[0-9]+");
    for (const auto& f : files) {
        TIter next(f->GetListOfKeys());
        while (auto* key = dynamic_cast<TKey*>(next())) {
            std::string name = key->GetName();
            if (std::regex_match(name, channelName)) ids.insert(std::stoi(name.substr(8)));
        }
    }
    if (o.channels != "all") {
        std::set<int> requested;
        std::istringstream list(o.channels);
        for (std::string item; std::getline(list, item, ',');) {
            int id = std::stoi(item);
            if (!ids.count(id)) throw std::runtime_error("Channel absent from both histogram files: " + item);
            requested.insert(id);
        }
        ids = std::move(requested);
    }
    if (ids.empty()) throw std::runtime_error("No channels to plot");
    if (o.layers != "all") {
        std::set<int> layers;
        std::istringstream list(o.layers);
        for (std::string item; std::getline(list, item, ',');) {
            if (!std::regex_match(item, std::regex("[0-9]+(-[0-9]+)?"))) throw std::runtime_error("Invalid layers");
            auto dash = item.find('-');
            int lo = std::stoi(item), hi = dash == std::string::npos ? lo : std::stoi(item.substr(dash + 1));
            if (lo < 0 || hi > 31 || hi < lo) throw std::runtime_error("Invalid layers");
            for (int l = lo; l <= hi; ++l) layers.insert(l);
        }
        for (auto it = ids.begin(); it != ids.end();) {
            if (!layers.count(*it / 100000)) it = ids.erase(it); else ++it;
        }
    }
    if (ids.empty()) { std::cout << "No observed channels in requested layers; skipped" << std::endl; return; }
    struct PlotSet { std::string name; std::vector<std::pair<std::string, std::string>> rows; };
    const std::vector<PlotSet> groups = {
        {"hg_lg", {{"hg_lg", "HG / LG: full range"}, {"hg_lg_zoom", "HG / LG: low ADC zoom"}}},
        {"adc_energy", {{"hg_energy", "HG / saved HG energy"}, {"lg_energy", "LG / saved LG energy"},
            {"hg_selected_energy", "HG / final energy (all matched hits)"},
            {"lg_selected_energy", "LG / final energy (all matched hits)"}}}
    };
    fs::path output(o.output), temporary(o.output + ".tmp");
    if (output.has_parent_path()) fs::create_directories(output.parent_path());
    // Destroy/close the file before removing an incomplete temporary on error.
    struct Cleanup {
        fs::path path;
        ~Cleanup() { if (!path.empty()) { std::error_code error; fs::remove(path, error); } }
    } cleanup{o.saveCanvases ? temporary : fs::path{}};
    std::unique_ptr<TFile> out;
    if (o.saveCanvases) {
        out = openFile(temporary.string(), "RECREATE");
        TNamed("data_histograms", (o.dataOrigin.empty() ? o.data : o.dataOrigin).c_str()).Write();
        TNamed("simulation_histograms", (o.mcOrigin.empty() ? o.mc : o.mcOrigin).c_str()).Write();
        TNamed("profile_definition", "Exact mean y and standard error for each x bin. All finite y values, including outside the density y display range. No calibration constants refitted.").Write();
    }
    const std::array<const char*, 2> labels{"Data", "Simulation"};
    const std::array<int, 2> colors{kBlue + 1, kOrange + 7}, markers{20, 24};
    size_t count = 0;
    size_t figureCount = 0;
    double rasterSeconds = 0, encodeSeconds = 0;
    for (int cell : ids) {
        int layer = cell / 100000, chip = cell % 100000 / 10000, channel = cell % 100;
        std::string directory = "channel_" + std::to_string(cell);
        std::string caption = o.energy + " GeV  |  layer " + std::to_string(layer) + ", chip " +
            std::to_string(chip) + ", channel " + std::to_string(channel) + "  |  CellID " + std::to_string(cell);
        fs::path dest = fs::path(o.figures) / ("layer" + std::to_string(layer)) / ("chip" + std::to_string(chip));
        fs::create_directories(dest);
        auto* rootDir = out ? out->mkdir(directory.c_str()) : nullptr;
        if (out && !rootDir) throw std::runtime_error("Cannot create output directory: " + directory);
        for (const auto& group : groups) {
            // Pad references must be cleared before these owned objects die.
            // Load/release only one group's dense histograms at a time.
            std::vector<std::unique_ptr<TObject>> keep;
            TCanvas canvas((group.name + "_" + std::to_string(cell)).c_str(), caption.c_str(), 1800, 470 * group.rows.size());
            canvas.Divide(3, group.rows.size(), .004, .007);
            for (size_t row = 0; row < group.rows.size(); ++row) {
                const auto& [name, title] = group.rows[row];
                const std::string prefix = directory + "/" + name;
                std::array<TH2D*, 2> hist{};
                for (size_t col = 0; col < files.size(); ++col) {
                    auto h = histogram<TH2D>(*files[col], prefix);
                    hist[col] = h.get();
                    if (h) keep.push_back(std::move(h));
                }
                const auto* model = hist[0] ? hist[0] : hist[1];
                if (!model) continue;
                double zmax = 2;
                for (const auto* h : hist) {
                    if (!h) continue;
                    if (!sameAxes(*h, *model)) throw std::runtime_error("Incompatible data/MC axes for " + prefix);
                    zmax = std::max(zmax, h->GetMaximum());
                }
                for (size_t col = 0; col < files.size(); ++col) {
                    canvas.cd(3 * row + col + 1); pad();
                    auto* h = hist[col];
                    if (h) {
                        h->GetXaxis()->SetTitleOffset(1.3); h->GetYaxis()->SetTitleOffset(1.75);
                        h->SetMinimum(.5); h->SetMaximum(zmax); gPad->SetLogz(true); h->Draw("COLZ");
                        auto* outside = dynamic_cast<TParameter<Long64_t>*>(files[col]->Get((prefix + "_outside").c_str()));
                        if (!outside) throw std::runtime_error("Missing outside count for " + prefix);
                        text(.16, .855, std::string(labels[col]) + ": N = " + std::to_string(llround(h->GetEntries())) +
                             ", outside = " + std::to_string(outside->GetVal()), .039);
                    } else {
                        frame(*model);
                        text(.22, .48, std::string(labels[col]) + ": no entries for this channel", .039);
                    }
                    text(.16, .95, caption, .031); text(.16, .90, title, .039);
                }
                canvas.cd(3 * row + 3); pad(); gPad->SetRightMargin(.06); gPad->SetGridy(true);
                frame(*model, true);
                auto legend = std::make_unique<TLegend>(.18, .66, .58, .79);
                legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(.037);
                for (size_t col = 0; col < files.size(); ++col) {
                    auto profile = histogram<TProfile>(*files[col], prefix + "_profile");
                    if (!profile) {
                        if (hist[col]) throw std::runtime_error("Missing profile for " + prefix);
                        continue;
                    }
                    auto graph = std::make_unique<TGraphErrors>();
                    for (int bin = 1; bin <= profile->GetNbinsX(); ++bin) {
                        if (profile->GetBinEntries(bin) <= 0) continue;
                        int point = graph->GetN();
                        graph->SetPoint(point, profile->GetBinCenter(bin), profile->GetBinContent(bin));
                        graph->SetPointError(point, 0, profile->GetBinError(bin));
                    }
                    graph->SetMarkerColor(colors[col]); graph->SetLineColor(colors[col]);
                    graph->SetMarkerStyle(markers[col]); graph->SetMarkerSize(.6);
                    graph->Draw("P SAME"); legend->AddEntry(graph.get(), labels[col], "pe");
                    keep.push_back(std::move(graph));
                }
                legend->Draw(); keep.push_back(std::move(legend));
                text(.16, .95, caption, .031);
                text(.16, .90, "Profile: mean and standard error", .039);
                text(.16, .845, "Includes y outside density display range", .032);
            }
            canvas.Modified();
            if (rootDir) {
                canvas.Update();
                rootDir->cd();
                if (canvas.Write(group.name.c_str()) <= 0) throw std::runtime_error("Cannot write comparison canvas");
            }
            auto figure = dest / ("channel" + std::to_string(channel) + "_cell" + std::to_string(cell) + "_" + group.name + ".png");
            const auto timing = fastPng(canvas, figure.c_str());
            rasterSeconds += timing.raster; encodeSeconds += timing.encode; ++figureCount;
            if (!fs::exists(figure) || fs::file_size(figure) == 0) throw std::runtime_error("Cannot save figure: " + figure.string());
            canvas.Close();
        }
        ++count;
        if (count % 50 == 0 || count == ids.size()) std::cout << "Plotted " << count << "/" << ids.size() << " channels" << std::endl;
    }
    if (out) {
        out->Close();
        if (out->TestBit(TFile::kWriteError)) throw std::runtime_error("Error writing comparison ROOT file");
        fs::rename(temporary, output);
    }
    std::cout << "PNG timing: " << figureCount << " figures; raster=" << rasterSeconds
              << " s, encode=" << encodeSeconds << " s" << std::endl;
}
int main(int argc, char** argv) {
    try { render(options(argc, argv)); return 0; }
    catch (const std::exception& error) { std::cerr << "ERROR: " << error.what() << std::endl; return 1; }
}
