#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLine.h"
#include "TTree.h"

namespace {

struct FitParam {
  double landauMPV = 0.0;
  double landauWidth = 0.0;
  double gausSigma = 0.0;
  double chi2ndf = 0.0;
  double chnEntries = 0.0;
};

bool loadReferenceTree(const char* path, std::map<int, FitParam>& outMap) {
  TFile* file = TFile::Open(path, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[ERROR] failed to open reference file: " << path << std::endl;
    return false;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("T_Event"));
  if (!tree) {
    std::cerr << "[ERROR] TTree 'T_Event' was not found in: " << path << std::endl;
    file->Close();
    delete file;
    return false;
  }

  int cellID = 0;
  int ndf = 0;
  int chnEntries = 0;
  double landauMPV = 0.0;
  double landauWidth = 0.0;
  double gausSigma = 0.0;
  double chiSqr = 0.0;

  tree->SetBranchAddress("cellIDs", &cellID);
  tree->SetBranchAddress("landauMPV", &landauMPV);
  tree->SetBranchAddress("landauWidth", &landauWidth);
  tree->SetBranchAddress("gausSigma", &gausSigma);
  tree->SetBranchAddress("ChiSqr", &chiSqr);
  tree->SetBranchAddress("ndf", &ndf);
  tree->SetBranchAddress("chnEntries", &chnEntries);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    FitParam p;
    p.landauMPV = landauMPV;
    p.landauWidth = landauWidth;
    p.gausSigma = gausSigma;
    p.chi2ndf = (ndf > 0) ? (chiSqr / ndf) : 0.0;
    p.chnEntries = chnEntries;
    outMap[cellID] = p;
  }

  file->Close();
  delete file;
  return true;
}

bool loadNewTree(const char* path, std::map<int, FitParam>& outMap) {
  TFile* file = TFile::Open(path, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[ERROR] failed to open new file: " << path << std::endl;
    return false;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("MIP_Fit"));
  if (!tree) {
    std::cerr << "[ERROR] TTree 'MIP_Fit' was not found in: " << path << std::endl;
    file->Close();
    delete file;
    return false;
  }

  int cellID = 0;
  int ndf = 0;
  double landauMPV = 0.0;
  double landauWidth = 0.0;
  double gausSigma = 0.0;
  double chiSquare = 0.0;

  tree->SetBranchAddress("CellID", &cellID);
  tree->SetBranchAddress("LandauMPV", &landauMPV);
  tree->SetBranchAddress("LandauWidth", &landauWidth);
  tree->SetBranchAddress("GauSigma", &gausSigma);
  tree->SetBranchAddress("ChiSquare", &chiSquare);
  tree->SetBranchAddress("NDF", &ndf);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    FitParam p;
    p.landauMPV = landauMPV;
    p.landauWidth = landauWidth;
    p.gausSigma = gausSigma;
    p.chi2ndf = (ndf > 0) ? (chiSquare / ndf) : 0.0;
    p.chnEntries = std::numeric_limits<double>::quiet_NaN();  // This tree has no chnEntries branch.
    outMap[cellID] = p;
  }

  file->Close();
  delete file;
  return true;
}

void drawScatterComparison(const std::string& name, const std::string& axisTitle,
                           const std::vector<double>& refValues, const std::vector<double>& simValues,
                           const std::string& outPrefix, TFile* outFile) {
  if (refValues.empty() || simValues.empty()) return;
  if (refValues.size() != simValues.size()) return;

  std::vector<double> x;
  std::vector<double> y;
  x.reserve(simValues.size());
  y.reserve(refValues.size());
  for (size_t i = 0; i < refValues.size(); ++i) {
    if (!std::isfinite(refValues[i]) || !std::isfinite(simValues[i])) continue;
    x.push_back(simValues[i]);  // x-axis: simulation
    y.push_back(refValues[i]);  // y-axis: reference
  }

  TCanvas* c = new TCanvas(Form("c_%s", name.c_str()), Form("%s scatter", name.c_str()), 900, 850);
  outFile->cd();

  if (x.empty()) {
    c->cd();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.12, 0.55, Form("%s: no valid points to draw", name.c_str()));
    latex.DrawLatex(0.12, 0.48, "simulation side may not contain this branch");
    c->Write();
    c->SaveAs(Form("%s_%s.png", outPrefix.c_str(), name.c_str()));
    return;
  }

  TGraph* g = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
  g->SetTitle(Form("%s scatter;simulation %s;reference %s", name.c_str(), axisTitle.c_str(), axisTitle.c_str()));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.65);
  g->SetMarkerColor(kBlue + 1);
  g->SetLineColor(kBlue + 1);
  g->Draw("AP");

  double xMin = *std::min_element(x.begin(), x.end());
  double xMax = *std::max_element(x.begin(), x.end());
  double yMin = *std::min_element(y.begin(), y.end());
  double yMax = *std::max_element(y.begin(), y.end());

  double lo = std::min(xMin, yMin);
  double hi = std::max(xMax, yMax);
  if (lo == hi) {
    lo -= 1.0;
    hi += 1.0;
  }
  g->GetXaxis()->SetLimits(lo, hi);
  g->SetMinimum(lo);
  g->SetMaximum(hi);

  TLine* diag = new TLine(lo, lo, hi, hi);
  diag->SetLineColor(kRed + 1);
  diag->SetLineStyle(2);
  diag->SetLineWidth(2);
  diag->Draw("same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.14, 0.88, Form("N(common channels) = %zu", x.size()));

  c->Write();
  c->SaveAs(Form("%s_%s.png", outPrefix.c_str(), name.c_str()));
}

}  // namespace

void compare_mip_fit_parameters(
    const char* refPath = "/home/murata_t/ScECAL_BeamTest/analysis/result/mip/MIP_simultaion_new_.root",
    const char* newPath =
        "/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/"
        "all_auto_muon_v4_trackfit_range0.7.root",
    const char* outRootPath = "compare_mip_fit_parameters.root",
    const char* outPrefix = "compare_mip_fit_parameters") {
  std::map<int, FitParam> refMap;
  std::map<int, FitParam> newMap;

  if (!loadReferenceTree(refPath, refMap)) return;
  if (!loadNewTree(newPath, newMap)) return;

  std::vector<int> commonCellIDs;
  commonCellIDs.reserve(refMap.size());
  for (std::map<int, FitParam>::const_iterator it = refMap.begin(); it != refMap.end(); ++it) {
    if (newMap.find(it->first) != newMap.end()) {
      commonCellIDs.push_back(it->first);
    }
  }

  if (commonCellIDs.empty()) {
    std::cerr << "[ERROR] no common CellID found between the two files." << std::endl;
    return;
  }
  std::sort(commonCellIDs.begin(), commonCellIDs.end());

  std::vector<double> ref_landauMPV, new_landauMPV;
  std::vector<double> ref_landauWidth, new_landauWidth;
  std::vector<double> ref_gausSigma, new_gausSigma;
  std::vector<double> ref_chi2ndf, new_chi2ndf;
  std::vector<double> ref_chnEntries, new_chnEntries;

  ref_landauMPV.reserve(commonCellIDs.size());
  new_landauMPV.reserve(commonCellIDs.size());
  ref_landauWidth.reserve(commonCellIDs.size());
  new_landauWidth.reserve(commonCellIDs.size());
  ref_gausSigma.reserve(commonCellIDs.size());
  new_gausSigma.reserve(commonCellIDs.size());
  ref_chi2ndf.reserve(commonCellIDs.size());
  new_chi2ndf.reserve(commonCellIDs.size());
  ref_chnEntries.reserve(commonCellIDs.size());
  new_chnEntries.reserve(commonCellIDs.size());

  for (size_t i = 0; i < commonCellIDs.size(); ++i) {
    const int cellID = commonCellIDs[i];
    const FitParam& a = refMap[cellID];
    const FitParam& b = newMap[cellID];

    ref_landauMPV.push_back(a.landauMPV);
    new_landauMPV.push_back(b.landauMPV);

    ref_landauWidth.push_back(a.landauWidth);
    new_landauWidth.push_back(b.landauWidth);

    ref_gausSigma.push_back(a.gausSigma);
    new_gausSigma.push_back(b.gausSigma);

    ref_chi2ndf.push_back(a.chi2ndf);
    new_chi2ndf.push_back(b.chi2ndf);

    ref_chnEntries.push_back(a.chnEntries);
    new_chnEntries.push_back(b.chnEntries);  // 0 if unavailable in new file
  }

  TFile* outFile = TFile::Open(outRootPath, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "[ERROR] failed to create output ROOT file: " << outRootPath << std::endl;
    if (outFile) delete outFile;
    return;
  }

  drawScatterComparison("landauMPV", "landauMPV", ref_landauMPV, new_landauMPV, outPrefix, outFile);
  drawScatterComparison("landauWidth", "landauWidth", ref_landauWidth, new_landauWidth, outPrefix, outFile);
  drawScatterComparison("gausSigma", "gausSigma", ref_gausSigma, new_gausSigma, outPrefix, outFile);
  drawScatterComparison("ChiSqrOverNdf", "ChiSqr/ndf", ref_chi2ndf, new_chi2ndf, outPrefix, outFile);
  drawScatterComparison("chnEntries", "chnEntries", ref_chnEntries, new_chnEntries, outPrefix, outFile);

  outFile->Write();
  outFile->Close();
  delete outFile;

  std::cout << "[INFO] comparison done." << std::endl;
  std::cout << "[INFO] common CellID count: " << commonCellIDs.size() << std::endl;
  std::cout << "[INFO] output ROOT: " << outRootPath << std::endl;
  std::cout << "[INFO] output PNG prefix: " << outPrefix << "_*.png" << std::endl;
  std::cout << "[NOTE] chnEntries does not exist in MIP_Fit, so chnEntries scatter has no valid points." << std::endl;
}
