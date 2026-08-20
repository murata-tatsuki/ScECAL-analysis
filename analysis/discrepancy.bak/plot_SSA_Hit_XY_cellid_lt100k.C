#include <algorithm>
#include <cstdio>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TTree.h"

// ========== ユーザー設定（ここを編集） ==========
/** 分析対象: ディレクトリなら直下の *.root をすべて、ファイルならその1本を使う */
static const char* kPlotSSA_Hit_XY_DefaultPath = "/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps/ssa/e-/20GeV";
/**
 * true: マクロを読み込んだ直後に kPlotSSA_Hit_XY_DefaultPath で自動実行
 *       （root plot_SSA_Hit_XY_cellid_lt100k.C 用。 .L で定義だけ欲しいときは false）
 */
static const bool kPlotSSA_Hit_XY_AutoRunOnLoad = true;
// ==============================================

namespace {

bool pathIsDirectory(const char* path) {
  if (!path || path[0] == '\0') {
    return false;
  }
  struct stat st;
  if (stat(path, &st) != 0) {
    return false;
  }
  return S_ISDIR(st.st_mode);
}

/** dirPath 直下の通常ファイルのうち名前が .root で終わるものをフルパスで列挙（ソート済み） */
void collectRootFilesInDirectory(const std::string& dirPath, std::vector<std::string>& out) {
  out.clear();
  DIR* dir = opendir(dirPath.c_str());
  if (!dir) {
    std::cerr << "[ERROR] opendir failed: " << dirPath << std::endl;
    return;
  }

  std::string base = dirPath;
  if (!base.empty() && base.back() != '/') {
    base += '/';
  }

  struct dirent* ent = nullptr;
  while ((ent = readdir(dir)) != nullptr) {
    const char* name = ent->d_name;
    if (name[0] == '.' && name[1] == '\0') {
      continue;
    }
    if (name[0] == '.' && name[1] == '.' && name[2] == '\0') {
      continue;
    }
    std::string fname(name);
    const size_t len = fname.size();
    if (len < 5 || fname.compare(len - 5, 5, ".root") != 0) {
      continue;
    }

    const std::string full = base + fname;
    struct stat st;
    if (stat(full.c_str(), &st) != 0 || S_ISDIR(st.st_mode)) {
      continue;
    }
    out.push_back(full);
  }
  closedir(dir);
  std::sort(out.begin(), out.end());
}

}  // namespace

/**
 * SSA_Hit ツリーの各ヒットについて ssaTag==1 かつ CellID < maxCellIDExclusive のものだけを
 * (Hit_X, Hit_Y) で1つの2次元ヒストグラムに詰めて表示する（複数 ROOT ファイル対応）。
 * 同一イベント内: 2D ヒストグラムの同じビンには2回以上 Fill しない（初回のみ）。
 * h2 に加え、全ビン内容の総和が 10000 になるようスケールした TH2 も別キャンバスに描画する。
 *
 * ROOT 上での例:
 *   root -l
 *   .L plot_SSA_Hit_XY_cellid_lt100k.C
 *   // ディレクトリ指定: 直下の *.root をすべて読む
 *   plot_SSA_Hit_XY_cellid_lt100k("/path/to/data_dir");
 *   // 単一ファイル（ディレクトリでなければファイルとして開く）
 *   plot_SSA_Hit_XY_cellid_lt100k("ECAL_Run97_20230430_002956.root");
 *   // ファイルを明示的に列挙
 *   std::vector<std::string> fs = {"a.root", "b.root"};
 *   plot_SSA_Hit_XY_cellid_lt100k(fs);
 *   // ファイル先頭の kPlotSSA_Hit_XY_DefaultPath を使う
 *   plot_SSA_Hit_XY_cellid_lt100k();
 *   // スクリプトをそのまま起動（先頭でパスを設定し AutoRunOnLoad=true のとき）
 *   root plot_SSA_Hit_XY_cellid_lt100k.C
 */
void plot_SSA_Hit_XY_cellid_lt100k(
    const std::vector<std::string>& rootPaths,
    const char* treeName = "SSA_Hit",
    int maxCellIDExclusive = 100000,
    int nxbins = 200,
    double xmin = -2500.0,
    double xmax = 2500.0,
    int nybins = 200,
    double ymin = -2500.0,
    double ymax = 2500.0) {
  if (rootPaths.empty()) {
    std::cerr << "[ERROR] no input ROOT files (empty path list)." << std::endl;
    return;
  }

  char htitle[384];
  std::snprintf(htitle,
                sizeof(htitle),
                "SSA_Hit Hit_X vs Hit_Y (ssaTag==1, CellID < %d; %zu file(s));Hit_X;Hit_Y",
                maxCellIDExclusive,
                rootPaths.size());

  TH2D* h2 = new TH2D(
      "hHitXY",
      "htitle",
      nxbins,
      xmin,
      xmax,
      nybins,
      ymin,
      ymax);

  Long64_t totalEntries = 0;
  Long64_t totalFilled = 0;
  int nFilesOk = 0;

  for (size_t fi = 0; fi < rootPaths.size(); ++fi) {
    const std::string& rootPath = rootPaths[fi];
    if (rootPath.empty()) {
      std::cerr << "[WARN] skip empty path at index " << fi << std::endl;
      continue;
    }

    TFile* file = TFile::Open(rootPath.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "[ERROR] failed to open: " << rootPath << std::endl;
      continue;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));
    if (!tree) {
      std::cerr << "[ERROR] TTree '" << treeName << "' not found in: " << rootPath << std::endl;
      file->Close();
      delete file;
      continue;
    }

    ++nFilesOk;

    std::vector<int>* cellID = nullptr;
    std::vector<int>* ssaTag = nullptr;
    std::vector<double>* hitX = nullptr;
    std::vector<double>* hitY = nullptr;

    tree->SetBranchAddress("CellID", &cellID);
    tree->SetBranchAddress("ssaTag", &ssaTag);
    tree->SetBranchAddress("Hit_X", &hitX);
    tree->SetBranchAddress("Hit_Y", &hitY);

    const Long64_t nEntries = tree->GetEntries();
    totalEntries += nEntries;

    for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      if (!cellID || !ssaTag || !hitX || !hitY) {
        continue;
      }
      const size_t n = cellID->size();
      if (ssaTag->size() != n || hitX->size() != n || hitY->size() != n) {
        std::cerr << "[WARN] " << rootPath << " entry " << i
                  << ": vector size mismatch (CellID=" << n << ", ssaTag=" << ssaTag->size()
                  << ", Hit_X=" << hitX->size() << ", Hit_Y=" << hitY->size() << "), skipped."
                  << std::endl;
        continue;
      }
      std::set<int> binsFilledThisEvent;
      for (size_t j = 0; j < n; ++j) {
        if ((*ssaTag)[j] != 1) {
          continue;
        }
        if ((*cellID)[j] >= maxCellIDExclusive) {
          continue;
        }
        const double x = (*hitX)[j];
        const double y = (*hitY)[j];
        const int gbin = h2->FindBin(x, y);
        if (!binsFilledThisEvent.insert(gbin).second) {
          continue;
        }
        h2->Fill(x, y);
        ++totalFilled;
      }
    }

    file->Close();
    delete file;
  }

  if (nFilesOk == 0) {
    std::cerr << "[ERROR] no file was read successfully; histogram not drawn." << std::endl;
    delete h2;
    return;
  }

  gStyle->SetOptStat(0);

  constexpr double kNormEntriesSum = 100000.0;
  constexpr double scaleFactor = kNormEntriesSum/97864;

  TCanvas* c = new TCanvas("cHitXY", "Hit XY", 900, 800);
  c->SetRightMargin(0.15);
  h2->SetXTitle("Hit_X [mm]");
  h2->SetYTitle("Hit_Y [mm]");
  h2->Draw("COLZ");

  TH2D* h2norm = dynamic_cast<TH2D*>(h2->Clone("hHitXY_norm"));
  if (h2norm) {
    const double sumEntries = h2->Integral();
    if (sumEntries > 0.0) {
      h2norm->Scale(kNormEntriesSum*scaleFactor / sumEntries);
    } else {
      std::cerr << "[WARN] h2 integral is 0; normalized TH2 is identical (no scaling)." << std::endl;
    }
    {
      std::string normTitle = h2->GetTitle();
      normTitle += " [bin contents scaled; sum=";
      normTitle += std::to_string(static_cast<int>(kNormEntriesSum));
      normTitle += "]";
      h2norm->SetTitle(normTitle.c_str());
    }
    h2norm->SetXTitle("Hit_X [mm]");
    h2norm->SetYTitle("Hit_Y [mm]");
    h2norm->SetZTitle("entries (normalized)");

    {
      const int nx = h2norm->GetNbinsX();
      const int ny = h2norm->GetNbinsY();
      std::vector<double> xCenters;
      std::vector<double> yCenters;
      std::vector<double> contents;
      xCenters.reserve(static_cast<size_t>(nx * ny));
      yCenters.reserve(static_cast<size_t>(nx * ny));
      contents.reserve(static_cast<size_t>(nx * ny));
      for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
          const double bc = h2norm->GetBinContent(ix, iy);
          if (bc < 10.0*scaleFactor) {
            continue;
          }
          xCenters.push_back(h2norm->GetXaxis()->GetBinCenter(ix));
          yCenters.push_back(h2norm->GetYaxis()->GetBinCenter(iy));
          contents.push_back(bc);
        }
      }
      const std::streamsize oldPrec = std::cout.precision();
      const std::ios::fmtflags oldFlags = std::cout.flags();
      std::cout << std::fixed << std::setprecision(2);
      for (size_t k = 0; k < xCenters.size(); ++k) {
        if (k > 0) {
          std::cout << ' ';
        }
        std::cout << xCenters[k];
      }
      std::cout << '\n';
      for (size_t k = 0; k < yCenters.size(); ++k) {
        if (k > 0) {
          std::cout << ' ';
        }
        std::cout << yCenters[k];
      }
      std::cout << '\n';
      for (size_t k = 0; k < contents.size(); ++k) {
        if (k > 0) {
          std::cout << ' ';
        }
        std::cout << (int)round(contents[k]);
      }
      std::cout << '\n';
      std::cout.precision(oldPrec);
      std::cout.flags(oldFlags);
    }

    TCanvas* cNorm = new TCanvas("cHitXY_norm", "Hit XY (normalized)", 900, 800);
    cNorm->SetRightMargin(0.15);
    h2norm->Draw("COLZ");
  }

  std::cout << "[INFO] files read OK: " << nFilesOk << " / " << rootPaths.size() << std::endl;
  std::cout << "[INFO] tree entries processed (sum over files): " << totalEntries << std::endl;
  std::cout << "[INFO] fills (ssaTag==1, per event: at most once per 2D bin; CellID < "
            << maxCellIDExclusive << "): " << totalFilled << std::endl;
  if (h2norm) {
    std::cout << "[INFO] h2 Integral() (sum of bin contents before norm): " << h2->Integral()
              << " -> scaled TH2 hHitXY_norm sum=" << kNormEntriesSum << std::endl;
  }
}

/**
 * パスがディレクトリならその直下の *.root をすべて対象にし、
 * そうでなければ1ファイルとして扱う。
 */
void plot_SSA_Hit_XY_cellid_lt100k(
    const char* rootPath,
    const char* treeName = "SSA_Hit",
    int maxCellIDExclusive = 100000,
    int nxbins = 42,
    double xmin = -111.3,
    double xmax = 111.3,
    int nybins = 42,
    double ymin = -111.3,
    double ymax = 111.3
  ) {
  std::vector<std::string> paths;
  if (!rootPath || rootPath[0] == '\0') {
    plot_SSA_Hit_XY_cellid_lt100k(
        paths, treeName, maxCellIDExclusive, nxbins, xmin, xmax, nybins, ymin, ymax);
    return;
  }

  if (pathIsDirectory(rootPath)) {
    collectRootFilesInDirectory(rootPath, paths);
    if (paths.empty()) {
      std::cerr << "[ERROR] no .root files found in directory: " << rootPath << std::endl;
      return;
    }
    std::cout << "[INFO] directory mode: " << paths.size() << " .root file(s) under " << rootPath
              << std::endl;
  } else {
    paths.emplace_back(rootPath);
  }

  plot_SSA_Hit_XY_cellid_lt100k(
      paths, treeName, maxCellIDExclusive, nxbins, xmin, xmax, nybins, ymin, ymax);
}

/** ファイル先頭の kPlotSSA_Hit_XY_DefaultPath を使用 */
void plot_SSA_Hit_XY_cellid_lt100k() {
  if (!kPlotSSA_Hit_XY_DefaultPath || kPlotSSA_Hit_XY_DefaultPath[0] == '\0') {
    std::cerr << "[ERROR] set kPlotSSA_Hit_XY_DefaultPath at the top of "
                 "plot_SSA_Hit_XY_cellid_lt100k.C (or call the const char* overload).\n";
    return;
  }
  plot_SSA_Hit_XY_cellid_lt100k(kPlotSSA_Hit_XY_DefaultPath);
}

namespace {

struct PlotSSAHitXYRunOnLoad {
  PlotSSAHitXYRunOnLoad() {
    if (!kPlotSSA_Hit_XY_AutoRunOnLoad) {
      return;
    }
    if (!kPlotSSA_Hit_XY_DefaultPath || kPlotSSA_Hit_XY_DefaultPath[0] == '\0') {
      std::cerr << "[WARN] kPlotSSA_Hit_XY_DefaultPath is empty; auto-run skipped. "
                   "Edit the macro or call plot_SSA_Hit_XY_cellid_lt100k(\"path\").\n";
      return;
    }
    plot_SSA_Hit_XY_cellid_lt100k();
  }
};

static PlotSSAHitXYRunOnLoad gPlotSSAHitXYRunOnLoad;

}  // namespace
