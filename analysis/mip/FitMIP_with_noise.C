#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TFitResult.h"
#include <iostream>

// =================================================================
// 1. フィット結果を保持する構造体
// =================================================================
struct MIPFitResult {
    double mpv;          // MIP Peak (Landau MPV)
    double mpv_err;      // MIP Peak の誤差
    double width;        // Landau Width
    double width_err;    // Landau Width の誤差
    double gsigma;       // Gaussian Sigma (分解能)
    double gsigma_err;   // Gaussian Sigma の誤差
    
    double chi2;         // Chi2
    int    ndf;          // NDF (自由度)
    int    status;       // ROOT/MINUIT フィットステータス (0: 成功)
};

// =================================================================
// 2. langaus.C から移植・高速化させた Landau x Gaussian 畳み込み関数
// =================================================================
double langaufun(double *x, double *par) {
    // Fit parameters:
    // par[0] = Width (scale) parameter of Landau density
    // par[1] = Most Probable (MP) parameter of Landau density
    // par[2] = Total area
    // par[3] = Width (sigma) of convoluted Gaussian function

    double invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;     // Landau maximum location shift

    // ★高速化のため convolution steps を 100.0 -> 20.0 に変更
    // (6000 チャンネル処理時の CPU 負荷を大幅軽減)
    double np = 20.0; 
    double sc = 5.0;  // convolution extends to +-sc Gaussian sigmas

    double xx, mpc, fland;
    double sum = 0.0;
    double xlow, xupp, step;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp - xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for (double i = 1.0; i <= np / 2.0; i++) {
        xx = xlow + (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);

        xx = xupp - (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

// =================================================================
// 3. 全領域結合モデル: Noise + Eff(Threshold) * langaufun
// =================================================================
double CombinedFitFunc(double *x, double *par) {
    // par[0..3]: langaufun 用 (Width, MP, Area, GSigma)
    // par[4..6]: Pedestal / Noise 用 (Area, Mean, Sigma)
    // par[7..8]: Threshold 効率用 (Threshold, Thresh_Width)

    double langau_val = langaufun(x, par);
    double noise_val  = par[4] * TMath::Gaus(x[0], par[5], par[6], kTRUE);
    double eff        = 0.5 * (1.0 + TMath::Erf((x[0] - par[7]) / (TMath::Sqrt(2.0) * par[8])));

    return noise_val + eff * langau_val;
}

// =================================================================
// 4. メインフィッティング関数
// =================================================================
MIPFitResult FitMIPHistogram(TH1* hist,
                             double threshold,
                             double thresh_width = 2.0,
                             double fixed_width  = -1.0, // > 0 で指定値を固定
                             double fixed_gsigma = -1.0, // > 0 で指定値を固定
                             double fit_min      = 0.0,
                             double fit_max      = 2000.0) // ※不要に大きくしない(1000~2000程度)
{
    MIPFitResult result = {0, 0, 0, 0, 0, 0, -1.0, -1, -1};

    if (!hist || hist->GetEntries() == 0) {
        std::cerr << "[Error] Histogram is null or empty!" << std::endl;
        return result;
    }

    // TF1 オブジェクトの動的生成
    TF1* fitFunc = new TF1("CombinedFitFunc", CombinedFitFunc, fit_min, fit_max, 9);
    fitFunc->SetParNames("Width", "MPV", "Area", "GSigma",
                          "N_Area", "N_Mean", "N_Sigma", "Threshold", "Thresh_Width");

    // --- A. Threshold 設定（expErfThre の結果で固定） ---
    fitFunc->FixParameter(7, threshold);
    fitFunc->FixParameter(8, thresh_width);

    // --- B. Pedestal / Noise 初期設定 ---
    double h_max = hist->GetMaximum();
    fitFunc->SetParameter(4, h_max * 5.0);
    fitFunc->SetParameter(5, 0.0);
    fitFunc->SetParameter(6, threshold * 0.3);

    // --- C. MIP MPV 初期設定（Threshold > MIP Peak に対応） ---
    double init_mpv = (hist->GetMean() > threshold) ? hist->GetMean() : threshold * 0.8;
    fitFunc->SetParameter(1, init_mpv);
    fitFunc->SetParLimits(1, 0.0, fit_max); // 下限を 0 にすることで隠れたピークも探索可能

    fitFunc->SetParameter(2, hist->Integral() * hist->GetBinWidth(1));
    fitFunc->SetParLimits(2, 0.0, hist->Integral() * hist->GetBinWidth(1) * 10.0);

    // --- D. Landau Width 設定（固定 or 自動探索） ---
    if (fixed_width > 0.0) {
        fitFunc->FixParameter(0, fixed_width);
    } else {
        fitFunc->SetParameter(0, init_mpv * 0.08);
        fitFunc->SetParLimits(0, 0.01, fit_max * 0.2);
    }

    // --- E. Gaussian Sigma 設定（固定 or 自動探索） ---
    if (fixed_gsigma > 0.0) {
        fitFunc->FixParameter(3, fixed_gsigma);
    } else {
        fitFunc->SetParameter(3, init_mpv * 0.15);
        fitFunc->SetParLimits(3, 0.01, fit_max * 0.3);
    }

    // --- フィット実行 ---
    // "Q": ログ出力を抑制, "R": フィット範囲指定, "S": 結果ポインタ取得
    TFitResultPtr fitResPtr = hist->Fit(fitFunc, "Q R S");

    // --- 結果の格納 ---
    result.status     = fitResPtr.Get() ? fitResPtr->Status() : -1;
    result.mpv        = fitFunc->GetParameter(1);
    result.mpv_err    = fitFunc->GetParError(1);
    result.width      = fitFunc->GetParameter(0);
    result.width_err  = fitFunc->GetParError(0);
    result.gsigma     = fitFunc->GetParameter(3);
    result.gsigma_err = fitFunc->GetParError(3);
    
    result.chi2       = fitFunc->GetChisquare();
    result.ndf        = fitFunc->GetNDF();

    // ★最重要: TF1 オブジェクトを解放してメモリリークを防止
    // delete fitFunc;

    return result;
}