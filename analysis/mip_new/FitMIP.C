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
    
    TF1*   fitFunc;      // ★追加: 描画・確認用のフィット関数
};

// =================================================================
// 2. Landau x Gaussian 畳み込み関数 (langaufun)
// =================================================================
double langaufun(double *x, double *par) {
    double invsq2pi = 0.3989422804014;
    double mpshift  = -0.22278298;

    double np = 20.0; // 高速化用の畳み込みステップ数
    double sc = 5.0;

    double xx, mpc, fland;
    double sum = 0.0;
    double xlow, xupp, step;

    mpc = par[1] - mpshift * par[0];

    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp - xlow) / np;

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
// 3. ノイズなし結合モデル: Eff(Threshold) * langaufun
// =================================================================
double CombinedFitFunc(double *x, double *par) {
    // par[0..3]: langaufun 用 (Width, MPV, Area, GSigma)
    // par[4..5]: Threshold 効率用 (Threshold, Thresh_Width)

    double langau_val = langaufun(x, par);
    double eff        = 0.5 * (1.0 + TMath::Erf((x[0] - par[4]) / (TMath::Sqrt(2.0) * par[5])));

    return eff * langau_val;
}

// =================================================================
// 4. メインフィッティング関数
// =================================================================
MIPFitResult FitMIPHistogram(TH1* hist,
                             int pix,                      // 0 or 1 (0なら10um, 1なら15um)
                             double threshold,
                             double thresh_width   = 2.0,
                             double fixed_width    = -1.0, // > 0 で指定値を固定
                             double fixed_gsigma   = -1.0, // > 0 で指定値を固定
                             double init_mpv_input = -1.0, 
                             double fit_min        = 0.0,
                             double fit_max        = 1200.0)
{
    MIPFitResult result = {0, 0, 0, 0, 0, 0, -1.0, -1, -1, nullptr};

    if (!hist || hist->GetEntries() == 0) {
        std::cerr << "[Error] Histogram is null or empty!" << std::endl;
        return result;
    }

    // ★ Rebinning の設定
    // 1. センサーごとの「理想的な1ビンのADC幅」を設定
    double target_bin_width = (pix == 0) ? 4.0 : 8.0; 
    // 2. 現在のヒストグラムのビン幅を取得
    double current_bin_width = hist->GetBinWidth(1);
    // 3. 理想のビン幅に近づくような Rebin 係数を計算
    int rebin_factor = std::max(1, (int)std::round(target_bin_width / current_bin_width));
    // 4. 統計量（イベント数）が極端に少ない場合の救済措置
    // 例: エントリー数が 200 未満のスカスカなデータなら、さらに2倍広げて滑らかにする
    if (hist->GetEntries() < 200) {
        rebin_factor *= 2;
    }
    // クローンを作成して実際に Rebin
    TH1* h_fit = (TH1*)hist->Clone(Form("%s_rebin", hist->GetName()));
    h_fit->Rebin(rebin_factor);
    // TF1 オブジェクトの動的生成
    TF1* fitFunc = new TF1(Form("CombinedFitFunc_%s", hist->GetName()), CombinedFitFunc, fit_min, fit_max, 6);
    fitFunc->SetParNames("Width", "MPV", "Area", "GSigma", "Threshold", "Thresh_Width");

    // --- A. Threshold 設定 ---
    fitFunc->FixParameter(4, threshold);
    fitFunc->FixParameter(5, thresh_width);

    // --- B. MIP MPV 初期設定 ---
    double init_mpv;
    if (init_mpv_input > 0.0) {
        init_mpv = init_mpv_input;
    } else {
        init_mpv = (h_fit->GetMean() > threshold) ? h_fit->GetMean() : threshold * 0.8;
    }
    fitFunc->SetParameter(1, init_mpv);
    fitFunc->SetParLimits(1, 0.0, fit_max);

    // --- C. Area (面積) 設定 ---
    // Rebin されたヒストグラムの面積情報を基に初期値を設定
    fitFunc->SetParameter(2, h_fit->Integral() * h_fit->GetBinWidth(1));
    fitFunc->SetParLimits(2, 0.0, h_fit->Integral() * h_fit->GetBinWidth(1) * 10.0);

    // --- D. Landau Width 設定 ---
    if (fixed_width > 0.0) {
        fitFunc->SetParameter(0, fixed_width);
        double width_range_low  = (pix == 0) ? 5.0  : 15.0;
        double width_range_high = (pix == 0) ? 80.0 : 200.0;
        fitFunc->SetParLimits(0, width_range_low, width_range_high);
    } else {
        fitFunc->SetParameter(0, init_mpv * 0.08);
        fitFunc->SetParLimits(0, 1.0, fit_max * 0.2);
    }

    // --- E. Gaussian Sigma 設定 ---
    if (fixed_gsigma > 0.0) {
        fitFunc->SetParameter(3, fixed_gsigma);
        double gaus_range_low  = (pix == 0) ? 5.0  : 10.0;
        double gaus_range_high = (pix == 0) ? 100.0 : 200.0;
        fitFunc->SetParLimits(3, gaus_range_low, gaus_range_high); // ★修正箇所 (3に変更)
    } else {
        fitFunc->SetParameter(3, init_mpv * 0.15);
        fitFunc->SetParLimits(3, 0.1, fit_max * 0.3);
    }

    // --- フィット実行 (Rebinされた h_fit に対して実行) ---
    TFitResultPtr fitResPtr = h_fit->Fit(fitFunc, "Q R S");

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

    // ★ 描画用のスケーリング調整
    // Rebinされたヒストグラムでフィットしたため、1ビンあたりの高さが rebin_factor 倍になっている。
    // 元のヒストグラム(Rebinなし)に重ねて描画するため、Area(面積)を rebin_factor で割る。
    double original_area = fitFunc->GetParameter(2);
    fitFunc->SetParameter(2, original_area / rebin_factor);

    // 元のヒストグラムの描画リストに関数を追加 (後で hist->Draw() した時に自動で赤線が出る)
    hist->GetListOfFunctions()->Add(fitFunc);
    
    // 戻り値の構造体に関数ポインタを格納
    result.fitFunc = fitFunc;

    // 一時的に作成した Rebin 用ヒストグラムをメモリから消去
    delete h_fit;

    return result;
}