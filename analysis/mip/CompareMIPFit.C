#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <iostream>
#include <map>
#include <algorithm>

// フィット結果保持用の構造体
struct FitData {
    double mpv;
    double width;
    double gsigma;
    double chi2;
    int    ndf;
    double chi2ndf;
};

// Layer 番号から 15um かどうか判定する関数
bool Is15umLayer(int layer) {
    return (layer <= 3 || layer >= 28);
}

void CompareMIPFit() {
    // 描画スタイルの設定
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // =================================================================
    // 1. ROOT ファイルおよび Tree の読み込み
    // =================================================================
    TFile *f_old = TFile::Open("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/all_auto_muon_v4_trackfit.root");
    TFile *f_new = TFile::Open("/home/murata_t/ScECAL_BeamTest/analysis/result/mip/mip.root");

    if (!f_old || f_old->IsZombie() || !f_new || f_new->IsZombie()) {
        std::cerr << "[Error] Cannot open input ROOT files!" << std::endl;
        return;
    }

    TTree *t_old = (TTree*)f_old->Get("MIP_Fit");
    TTree *t_new = (TTree*)f_new->Get("MIP_Fit");

    if (!t_old || !t_new) {
        std::cerr << "[Error] Cannot find TTree 'MIP_Fit'!" << std::endl;
        return;
    }

    // --- 旧 Tree の Branch 設定 ---
    int old_cellid, old_ndf;
    double old_mpv, old_width, old_gsigma, old_chi2;
    t_old->SetBranchAddress("CellID", &old_cellid);
    t_old->SetBranchAddress("LandauMPV", &old_mpv);
    t_old->SetBranchAddress("LandauWidth", &old_width);
    t_old->SetBranchAddress("GauSigma", &old_gsigma);
    t_old->SetBranchAddress("ChiSquare", &old_chi2);
    t_old->SetBranchAddress("NDF", &old_ndf);

    // --- 新 Tree の Branch 設定 ---
    int new_cellid, new_ndf;
    double new_mpv, new_width, new_gsigma, new_chi2;
    t_new->SetBranchAddress("CellID", &new_cellid);
    t_new->SetBranchAddress("LandauMPV", &new_mpv);
    t_new->SetBranchAddress("LandauWidth", &new_width);
    t_new->SetBranchAddress("GauSigma", &new_gsigma);
    t_new->SetBranchAddress("ChiSquare", &new_chi2);
    t_new->SetBranchAddress("NDF", &new_ndf);

    // =================================================================
    // 2. データロードと CellID マッチング
    // =================================================================
    std::map<int, FitData> map_old;
    std::map<int, FitData> map_new;

    for (Long64_t i = 0; i < t_old->GetEntries(); ++i) {
        t_old->GetEntry(i);
        double c2ndf = (old_ndf > 0) ? old_chi2 / old_ndf : -1.0;
        map_old[old_cellid] = {old_mpv, old_width, old_gsigma, old_chi2, old_ndf, c2ndf};
    }

    for (Long64_t i = 0; i < t_new->GetEntries(); ++i) {
        t_new->GetEntry(i);
        double c2ndf = (new_ndf > 0) ? new_chi2 / new_ndf : -1.0;
        map_new[new_cellid] = {new_mpv, new_width, new_gsigma, new_chi2, new_ndf, c2ndf};
    }

    // =================================================================
    // 3. ヒストグラムの定義 (10um / 15um 分離)
    // =================================================================
    // A. Chi2/NDF
    TH1D *h_c2ndf_old_10 = new TH1D("h_c2ndf_old_10", "#chi^{2}/NDF (10 #mum);#chi^{2}/NDF;Entries", 100, 0, 5);
    TH1D *h_c2ndf_new_10 = new TH1D("h_c2ndf_new_10", "#chi^{2}/NDF (10 #mum);#chi^{2}/NDF;Entries", 100, 0, 5);
    TH1D *h_c2ndf_old_15 = new TH1D("h_c2ndf_old_15", "#chi^{2}/NDF (15 #mum);#chi^{2}/NDF;Entries", 100, 0, 5);
    TH1D *h_c2ndf_new_15 = new TH1D("h_c2ndf_new_15", "#chi^{2}/NDF (15 #mum);#chi^{2}/NDF;Entries", 100, 0, 5);

    TH2D *h2_c2ndf_10 = new TH2D("h2_c2ndf_10", "#chi^{2}/NDF Improvement (10 #mum);Old #chi^{2}/NDF;New #chi^{2}/NDF", 100, 0, 5, 100, 0, 5);
    TH2D *h2_c2ndf_15 = new TH2D("h2_c2ndf_15", "#chi^{2}/NDF Improvement (15 #mum);Old #chi^{2}/NDF;New #chi^{2}/NDF", 100, 0, 5, 100, 0, 5);

    // B. Landau Width
    TH1D *h_width_old_10 = new TH1D("h_width_old_10", "Landau Width (10 #mum);Landau Width [ADC];Entries", 100, 0, 150);
    TH1D *h_width_new_10 = new TH1D("h_width_new_10", "Landau Width (10 #mum);Landau Width [ADC];Entries", 100, 0, 150);
    TH1D *h_width_old_15 = new TH1D("h_width_old_15", "Landau Width (15 #mum);Landau Width [ADC];Entries", 100, 0, 300);
    TH1D *h_width_new_15 = new TH1D("h_width_new_15", "Landau Width (15 #mum);Landau Width [ADC];Entries", 100, 0, 300);

    // C. MPV (1D 重ね合わせ + 2D 相関 + 差分)
    TH1D *h_mpv_old_10 = new TH1D("h_mpv_old_10", "Landau MPV (10 #mum);Landau MPV [ADC];Entries", 100, 0, 300);
    TH1D *h_mpv_new_10 = new TH1D("h_mpv_new_10", "Landau MPV (10 #mum);Landau MPV [ADC];Entries", 100, 0, 300);
    TH1D *h_mpv_old_15 = new TH1D("h_mpv_old_15", "Landau MPV (15 #mum);Landau MPV [ADC];Entries", 100, 0, 1000);
    TH1D *h_mpv_new_15 = new TH1D("h_mpv_new_15", "Landau MPV (15 #mum);Landau MPV [ADC];Entries", 100, 0, 1000);

    TH2D *h2_mpv_10 = new TH2D("h2_mpv_10", "Landau MPV Correlation (10 #mum);Old MPV [ADC];New MPV [ADC]", 100, 0, 300, 100, 0, 300);
    TH2D *h2_mpv_15 = new TH2D("h2_mpv_15", "Landau MPV Correlation (15 #mum);Old MPV [ADC];New MPV [ADC]", 100, 0, 1000, 100, 0, 1000);
    TH1D *h_dmpv_10 = new TH1D("h_dmpv_10", "#Delta MPV (New - Old) (10 #mum);#Delta MPV [ADC];Entries", 100, -50, 50);
    TH1D *h_dmpv_15 = new TH1D("h_dmpv_15", "#Delta MPV (New - Old) (15 #mum);#Delta MPV [ADC];Entries", 100, -100, 100);

    // D. Gaussian Sigma
    TH1D *h_gsigma_old_10 = new TH1D("h_gsigma_old_10", "Gaussian Sigma (10 #mum);GauSigma [ADC];Entries", 100, 0, 100);
    TH1D *h_gsigma_new_10 = new TH1D("h_gsigma_new_10", "Gaussian Sigma (10 #mum);GauSigma [ADC];Entries", 100, 0, 100);
    TH1D *h_gsigma_old_15 = new TH1D("h_gsigma_old_15", "Gaussian Sigma (15 #mum);GauSigma [ADC];Entries", 100, 0, 200);
    TH1D *h_gsigma_new_15 = new TH1D("h_gsigma_new_15", "Gaussian Sigma (15 #mum);GauSigma [ADC];Entries", 100, 0, 200);

    // =================================================================
    // 4. データ Fill
    // =================================================================
    for (auto const& [cellID, d_new] : map_new) {
        int layer = cellID / 100000;
        bool is15 = Is15umLayer(layer);

        if (!is15) {
            h_c2ndf_new_10->Fill(d_new.chi2ndf);
            h_width_new_10->Fill(d_new.width);
            h_mpv_new_10->Fill(d_new.mpv);
            h_gsigma_new_10->Fill(d_new.gsigma);
        } else {
            h_c2ndf_new_15->Fill(d_new.chi2ndf);
            h_width_new_15->Fill(d_new.width);
            h_mpv_new_15->Fill(d_new.mpv);
            h_gsigma_new_15->Fill(d_new.gsigma);
        }

        if (map_old.count(cellID)) {
            auto d_old = map_old[cellID];
            if (!is15) {
                h2_c2ndf_10->Fill(d_old.chi2ndf, d_new.chi2ndf);
                h2_mpv_10->Fill(d_old.mpv, d_new.mpv);
                h_dmpv_10->Fill(d_new.mpv - d_old.mpv);
            } else {
                h2_c2ndf_15->Fill(d_old.chi2ndf, d_new.chi2ndf);
                h2_mpv_15->Fill(d_old.mpv, d_new.mpv);
                h_dmpv_15->Fill(d_new.mpv - d_old.mpv);
            }
        }
    }

    for (auto const& [cellID, d_old] : map_old) {
        int layer = cellID / 100000;
        bool is15 = Is15umLayer(layer);
        if (!is15) {
            h_c2ndf_old_10->Fill(d_old.chi2ndf);
            h_width_old_10->Fill(d_old.width);
            h_mpv_old_10->Fill(d_old.mpv);
            h_gsigma_old_10->Fill(d_old.gsigma);
        } else {
            h_c2ndf_old_15->Fill(d_old.chi2ndf);
            h_width_old_15->Fill(d_old.width);
            h_mpv_old_15->Fill(d_old.mpv);
            h_gsigma_old_15->Fill(d_old.gsigma);
        }
    }

    // =================================================================
    // 5. 描画と PNG 保存処理
    // =================================================================
    auto Set1DStyle = [](TH1D* h_old, TH1D* h_new) {
        h_old->SetLineColor(kRed+1);
        h_old->SetLineWidth(2);
        h_old->SetLineStyle(2); // 破線
        h_new->SetLineColor(kBlue+1);
        h_new->SetLineWidth(2);
        h_new->SetLineStyle(1); // 実線
    };

    Set1DStyle(h_c2ndf_old_10, h_c2ndf_new_10);
    Set1DStyle(h_c2ndf_old_15, h_c2ndf_new_15);
    Set1DStyle(h_width_old_10, h_width_new_10);
    Set1DStyle(h_width_old_15, h_width_new_15);
    Set1DStyle(h_mpv_old_10, h_mpv_new_10);
    Set1DStyle(h_mpv_old_15, h_mpv_new_15);
    Set1DStyle(h_gsigma_old_10, h_gsigma_new_10);
    Set1DStyle(h_gsigma_old_15, h_gsigma_new_15);

    // --- Plot 1: Chi2/NDF 比較 ---
    TCanvas *c1 = new TCanvas("c1", "Chi2/NDF Comparison", 1200, 1000);
    c1->Divide(2, 2);

    c1->cd(1); gPad->SetGrid();
    h_c2ndf_old_10->Draw("HIST");
    h_c2ndf_new_10->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg1->AddEntry(h_c2ndf_old_10, "Old Fit", "l");
    leg1->AddEntry(h_c2ndf_new_10, "New Fit", "l");
    leg1->Draw();

    c1->cd(2); gPad->SetGrid();
    h_c2ndf_old_15->Draw("HIST");
    h_c2ndf_new_15->Draw("HIST SAME");
    leg1->Draw();

    c1->cd(3); gPad->SetGrid();
    h2_c2ndf_10->Draw("COLZ");
    TLine *line1 = new TLine(0, 0, 5, 5);
    line1->SetLineColor(kRed); line1->SetLineStyle(2); line1->Draw();

    c1->cd(4); gPad->SetGrid();
    h2_c2ndf_15->Draw("COLZ");
    line1->Draw();

    c1->SaveAs("../result/mip/figures/compare/Compare_Chi2NDF.png");

    // --- Plot 2: Landau Width 比較 ---
    TCanvas *c2 = new TCanvas("c2", "Landau Width Comparison", 1200, 500);
    c2->Divide(2, 1);

    c2->cd(1); gPad->SetGrid();
    h_width_old_10->Draw("HIST");
    h_width_new_10->Draw("HIST SAME");
    leg1->Draw();

    c2->cd(2); gPad->SetGrid();
    h_width_old_15->Draw("HIST");
    h_width_new_15->Draw("HIST SAME");
    leg1->Draw();

    c2->SaveAs("../result/mip/figures/compare/Compare_LandauWidth.png");

    // --- Plot 3: MPV 相関・差分 ---
    TCanvas *c3 = new TCanvas("c3", "MPV Correlation and Diff", 1200, 1000);
    c3->Divide(2, 2);

    c3->cd(1); gPad->SetGrid();
    h2_mpv_10->Draw("COLZ");
    TLine *line2 = new TLine(0, 0, 300, 300);
    line2->SetLineColor(kRed); line2->SetLineStyle(2); line2->Draw();

    c3->cd(2); gPad->SetGrid();
    h2_mpv_15->Draw("COLZ");
    TLine *line3 = new TLine(0, 0, 1000, 1000);
    line3->SetLineColor(kRed); line3->SetLineStyle(2); line3->Draw();

    c3->cd(3); gPad->SetGrid();
    h_dmpv_10->SetLineColor(kBlue+1); h_dmpv_10->SetLineWidth(2);
    h_dmpv_10->Draw("HIST");

    c3->cd(4); gPad->SetGrid();
    h_dmpv_15->SetLineColor(kBlue+1); h_dmpv_15->SetLineWidth(2);
    h_dmpv_15->Draw("HIST");

    c3->SaveAs("../result/mip/figures/compare/Compare_MPV_Correlation.png");

    // --- Plot 4: GauSigma 比較 (Log-y スケール適用) ---
    TCanvas *c4 = new TCanvas("c4", "Gaussian Sigma Comparison", 1200, 500);
    c4->Divide(2, 1);

    c4->cd(1); gPad->SetGrid(); gPad->SetLogy();
    h_gsigma_old_10->SetMinimum(0.5);
    h_gsigma_old_10->Draw("HIST");
    h_gsigma_new_10->Draw("HIST SAME");
    leg1->Draw();

    c4->cd(2); gPad->SetGrid(); gPad->SetLogy();
    h_gsigma_old_15->SetMinimum(0.5);
    h_gsigma_old_15->Draw("HIST");
    h_gsigma_new_15->Draw("HIST SAME");
    leg1->Draw();

    c4->SaveAs("../result/mip/figures/compare/Compare_GauSigma.png");

    // --- Plot 5: Landau MPV 1D 重ね合わせ比較 ---
    TCanvas *c5 = new TCanvas("c5", "Landau MPV 1D Comparison", 1200, 500);
    c5->Divide(2, 1);

    c5->cd(1); gPad->SetGrid();
    double max_10 = std::max(h_mpv_old_10->GetMaximum(), h_mpv_new_10->GetMaximum());
    h_mpv_old_10->SetMaximum(max_10 * 1.15);
    h_mpv_old_10->Draw("HIST");
    h_mpv_new_10->Draw("HIST SAME");
    leg1->Draw();

    c5->cd(2); gPad->SetGrid();
    double max_15 = std::max(h_mpv_old_15->GetMaximum(), h_mpv_new_15->GetMaximum());
    h_mpv_old_15->SetMaximum(max_15 * 1.15);
    h_mpv_old_15->Draw("HIST");
    h_mpv_new_15->Draw("HIST SAME");
    leg1->Draw();

    c5->SaveAs("../result/mip/figures/compare/Compare_MPV.png");

    std::cout << "[Success] All 5 comparison plots exported successfully as PNG!" << std::endl;
}