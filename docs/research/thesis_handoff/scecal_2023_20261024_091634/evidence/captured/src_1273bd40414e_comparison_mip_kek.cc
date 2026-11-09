#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include "TCanvas.h"
#include "langaus_newConvoluted.C"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;

// ssa後にcutをかけてMPVを求めるマクロ
// pedestal を引いている

// const bool cut = true;  // cutをするかどうかのbool

// const bool fitting = false;




int main(int argc, char* argv[])
{ 
    if(argc < 2){                                                     //エラー処理
        cout << "usage: ./comparison_mip_kek output.root input.root input_1hit.root input_trackFit.root " << endl;
        return 1;
    }

    TFile fileout(argv[1],"RECREATE");

    const int layerNu = 32;
    const int chipNu = 6;
    const int channelNu = 36;
    const int rowNu = 5;
    const int colNu = 42;
    const int doublelayerNu = 2;


    TFile *fin_all = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/kek/fit/allRun.root");
    TFile *fin_1hit = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/kek/fit_1hit/allRun.root");
    // TFile *fin_trackFit = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/kek/fit_trackFit/ECAL_Run171_20250605_162329.root");
    // TFile *fin_trackFit = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/kek/fit_trackFit_layer16/ECAL_Run171_20250605_162329.root");
    TFile *fin_trackFit = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/kek/fit_trackFit_layer16/allRun.root");
    TH1D *MIP_histo_all[layerNu][chipNu][channelNu];
    TH1D *MIP_histo_1hit[layerNu][chipNu][channelNu];
    TH1D *MIP_histo_trackFit[layerNu][chipNu][channelNu];
    fin_all->cd();
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                MIP_histo_all[i][j][k] = (TH1D*)fin_all->Get(Form("chip_canvas/Layer%d_Chip%d_Chn%d",i,j,k));
                if(MIP_histo_all[i][j][k] == nullptr) continue;
            }
        }
    }
    cout << "all hit reading finished" << endl;

    fin_1hit->cd();
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                MIP_histo_1hit[i][j][k] = (TH1D*)fin_1hit->Get(Form("chip_canvas/Layer%d_Chip%d_Chn%d",i,j,k));
                if(MIP_histo_1hit[i][j][k] == nullptr) continue;
            }
        }
    }
    cout << "1hit reading finished" << endl;
    
    fin_trackFit->cd();
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                MIP_histo_trackFit[i][j][k] = (TH1D*)fin_trackFit->Get(Form("chip_canvas/Layer%d_Chip%d_Chn%d",i,j,k));
                if(MIP_histo_trackFit[i][j][k] == nullptr) continue;
            }
        }
    }
    cout << "trackFit reading finished" << endl;




    fileout.cd();   

    char char_tmp[100];
    TCanvas* C_MIP[layerNu][chipNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
        cout << i_layer << "/" << layerNu << endl;
        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
            sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
            C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
            C_MIP[i_layer][i_chip]->Divide(6, 6);

            for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
                if (i_chip == 5 && i_channel > 29) continue;
                if(MIP_histo_all[i_layer][i_chip][i_channel] == nullptr) continue;

                C_MIP[i_layer][i_chip]->cd(i_channel + 1);
                MIP_histo_all[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
                MIP_histo_all[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
                MIP_histo_all[i_layer][i_chip][i_channel]->SetLineColor(kGreen);
                MIP_histo_all[i_layer][i_chip][i_channel]->Rebin(5);
                MIP_histo_all[i_layer][i_chip][i_channel]->Draw("HIST");

                MIP_histo_1hit[i_layer][i_chip][i_channel]->SetLineColor(kBlue);
                MIP_histo_1hit[i_layer][i_chip][i_channel]->Rebin(5);
                MIP_histo_1hit[i_layer][i_chip][i_channel]->Draw("HIST same");

                MIP_histo_trackFit[i_layer][i_chip][i_channel]->SetLineColor(kBlack);
                MIP_histo_trackFit[i_layer][i_chip][i_channel]->Rebin(5);
                MIP_histo_trackFit[i_layer][i_chip][i_channel]->Draw("HIST same");

            }
            C_MIP[i_layer][i_chip]->Write();
        }
    }




}