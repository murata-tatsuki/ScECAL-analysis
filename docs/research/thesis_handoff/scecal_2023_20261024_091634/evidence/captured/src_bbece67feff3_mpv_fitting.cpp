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
#include "langaus.C"				// langaus(chargeH, &fitFunc, &peakP, &peakPError);
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/src/EBUdecode.cxx"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;


// 作成されたmip spectrumからfitするマクロ
// error function * landau

void mpv_fitting()
{
    const int LayerNo=30;
    TGraph* gEff[LayerNo];
    TGraph* gTemp[LayerNo];
    for(int i=0; i!=LayerNo; ++i)
    {
        gEff[i] = new TGraph();
        gTemp[i] = new TGraph();
    }
    // hit chn=0
   // double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,
   //                           2,1,4,2,1,2,3,3,3,1,
   //                           2,4,1,1,0,0,0,0,9,3};
    // hit chn<=10
    double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,
                              2,1,4,2,1,2,3,3,3,1,
                              2,4,1,2,0,1,0,0,9,3};

    // ifstream filelist("./datalist_CR.txt");
    // int fileNumber = 105;
    // ifstream filelist("./datalist_sps.txt");
    // int fileNumber = 28;
    // ifstream filelist("./datalist_pi_tilt.txt");
    // int fileNumber = 2;
    // if(!filelist.good()){cout<<"open datalist file failed !"<<endl; exit(0);}
    int nbin=0;
    double effmean = 0.;
    double NhitLayer[30] = {0};
    double totalEntry = 0;


    const int layerNu = 32;
    const int chipNu = 6;
    const int channelNu = 36;
    const double _xInterval = 5.3;  // 300 um gap in width direction
    const double _yInterval = 45.4; // 400 um gap in length direction
    const int rowNu = 42;
    const int columnNu = 5;
    const double _zInitialPosition = 1.;
    const double _zIntervalEven = 12.7; // mm
    const double _zIntervalOdd = 19.9; // mm
    double _xMin = -rowNu*_xInterval/2.;
    double _xMax = rowNu*_xInterval/2.;
    double _yMin = -columnNu*_yInterval/2.;
    double _yMax = columnNu*_yInterval/2.;
    double _zMin = -10; double _zMax = 320; int _zBin = 165;
    double coormin = -(rowNu-1)*_xInterval/2.;
    double accurateBinY[10] = {-113.3,-68.3, -67.9,-22.9, -22.5,22.5, 22.9,67.9, 68.3,113.3};
    double accurateBinX[84] = {0};
    for(int i=0;i<rowNu;i++){
        accurateBinX[2*i] = -111.15 + i*_xInterval;
        accurateBinX[2*i+1] = accurateBinX[2*i] + 5.0;
    }
    TH2F *gainMap[30];      // channelごとのgain
    TH2F *mipMap[30];       // channelごとのMIP
    TH2F *eff_2d[30];       // channelごとのdetection efficiency
    TH2F *hit_2d[30];       // channelごとのhit数
    TH2F *fit_2d[30];       // track fit のhit数 (trackが通過しているかどうか)
    TH2F *fit_hit_2d[30];   // trackが通過しているchannelに実際にhitがあった数
    for(int i=0;i<30;i++){
        double xNu_ = i%2==0 ? columnNu : rowNu;
        double yNu_ = i%2==0 ? rowNu : columnNu;
        double xMax_ = i%2==0 ? _yMax : _xMax;
        double yMax_ = i%2==0 ? _xMax : _yMax;

        double accuratexNu_ = xNu_ * 2 - 1;
        double accurateyNu_ = yNu_ * 2 - 1;
        double* accurateXBin = i%2==0 ? accurateBinY : accurateBinX;
        double* accurateYBin = i%2==0 ? accurateBinX : accurateBinY;
        eff_2d[i] = new TH2F(Form("eff_2d_%d",i),Form("efficiency at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        gainMap[i] = new TH2F(Form("gainMap_%d",i),Form("gain at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        mipMap[i] = new TH2F(Form("mipMap_%d",i),Form("mip at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        hit_2d[i] = new TH2F(Form("hit_2d_%d",i),Form("hit map layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        // fit_2d[i] = new TH2F(Form("fit_2d_%d",i),Form("position of fit at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        // fit_hit_2d[i] = new TH2F(Form("fit_hit_2d_%d",i),Form("number of fit match hit at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        fit_2d[i] = new TH2F(Form("fit_2d_%d",i),Form("position of fit at layer %d",i), accuratexNu_,accurateXBin, accurateyNu_,accurateYBin);
        fit_hit_2d[i] = new TH2F(Form("fit_hit_2d_%d",i),Form("number of fit match hit at layer %d",i), accuratexNu_,accurateXBin, accurateyNu_,accurateYBin);
    }
    double NhitLayer_channel[30][6][36] = {0};
    double totalEntry_channel[30][6][36] = {0};
    double layer_z[30][2] = { { 0, 1},{ 1,  12.2},{ 2,  20.9},{ 3,  32.1},{ 4,  40.8},{ 5,    52},{ 6,  60.7},{ 7,  71.9},{ 8,  80.6},{ 9,  91.8},{10, 100.5},
                                {11, 111.7},{12, 120.4},{13, 131.6},{14, 140.3},{15, 151.5},{16, 160.2},{17, 171.4},{18, 180.1},{19, 191.3},{20,   200},
                                {21, 211.2},{22, 219.9},{23, 231.1},{24, 239.8},{25,   251},{26, 259.7},{27, 270.9},{28, 279.6},{29, 290.8} };


  /////////////////////////////////////////////////// mip 
    double MIP[layerNu][chipNu][channelNu];
    double MIP_entry[layerNu][chipNu][channelNu];
    TH1F *MIP_entry_check = new TH1F("MIP_entry_check","MIP_entry_check",1000,0,10000);
    // TGraphErrors *gain_graph = new TGraphErrors();
    // int count_gain_graph = 0;
    int _MIPcellID, _MIPndf, _MIPchnEntries;
    double _landauMPV, _landauWidth, _gausSigma, _ChiSqr;
    // TFile *gainfile = new TFile("/megraid01/users/tsuji/scecal/ScECAL_CR/lightCalib/analyseCode/CaliSPS/macro/SPSgain.root");
    TFile *mipCRfile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/MIP_normalgain_all.root");
    TTree *t_mipCRfile = (TTree*) mipCRfile->Get("T_Event");
    t_mipCRfile->SetBranchAddress("cellIDs",&_MIPcellID);
    t_mipCRfile->SetBranchAddress("chnEntries",&_MIPchnEntries);
    t_mipCRfile->SetBranchAddress("landauMPV",&_landauMPV);
    t_mipCRfile->SetBranchAddress("landauWidth",&_landauWidth);
    t_mipCRfile->SetBranchAddress("gausSigma",&_gausSigma);
    t_mipCRfile->SetBranchAddress("ChiSqr",&_ChiSqr);
    t_mipCRfile->SetBranchAddress("ndf",&_MIPndf);
    int MIPentry = 0;
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                MIP[i][j][k] = 0;
                MIP_entry[i][j][k] = 0;
                if(j==5&&k>29) continue;
                t_mipCRfile->GetEntry(MIPentry);
                MIPentry++;
                int cellid_now = (i)*1e6 + j*10000 + k;
                if(cellid_now!=_MIPcellID){
                    cout << i << ", " << j << ", " << k << "   " << _MIPcellID << " " << cellid_now;
                    cout << "    something is wrong at mip search" << endl;
                }
                if(_MIPchnEntries<200) continue;
                MIP[i][j][k] = _landauMPV;
                MIP_entry[i][j][k] = _MIPchnEntries;
                MIP_entry_check->Fill(_MIPchnEntries);
                //cout << i << ", " << j << ", " << k << "   gain : "  << _SPSgain << "  error : " << _SPSgainError << endl;
                // gain_graph->SetPoint(count_gain_graph,count_gain_graph,_SPSgain);
                // gain_graph->SetPointError(count_gain_graph,0,_SPSgainError);
                // count_gain_graph++;
            }
      }
    }
    mipCRfile->Close();
    cout << "sps mip reading finished" << endl;
  ///////////////////////////////////////////////////

  /////////////////////////////////////////////////// gain 
    double SPSgain[layerNu][chipNu][channelNu];
    double SPSgain_error[layerNu][chipNu][channelNu];
    // TH1F *SPSgain_histo = new TH1F("SPSgain_histo","gain",120,10,40);
    // TGraphErrors *gain_graph = new TGraphErrors();
    // int count_gain_graph = 0;
    int _SPScellID;
    double _SPSgain, _SPSgainError;
    // TFile *gainfile = new TFile("/megraid01/users/tsuji/scecal/ScECAL_CR/lightCalib/analyseCode/CaliSPS/macro/SPSgain.root");
    TFile *gainfile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/led/20230504_1425_Calib_LED_all.root");
    TTree *t_gainfile = (TTree*) gainfile->Get("T_Event");
    t_gainfile->SetBranchAddress("cellIDs",&_SPScellID);
    t_gainfile->SetBranchAddress("gain",&_SPSgain);
    t_gainfile->SetBranchAddress("gainError",&_SPSgainError);
    int SPSentry = 0;
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                t_gainfile->GetEntry(SPSentry);
                int cellid_now = (i)*1e6 + j*10000 + k;
                if(cellid_now!=_SPScellID){
                	cout << i << ", " << j << ", " << k << "   " << _SPScellID << " " << cellid_now;
                	cout << "    something is wrong at gain search" << endl;
                }
                SPSgain[i][j][k] = _SPSgain;
                SPSgain_error[i][j][k] = _SPSgainError;
                SPSentry++;
                //cout << i << ", " << j << ", " << k << "   gain : "  << _SPSgain << "  error : " << _SPSgainError << endl;
                // SPSgain_histo->Fill(_SPSgain);
                // gain_graph->SetPoint(count_gain_graph,count_gain_graph,_SPSgain);
                // gain_graph->SetPointError(count_gain_graph,0,_SPSgainError);
                // count_gain_graph++;
            }
        }
    }
    gainfile->Close();
    cout << "sps gain reading finished" << endl;
  ///////////////////////////////////////////////////

  /////////////////////////////////////////////////// pedestal
    TFile *fileinPed;
    fileinPed = new TFile("/megraid01/users/data_beamtest/ECAL_data/analysed_oldFormat/2023/ps/ped/fit/pedestal/ECAL_Run23_20230518_153022.root");
    TTree *treePed = (TTree*) fileinPed->Get("ChnLevel");
    int entry_max_Ped = treePed->GetEntries();
    vector<int> *cellIDs_ped = nullptr;
    vector<float> *pedMeanTimes = nullptr, *pedSigmaTimes = nullptr, *pedMeanCharges = nullptr, *pedSigmaCharges = nullptr;
    double HG_ped[layerNu][chipNu][channelNu];
    double HG_pedSigma[layerNu][chipNu][channelNu];
    double LG_ped[layerNu][chipNu][channelNu];
    double LG_pedSigma[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
            for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
                HG_ped[i_layer][i_chip][i_channel] = 0;
                HG_pedSigma[i_layer][i_chip][i_channel] = 0;
                LG_ped[i_layer][i_chip][i_channel] = 0;
                LG_pedSigma[i_layer][i_chip][i_channel] = 0;
            }
        }
    }
    treePed->SetBranchAddress("cellIDs", &cellIDs_ped);
    treePed->SetBranchAddress("pedMeanTimes", &pedMeanTimes);
    treePed->SetBranchAddress("pedSigmaTimes", &pedSigmaTimes);
    treePed->SetBranchAddress("pedMeanCharges", &pedMeanCharges);
    treePed->SetBranchAddress("pedSigmaCharges", &pedSigmaCharges);
    for(int ientry=0; ientry<entry_max_Ped; ientry++){
        treePed->GetEntry(ientry);
        for(int i=0; i<(int)cellIDs_ped->size(); i++){
            int _layer = cellIDs_ped->at(i)/1000000;
            int _chip = (cellIDs_ped->at(i)/10000)%100;
            int _channel = cellIDs_ped->at(i)%100;

            HG_ped[_layer][_chip][_channel] = pedMeanTimes->at(i);
            HG_pedSigma[_layer][_chip][_channel] = pedSigmaTimes->at(i);
            LG_ped[_layer][_chip][_channel] = pedMeanCharges->at(i);
            LG_pedSigma[_layer][_chip][_channel] = pedSigmaCharges->at(i);
        }
    }
    fileinPed->Close();
    cout << "sps pedestal reading finished" << endl;
  ///////////////////////////////////////////////////

  /////////////////////////////////////////////////// threshold
    const double constantFraction = 0.05;
    const int nbin_above_pedestal = 0;
    TF1 *f1 = new TF1("f1","(1+TMath::Erf((x-[0])/[1]))/2 * landau(2)",0,2000);
    double channel_thre[layerNu][chipNu][channelNu];
    TH1F *mip_thre_histo[layerNu][chipNu][channelNu];
    TCanvas *mipthrecanvas = new TCanvas("mipthrecanvas","mipthrecanvas",1);
    mipthrecanvas->Divide(6, 6);
    const bool cehck_canvas_mip = false;
    TH1F *threshold_check = new TH1F("threshold_check","threshold_check",100,-50,350);
    // TH1F *SPSgain_histo = new TH1F("SPSgain_histo","gain",120,10,40);
    // TGraphErrors *gain_graph = new TGraphErrors();
    // int count_gain_graph = 0;
    int _threcellID;
    TH1F *_thre_histo = new TH1F();
    TH1F *_thre_histo_tempCorrected = new TH1F();
    TFile *threfile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip/withPedestal/MIP_gain_all.root");
    TTree *t_threfile = (TTree*) threfile->Get("T_Event");
    t_threfile->SetBranchAddress("cellIDs",&_threcellID);
    t_threfile->SetBranchAddress("histo",&_thre_histo);
    t_threfile->SetBranchAddress("histo_tempCorrected",&_thre_histo_tempCorrected);
    int threEntry = 0;
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                channel_thre[i][j][k] = 0;
                t_threfile->GetEntry(threEntry);
                int cellid_now = (i)*1e6 + j*10000 + k;
                if(cellid_now!=_threcellID){
                	cout << i << ", " << j << ", " << k << "   " << _threcellID << " " << cellid_now;
                	cout << "    something is wrong at thredhold search" << endl;
                }

                if(i==0 && j==2 && cehck_canvas_mip){
                    cout << i << ", " << j << ", " << k << endl;
                    mipthrecanvas->cd();
                    mipthrecanvas->cd(k+1);
                    _thre_histo_tempCorrected->Draw();
                    // double thre_abso = channel_thre[i][j][k] + HG_ped[i][j][k];
                    // TLine *v_line= new TLine(thre_abso,0,thre_abso,100); //declare the vertical line 
                    // v_line->SetLineColor(kRed);
                    // v_line->SetLineWidth(2);
                    // v_line->Draw("same");
                }
                _thre_histo_tempCorrected->Rebin(5);
                int nbin = _thre_histo_tempCorrected->GetXaxis()->GetNbins();
                double histo_min = _thre_histo_tempCorrected->GetXaxis()->GetXmin();
                double histo_max = _thre_histo_tempCorrected->GetXaxis()->GetXmax();
                int binWidth = (histo_max - histo_min) / nbin;

                double maximum = _thre_histo_tempCorrected->GetBinContent(_thre_histo_tempCorrected->GetMaximumBin());
                int pedestalBin = _thre_histo_tempCorrected->FindBin(HG_ped[i][j][k]);
                int threBin = _thre_histo_tempCorrected->FindFirstBinAbove(maximum * constantFraction, 1, pedestalBin+nbin_above_pedestal, -1);
                double thre = binWidth * threBin + histo_min;

                channel_thre[i][j][k] = thre - HG_ped[i][j][k];
                if(i>=4 && i<28) threshold_check->Fill(thre - HG_ped[i][j][k]);
                threEntry++;
                //cout << i << ", " << j << ", " << k << "   gain : "  << _SPSgain << "  error : " << _SPSgainError << endl;
                // SPSgain_histo->Fill(_SPSgain);
                // gain_graph->SetPoint(count_gain_graph,count_gain_graph,_SPSgain);
                // gain_graph->SetPointError(count_gain_graph,0,_SPSgainError);
                // count_gain_graph++;
                mip_thre_histo[i][j][k] = (TH1F*)_thre_histo_tempCorrected->Clone();
                if(i==0 && j==2 && cehck_canvas_mip){
                    // cout << i << ", " << j << endl;
                    // mipthrecanvas->cd();
                    mipthrecanvas->cd(k+1);
                    // _thre_histo_tempCorrected->Draw();
                    double thre_abso = thre;
                    TLine *v_line= new TLine(thre_abso,0,thre_abso,100); //declare the vertical line 
                    v_line->SetLineColor(kRed);
                    v_line->SetLineWidth(2);
                    v_line->Draw("same");
                }
                f1->SetParameter(0,channel_thre[i][j][k]);
                f1->SetParameter(1,1);
                f1->SetParameter(2,maximum);
                f1->SetParameter(3,HG_ped[i][j][k]+MIP[i][j][k]);
                f1->SetParameter(4,1);
                _thre_histo_tempCorrected->Fit("f1","NQ","",HG_ped[i][j][k],HG_ped[i][j][k]+MIP[i][j][k]*1.5);
                if(i==0 && j==2 && k==5){
                    cout << i << ", " << j << ", " << k << endl;
                    mipthrecanvas->cd();
                    mipthrecanvas->cd(k+1);
                    _thre_histo_tempCorrected->Draw();
                    f1->Draw("same");
                }
            }
        }
    }
    threfile->Close();
    cout << "sps threshold reading finished" << endl;
  ///////////////////////////////////////////////////







}
