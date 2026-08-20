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
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLine.h"
#include "TGaxis.h"
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include "TLegend.h"
#include "TCanvas.h"
#include "EBUdecode.h"
#include "EBUdecode.cxx"


double* EBUdecode_itr(int LayerIDs, int ChipIDs, int ChannelIDs){
  // static const int layerNu = 32;
  static const int chipNu = 6;
  static const int chnNu = 36;
  int decodeID[chipNu][chnNu] = { 0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13, 55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,
                                  24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,
                                  149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                                  150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
                                  191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181, 139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
                                  137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, 210, 210, 210, 210, 210, 210 };

  int ScintillatorIDs = decodeID[ChipIDs][ChannelIDs];
  double layerZ;

  const double _xInterval = 5.3;  // 300 um gap in width direction
  const double _yInterval = 45.4; // 400 um gap in length direction
  const int rowNu = 42;
  const int columnNu = 5;
  int _yID = ScintillatorIDs / rowNu;
  int _xID = ScintillatorIDs % rowNu;
  static double _position[3] = { 0 };
  double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
  double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;

  // for prototype test
  if (LayerIDs % 2 == 0) {
    _position[0] = -y0;
    _position[1] = -x0;
  }
  if (LayerIDs % 2 == 1) {
    _position[0] = -x0;
    _position[1] = -y0;
  }
  if (LayerIDs % 2 == 0)
    layerZ = 1 + LayerIDs / 2 * 19.9;
  //else layerZ = 12.95+(LayerIDs-1)/2*19.9;
  else
    layerZ = 12.2 + (LayerIDs - 1) / 2 * 19.9;
  _position[2] = layerZ;
  // for grouped CR test
  // if(LayerIDs%4==0) {
  //   _position[0] = -y0;
  //   _position[1] = -x0;
  // }
  // if(LayerIDs%4==1) {
  //   _position[0] = x0;
  //   _position[1] = y0;
  // }
  // if(LayerIDs%4==2) {
  //   _position[0] = y0;
  //   _position[1] = x0;
  // }
  // if(LayerIDs%4==3) {
  //   _position[0] = -x0;
  //   _position[1] = -y0;
  // }
  // _position[2] = layerZ[LayerIDs];

  return _position;
}

const string beamType = "sps";
// const string beamType = "ps";

int main(int argc, char* argv[])
{ 
  // if(argc < 3){                           //エラー処理
  //   cout << "usage: ./between_files  output.root Nenegry E1...En E1Nfiles...EnNfiles input1.root input2.root  ...  cog_range b_only_best figure_path " << endl;
  //   return 1;
  // }
  gROOT->SetBatch(kTRUE);

  //cout << argv[1] << endl;
  // cout << "=====>  " << argv[1] << endl;


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
  //               2,1,4,2,1,2,3,3,3,1,
  //               2,4,1,1,0,0,0,0,9,3};
  // hit chn<=10
  double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,2,1,4,2,1,2,3,3,3,1,2,4,1,2,0,1,0,0,9,3};

  TH1F* hEff0 = new TH1F("hEff0","",30,0,30);
  ifstream filelist(Form("./datalist_%s.txt",bramType.c_str()));
  // ifstream filelist("./datalist_CR.txt");
  // ifstream filelist("./datalist_pi_tilt.txt");
  if(!filelist.good()){cout<<"open datalist file failed !"<<endl; exit(0);}
  int fileNumber = 0;
  string runName = "all";
  for(int j=0; !filelist.eof(); j++){
    fileNumber++;
    string dname;
    filelist>>dname;
    if(filelist.eof()) break;
  }
  if(fileNumber==1){
    for(int j=0; !filelist.eof(); j++){
      string dname;
      filelist>>dname;
      runName = dname;
      if(filelist.eof()) break;
    }
    runName.erase(0,runName.find_last_of("/")+1);
    runName.erase(runName.find_last_of("."));
  }

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
  TH2F *gainMap[30];    // channelごとのgain
  TH2F *mipMap[30];     // channelごとのMIP
  TH2F *eff_2d[30];     // channelごとのdetection efficiency
  TH2F *hit_2d[30];     // channelごとのhit数
  TH2F *fit_2d[30];     // track fit のhit数 (trackが通過しているかどうか)
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
  double layer_z[30][2] = { { 0, 1},{ 1,  12.2},{ 2,  20.9},{ 3,  32.1},{ 4,  40.8},{ 5,  52},{ 6,  60.7},{ 7,  71.9},{ 8,  80.6},{ 9,  91.8},{10, 100.5},
                {11, 111.7},{12, 120.4},{13, 131.6},{14, 140.3},{15, 151.5},{16, 160.2},{17, 171.4},{18, 180.1},{19, 191.3},{20,   200},
                {21, 211.2},{22, 219.9},{23, 231.1},{24, 239.8},{25,   251},{26, 259.7},{27, 270.9},{28, 279.6},{29, 290.8} };


  double SPSgain[layerNu][chipNu][channelNu];
  double SPSgain_error[layerNu][chipNu][channelNu];
  double MIP[layerNu][chipNu][channelNu];
  double MIP_entry[layerNu][chipNu][channelNu];
  double channel_thre[layerNu][chipNu][channelNu];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        SPSgain[i][j][k] = 0;
        SPSgain_error[i][j][k] = 0;
        MIP[i][j][k] = 0;
        MIP_entry[i][j][k] = 0;
        channel_thre[i][j][k] = 0;
      }
    }
  }
  /*
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
            cout << "  something is wrong at mip search" << endl;
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
            cout << "  something is wrong at gain search" << endl;
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
    TFile *threfile = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root");
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
            cout << "  something is wrong at thredhold search" << endl;
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
          mip_thre_histo[i][j][k]->SetName(Form("mip_thre_histo_%d_%d_%d",i,j,k));
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
        }
      }
    }
    threfile->Close();
    cout << "sps threshold reading finished" << endl;
  ///////////////////////////////////////////////////

  */




  int deadCells[65][3] = { {0, 1, 3}, {0, 1, 33}, {0, 1, 34}, {1, 27, 0}, {2, 1, 1}, {2, 1, 28}, {2, 1, 30}, {2, 4, 37}, {3, 36, 4}, {4, 4, 11}, {5, 29, 2}, {5, 0, 3}, {5, 26, 3}, {6, 0, 7}, {6, 0, 8}, 
          {7, 4, 0}, {8, 2, 31}, {9, 29, 4}, {10, 1, 20}, {10, 2, 11}, {11, 3, 0}, {12, 1, 31}, {12, 1, 32}, {12, 2, 9}, {12, 2, 11}, {13, 8, 0}, {13, 9, 0}, {14, 0, 28}, {15, 25, 0}, 
          {15, 17, 1}, {16, 2, 40}, {16, 3, 15}, {16, 3, 19}, {17, 13, 3}, {17, 38, 3}, {17, 40, 3}, {18, 2, 40}, {18, 3, 10}, {18, 3, 12}, {19, 2, 4}, {20, 2, 26}, {20, 3, 0}, {21, 16, 1}, 
          {21, 4, 3}, {21, 1, 4}, {21, 40, 4}, {22, 0, 24}, {23, 3, 0}, {23, 4, 0}, {24, 3, 20}, {25, 30, 4}, {26, 4, 15}, {27, 35, 4}, {28, 0, 2}, {28, 0, 27}, {28, 0, 32}, {28, 2, 9}, 
          {28, 2, 26}, {28, 3, 3}, {28, 3, 5}, {28, 3, 21}, {28, 4, 26}, {29, 8, 0}, {29, 12, 2}, {29, 38, 2} };
          




  // ==============================================================
  //for(int j=0; j!=100; j++){
  for(int j=0; !filelist.eof(); j++){
    cout << j << endl;
    nbin++;
    string dname;
    filelist>>dname;
    if(filelist.eof()) break;
    std::vector<int>* _cell = nullptr;
    std::vector<double>* _posX = nullptr;
    std::vector<double>* _posY = nullptr;
    std::vector<double>* _posZ = nullptr;
    std::vector<int>* _newCell = nullptr;
    std::vector<double>* _charge = nullptr;
    std::vector<double>* _trackFitPars = nullptr;
    std::vector<double>* _residualX = nullptr;
    std::vector<double>* _residualY = nullptr;
    std::vector<double>* _newPosX = nullptr;
    std::vector<double>* _newPosY = nullptr;
    std::vector<double>* _newPosZ = nullptr;
    std::vector<double>* _temp = nullptr;  

    TFile* file = new TFile(dname.c_str());
    TTree* fNtuple1 = (TTree*)file->Get("T_Event");
    int _totalEntries = fNtuple1->GetEntries();
    if(!fNtuple1) cout<<"open RawtoRoot file failed "<<endl;
    fNtuple1 ->SetBranchAddress("CellID",&_cell);
    fNtuple1 ->SetBranchAddress("Hit_X",&_posX);
    fNtuple1 ->SetBranchAddress("Hit_Y",&_posY);
    fNtuple1 ->SetBranchAddress("Hit_Z",&_posZ);
    fNtuple1 ->SetBranchAddress("hitCellnew",&_newCell);
    fNtuple1 ->SetBranchAddress("energyDep",&_charge);
    fNtuple1 ->SetBranchAddress("trackFitPars",&_trackFitPars);
    fNtuple1 ->SetBranchAddress("residualX",&_residualX);
    fNtuple1 ->SetBranchAddress("residualY",&_residualY);
    fNtuple1 ->SetBranchAddress("hitPosXnew",&_newPosX);
    fNtuple1 ->SetBranchAddress("hitPosYnew",&_newPosY);
    fNtuple1 ->SetBranchAddress("hitPosZnew",&_newPosZ);
    fNtuple1 ->SetBranchAddress("Temperature",&_temp);
    // fNtuple1 ->SetBranchAddress("TemperatureNew",&_temp);

    TH1F* htemp[LayerNo];
    int totalHitLayer[LayerNo];
    for(int i=0; i<LayerNo; i++){
      char hname[50];
      sprintf(hname,"layer_%d",i);
      // htemp[i] = new TH1F(hname,hname,100,12,30);
      totalHitLayer[i] = 0;
    }

    cout<<j<<" : "<<dname<<" : "<<_totalEntries<<endl;
    for(int entry=0; entry!=_totalEntries; ++entry)
    {
      if(entry%10000==0) 
        cout<<" Event : "<<entry<<endl;
      fNtuple1->GetEntry(entry);

      vector<tuple<int, double, double>> channelCoordinates(_cell->size());
      for(int i=0;i<_cell->size();i++){
        int _layer = _cell->at(i)/1e5;
        channelCoordinates.at(i) = make_tuple(_cell->at(i), _posX->at(i), _posY->at(i));
      }
      // sort(channelCoordinates.begin(), channelCoordinates.end());
      vector< int > hitLayer;
      double slopeX = _trackFitPars->at(0);
      double slopeY = _trackFitPars->at(4);
      double interceptX = _trackFitPars->at(2);
      double interceptY = _trackFitPars->at(6);
      for(int hit=0; hit!=(int)_newCell->size(); ++hit)
      {
        if(_newCell->at(hit)==-1) continue;
        int cellIDnow = _newCell->at(hit);
        int _layerID = _newCell->at(hit)/1e5;
        if(_layerID>29) continue;

        vector<int>::iterator result = find(hitLayer.begin(),hitLayer.end(),_layerID);
        if(result==hitLayer.end()) {
          totalHitLayer[_layerID]++;
          hitLayer.push_back(_layerID);
          // htemp[_layerID]->Fill(_temp->at(hit));
        }

        // double x_now = _newPosX->at(hit);
        // double y_now = _newPosY->at(hit);
        double x_now=-200, y_now=-200;
        for(int ii=0;ii<_cell->size();ii++){
          const int& cellIDcheck = get<0>(channelCoordinates.at(ii));
          if(cellIDnow==cellIDcheck){
            x_now = get<1>(channelCoordinates.at(ii));
            y_now = get<2>(channelCoordinates.at(ii));
            break;
          }
        }
        hit_2d[_layerID]->Fill(x_now,y_now);
      }
      for(int ilayer=0; ilayer<30; ++ilayer){
        double z_now = layer_z[ilayer][1];
        double x_fit = slopeX * z_now + interceptX;
        double y_fit = slopeY * z_now + interceptY;
        fit_2d[ilayer]->Fill(x_fit,y_fit);

        for(int hit=0; hit!=(int)_newCell->size(); ++hit)
        {
          if(_newCell->at(hit)==-1) continue;
          int _layerID = _newCell->at(hit)/1e5;
          int _chipID = _newCell->at(hit)%100000/1e4;
          int _channelID = _newCell->at(hit)%100;
          if(_layerID>29 || _layerID!=ilayer) continue;

          // double x_now = _newPosX->at(hit);   // 今はpositionが修正されているから、元の座標を使うように変更   trackのsrcを変更する
          // double y_now = _newPosY->at(hit);
          double* position = EBUdecode_itr(_layerID,_chipID,_channelID);
          double x_now = position[0];
          double y_now = position[1];
          
          double x_int = ilayer%2==0 ? 45.0/2 : 5.0/2;
          double y_int = ilayer%2==0 ? 5.0/2 : 45.0/2;
          bool condition = (x_now-x_int<x_fit && x_now+x_int>x_fit) && (y_now-y_int<y_fit && y_now+y_int>y_fit);
          if(condition){
            fit_hit_2d[ilayer]->Fill(x_fit,y_fit);
            if(true){
              int nbinX = fit_hit_2d[ilayer]->GetXaxis()->FindBin(x_fit);
              int nbinY = fit_hit_2d[ilayer]->GetYaxis()->FindBin(y_fit);
              if(nbinX%2==0 && nbinY%2==0) cout << "   Event : " << entry << "  " <<  _layerID << ", " << x_fit << ", " << y_fit << ", " << nbinX << ", " << nbinY << endl;
            }

            bool deadCondition = false;
            for(int idead=0;idead<65;idead++){
              if(deadCells[idead][0]!=_layerID) continue;
              double x_dead = ilayer%2==0 ? (deadCells[idead][1]-2.0)*_yInterval : deadCells[idead][1]*_xInterval - 108.5;
              double y_dead = ilayer%2==0 ? deadCells[idead][2]*_xInterval - 108.5 : (deadCells[idead][2]-2.0)*_yInterval;
              deadCondition = deadCondition || (x_now-x_int<x_dead && x_now+x_int>x_dead) && (y_now-y_int<y_dead && y_now+y_int>y_dead);
            }
            if(deadCondition) cout << "   Event : " << entry << " has dead cell hit  " << _newCell->at(hit) << endl;
          }
        }
      }
    }
    totalEntry += _totalEntries;
    for(int i=0; i!=LayerNo; ++i)
    {
      NhitLayer[i] += totalHitLayer[i];
      double eff = (double)totalHitLayer[i]/(double)_totalEntries;
      // cout << i << "  " << totalHitLayer[i] << ", " << _totalEntries << ", " << eff << endl;
      if(i!=1&&i!=29) eff += blindChn[i]*0.5/100.;
      gEff[i]->SetPoint(j,j,eff);
      // if(j==105) hEff0->SetBinContent(i+1,eff);
      if(j==fileNumber) {
        eff = NhitLayer[i] / totalEntry;
        hEff0->SetBinContent(i+1,eff);
        if(i!=1&&i!=29) effmean+=eff;
      }
    }
    // for(int i=0; i<LayerNo; i++){
    //   //double temp = htemp[i]->GetMean();
    //   double temp = htemp[i]->GetMean()/20.;
    //   gTemp[i]->SetPoint(j,j,temp);
    // }
  }

  if(false){  // 時間変動を表示させるかどうか
    // ==============================================================
    TCanvas* c = new TCanvas("c","c",10,10,800,600);
    // TLine* line = new TLine(0,0.935,30,0.935);
    TLine* line = new TLine(0,effmean/28,30,effmean/28);
    line->SetLineWidth(2);
    line->SetLineColor(2);
    hEff0->GetYaxis()->SetTitle("Efficiency");
    hEff0->GetXaxis()->SetTitle("LayerID");
    hEff0->GetXaxis()->CenterTitle();
    hEff0->GetYaxis()->CenterTitle();
    hEff0->GetYaxis()->SetRangeUser(0,1.1);
    hEff0->SetLineColor(4);
    hEff0->SetLineWidth(2);
    hEff0->SetFillStyle(3004);
    hEff0->SetFillColor(4);
    hEff0->Draw("HIST");
    line->Draw("same");
    TLatex* latex = new TLatex(.4,.7,"Preliminary");
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.06);
    latex->SetTextAlign(33);
    latex->SetTextColor(15);
    latex->Draw();
    
    // ==============================================================
    TCanvas* cc = new TCanvas("cc","cc",10,10,1800,800);
    gStyle->SetOptStat(false);
    TH2D* hTop = new TH2D("hTop","",10,0,nbin,10,0.75,1.01);
    hTop->GetYaxis()->SetTitle("Efficiency");
    hTop->GetXaxis()->SetTitle("Date");
    hTop->GetXaxis()->CenterTitle();
    hTop->GetYaxis()->CenterTitle();
    hTop->GetXaxis()->SetLabelSize(0);
    hTop->GetYaxis()->SetTitleSize(0.06);
    hTop->GetYaxis()->SetTitleOffset(.41);
    cc->cd();
    TPad* pad1 = new TPad("pad1","pad1",.01,.55,1,1);
    pad1->SetMargin(.07,.1,.03,.1);
    pad1->Draw();
    pad1->cd();
    hTop->Draw();
    TLegend* lg = new TLegend(.9,.03,.98,.9);
    char lname[50];
    for(int i=0; i<LayerNo; i++)
    {
      if(i<=3 || i >=28){
        gEff[i]->SetMarkerStyle(21);
        gEff[i]->SetMarkerColor(1+i);
        gEff[i]->SetLineColor(1+i);
        gEff[i]->SetMarkerSize(1);
        gEff[i]->SetLineWidth(1);
        gEff[i]->SetFillColor(false);
        gEff[i]->Draw("sameLP");
        sprintf(lname,"Layer_%d",i);
        lg->AddEntry(gEff[i],lname);
      }
    }
    lg->Draw("same");
    // -----------------------------------------------------------
    TH2D* hBottom = new TH2D("hBottom","",9,0,nbin,10,0.75,1.01);
    hBottom->GetYaxis()->SetTitle("Efficiency");
    hBottom->GetXaxis()->SetTitle("Date");
    hBottom->GetXaxis()->CenterTitle();
    hBottom->GetYaxis()->CenterTitle();
    hBottom->GetXaxis()->SetLabelSize(0.06);
    hBottom->GetXaxis()->SetTitleSize(0.06);
    hBottom->GetYaxis()->SetTitleSize(0.06);
    hBottom->GetXaxis()->SetTitleOffset(1.2);
    hBottom->GetYaxis()->SetTitleOffset(.41);
    for(int bin=0; bin<9; bin++){
      char hname[50];
      sprintf(hname,"%d",(bin+1)*10);
      hBottom->GetXaxis()->SetBinLabel(bin+1,hname);
    }
    cc->cd();
    TPad* pad2 = new TPad("pad2","pad2",.01,.05,1,.55);
    pad2->SetMargin(.07,.1,.15,.05);
    pad2->Draw();
    pad2->cd();
    hBottom->Draw();
    TLegend* lg1 = new TLegend(.9,.15,.98,.95);
    for(int i=0; i<LayerNo; i++)
    {
      if(i>3 &&i <28){
        gEff[i]->SetMarkerStyle(20);
        gEff[i]->SetMarkerColor(1+i);
        gEff[i]->SetLineColor(1+i);
        gEff[i]->SetMarkerSize(1);
        gEff[i]->SetLineWidth(1);
        gEff[i]->SetFillColor(false);
        gEff[i]->Draw("sameLP");
        sprintf(lname,"Layer_%d",i);
        lg1->AddEntry(gEff[i],lname);
      }
    }
    lg1->Draw("same");

    // -----------------------------------------------------------
    TCanvas* ctemp = new TCanvas("ctemp","ctemp",10,10,1800,600);
    int _MaxTemp=30;
    int _MinTemp=12;
    //TH2D* h = new TH2D("h","",10,0,nbin,10,_MinTemp,_MaxTemp);
    TH2D* h = new TH2D("h","",10,0,nbin,10,0.6,1.5);
    h->Draw();
    h->GetYaxis()->SetTitle("Temperature [ {}^{o}C ]");
    h->GetYaxis()->SetTitleOffset(.55);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitle("Date");
    h->GetXaxis()->SetTitleOffset(1.15);
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);

    TLegend* lg2 = new TLegend(.9,.1,.97,.9);
    for(int i=0; i<LayerNo; i++){
      gTemp[i]->SetMarkerColor(1+i);
      gTemp[i]->SetMarkerSize(1.5);
      gTemp[i]->SetMarkerStyle(20);
      gTemp[i]->SetFillColor(0);
      gTemp[i]->Draw("sameP");
      char lname[50];
      sprintf(lname,"Layer_%d",i);
      lg2->AddEntry(gTemp[i],lname);
    }
    lg2->SetFillColor(0);
    lg2->SetLineWidth(0);
    lg2->Draw();

    // -----------------------------------------------------------
    TCanvas* cet = new TCanvas("cet","cet",10,10,1800,600);
    gPad->SetGridy();
    const int layer = 14;
    const int layer2 = 15;
    TH2D* hTempEff = new TH2D("hTempEff","",9,0,nbin,10,0.55,1.3);
    hTempEff->GetYaxis()->SetTitle("Efficiency");
    hTempEff->GetXaxis()->SetTitle("Date");
    hTempEff->GetXaxis()->CenterTitle();
    hTempEff->GetYaxis()->CenterTitle();
    hTempEff->GetXaxis()->SetLabelSize(0.04);
    hTempEff->GetXaxis()->SetTitleSize(0.04);
    hTempEff->GetYaxis()->SetLabelSize(0.04);
    hTempEff->GetYaxis()->SetTitleSize(0.04);
    hTempEff->GetXaxis()->SetTitleOffset(1.15);
    hTempEff->GetYaxis()->SetTitleOffset(.6);
    for(int bin=0; bin<9; bin++){
      char hname[50];
      sprintf(hname,"%d",(bin+1)*10);
      hTempEff->GetXaxis()->SetBinLabel(bin+1,hname);
    }
    hTempEff->Draw();

    //TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
    //              gPad->GetUxmax(),gPad->GetUymax(),_MinTemp,_MaxTemp,510,"+L");
    TGaxis* axis = new TGaxis(nbin,.75,
                  nbin,1.3,_MinTemp,_MaxTemp,510,"+L");
    axis->SetLineColor(2);
    axis->SetLabelColor(2);
    axis->SetTitle("Temperature [ {}^{o}C ]");
    axis->SetTitleColor(2);
    axis->SetTitleOffset(.6);
    axis->SetLabelSize(.04);
    axis->CenterTitle();
    axis->Draw();

    char hname[50];
    sprintf(hname,"Layer_%d",layer);
    hTempEff->SetTitle(hname);
    TLegend* lg3 = new TLegend(.6,.75,.85,.9);
    gEff[layer]->SetMarkerColor(3);
    gEff[layer]->SetMarkerSize(1.5);
    gEff[layer]->SetMarkerStyle(20);
    gEff[layer]->SetLineWidth(2);
    gEff[layer]->SetLineColor(3);
    gEff[layer]->Draw("sameLP");
    gEff[layer2]->SetMarkerColor(4);
    gEff[layer2]->SetMarkerSize(1.5);
    gEff[layer2]->SetMarkerStyle(21);
    gEff[layer2]->SetLineWidth(2);
    gEff[layer2]->SetLineColor(4);
    gEff[layer2]->Draw("sameLP");
    gTemp[layer]->SetMarkerColor(2);
    gTemp[layer]->SetMarkerSize(1.5);
    gTemp[layer]->SetMarkerStyle(22);
    gTemp[layer]->SetLineWidth(2);
    gTemp[layer]->SetLineColor(2);
    gTemp[layer]->Draw("sameLP");
    gTemp[layer2]->SetMarkerColor(5);
    gTemp[layer2]->SetMarkerSize(1.5);
    gTemp[layer2]->SetMarkerStyle(23);
    gTemp[layer2]->SetLineWidth(2);
    gTemp[layer2]->SetLineColor(5);
    gTemp[layer2]->Draw("sameLP");
    sprintf(hname,"Layer_%d Detection Efficiency",layer);
    lg3->AddEntry(gEff[layer],hname);
    sprintf(hname,"Layer_%d Detection Efficiency",layer2);
    lg3->AddEntry(gEff[layer2],hname);
    sprintf(hname,"Layer_%d Average Temperature",layer);
    lg3->AddEntry(gTemp[layer],hname);
    sprintf(hname,"Layer_%d Average Temperature",layer2);
    lg3->AddEntry(gTemp[layer2],hname);
    lg3->Draw();
  }


  double channel_efficiency[30][columnNu][rowNu];
  double channel_efficiency_entry[30][columnNu][rowNu];
  double channel_fit_hit[30][columnNu][rowNu];
  double channel_hit[30][columnNu][rowNu];
  for(int i=0;i<30;i++){
    double denominator=0, numerator=0;
    for(int _x=1;_x<columnNu+1;_x++){
      for(int _y=1;_y<rowNu+1;_y++){
        int x_ = i%2==0 ? _x : _y;
        int y_ = i%2==0 ? _y : _x;
        // double efficiency = fit_hit_2d[i]->GetBinContent(x_,y_) / fit_2d[i]->GetBinContent(x_,y_);
        // if(fit_hit_2d[i]->GetBinContent(x_,y_)!=0) eff_2d[i]->SetBinContent(x_,y_, efficiency);
        double efficiency = fit_hit_2d[i]->GetBinContent(2*x_-1,2*y_-1) / fit_2d[i]->GetBinContent(2*x_-1,2*y_-1);
        denominator += fit_2d[i]->GetBinContent(2*x_-1,2*y_-1);
        numerator += fit_hit_2d[i]->GetBinContent(2*x_-1,2*y_-1);
        if(fit_hit_2d[i]->GetBinContent(2*x_-1,2*y_-1)!=0) eff_2d[i]->SetBinContent(x_,y_, efficiency);
        channel_efficiency[i][_x-1][_y-1] = efficiency;
        channel_efficiency_entry[i][_x-1][_y-1] = fit_2d[i]->GetBinContent(2*x_-1,2*y_-1);
        channel_fit_hit[i][_x-1][_y-1] = fit_hit_2d[i]->GetBinContent(2*x_-1,2*y_-1);
        channel_hit[i][_x-1][_y-1] = hit_2d[i]->GetBinContent(x_,y_);
        if(efficiency==0) cout << " ";
        if(efficiency<0.05) cout << i << ", " << x_-1 << ", " << y_-1 << endl;
      }
    }
    cout << " layer " << i << " efficiency " << numerator / denominator << endl;
  }

  bool check_relations = false;
  if(check_relations){
    TCanvas *check = new TCanvas("check","check",1);
    TCanvas *check_mip = new TCanvas("check_mip","check_mip",1);
    TCanvas *check_mip_15 = new TCanvas("check_mip_15","check_mip_15",1);
    TCanvas *check_mip_10 = new TCanvas("check_mip_10","check_mip_10",1);
    TCanvas *check_gain = new TCanvas("check_gain","check_gain",1);
    TCanvas *check_eff_vs_fit = new TCanvas("check_eff_vs_fit","check_eff_vs_fit",1);
    TCanvas *check_eff_vs_hit = new TCanvas("check_eff_vs_hit","check_eff_vs_hit",1);
    TCanvas *check_eff_vs_hit_fit = new TCanvas("check_eff_vs_hit_fit","check_eff_vs_hit_fit",1);
    TCanvas *check_eff_vs_thre = new TCanvas("check_eff_vs_thre","check_eff_vs_thre",1);
    TCanvas *check_eff_vs_thre_15 = new TCanvas("check_eff_vs_thre_15","check_eff_vs_thre_15",1);
    TCanvas *check_eff_vs_thre_10 = new TCanvas("check_eff_vs_thre_10","check_eff_vs_thre_10",1);
    // TCanvas *mipthrecanvas = new TCanvas("mipthrecanvas","mipthrecanvas",1);
    // mipthrecanvas->Divide(6, 6);
    // if(i==0 && j==2  cehck_canvas_mip){
    //         cout << i << ", " << j << endl;
    //         mipthrecanvas->cd();
    //         mipthrecanvas->cd(k+1);
    //         mip_thre_histo[i][j][k]->Draw();
    //         double thre_abso = channel_thre[i][j][k] + HG_ped[i][j][k];
    //         TLine *v_line= new TLine(thre_abso,0,thre_abso,100); //declare the vertical line 
    //         v_line->SetLineColor(kRed);
    //         v_line->SetLineWidth(2);
    //         v_line->Draw("same");
    //       }
    TGraph *relation = new TGraph();
    relation->GetXaxis()->SetTitle("efficiency");
    relation->GetYaxis()->SetTitle("gain");
    int relationCount = 0;
    TH2F *relation_mip = new TH2F("relation_mip","relation_mip;efficiency;mip", 100,0,1, 100,0,900);
    TH2F *relation_mip_15 = new TH2F("relation_mip_15","15 um SiPM;efficiency;mip", 100,0,1, 100,0,900);
    TH2F *relation_mip_10 = new TH2F("relation_mip_10","10 um SiPM;efficiency;mip", 100,0,1, 100,0,900);
    TH2F *relation_gain = new TH2F("relation_gain","relation_gain;efficiency;gain", 100,0,1, 80,5,35);
    TH2F *relation_eff_fit = new TH2F("relation_eff_fit","relation_eff_fit;efficiency;# of fit", 100,0,1, 100,0,4000);
    TH2F *relation_eff_fit_hit = new TH2F("relation_eff_fit_hit","relation_eff_fit_hit;efficiency;# of fit_hit", 100,0,1, 100,0,3500);
    TH2F *relation_eff_hit = new TH2F("relation_eff_hit","relation_eff_hit;efficiency;# of hit", 100,0,1, 100,0,3500);
    // TH2F *relation_eff_thre = new TH2F("relation_eff_thre","relation_eff_thre;efficiency;MPV-threshold [ADC]", 100,0,1, 160,0,800);
    // TH2F *relation_eff_thre_15 = new TH2F("relation_eff_thre_15","15 um SiPM;efficiency;MPV-threshold [ADC]", 100,0,1, 160,0,800);
    // TH2F *relation_eff_thre_10 = new TH2F("relation_eff_thre_10","10 um SiPM;efficiency;MPV-threshold [ADC]", 100,0,1, 100,0,100);
    // TH2F *relation_eff_thre = new TH2F("relation_eff_thre","relation_eff_thre;efficiency;threshold - pedestal [ADC]", 100,0,1, 200,0,200);
    // TH2F *relation_eff_thre_15 = new TH2F("relation_eff_thre_15","15 um SiPM;efficiency;threshold - pedestal [ADC]", 100,0,1, 200,0,200);
    // TH2F *relation_eff_thre_10 = new TH2F("relation_eff_thre_10","10 um SiPM;efficiency;threshold - pedestal [ADC]", 100,0,1, 200,0,200);
    TH2F *relation_eff_thre = new TH2F("relation_eff_thre","relation_eff_thre;efficiency;threshold / MIP", 100,0,1, 200,0,1.1);
    TH2F *relation_eff_thre_15 = new TH2F("relation_eff_thre_15","15 um SiPM;efficiency;threshold / MIP", 100,0,1, 200,0,1.1);
    TH2F *relation_eff_thre_10 = new TH2F("relation_eff_thre_10","10 um SiPM;efficiency;threshold / MIP", 100,0,1, 200,0,1.1);
    for(int i=0; i<30;i++){
      for(int j=0; j<chipNu;j++){
        for(int k=0; k<channelNu;k++){
          if(j==5&&k>29) continue;

          // if(i==0 && j==2 && cehck_canvas_mip){
          //   cout << i << ", " << j << ", " << k << endl;
          //   mipthrecanvas->cd();
          //   mipthrecanvas->cd(k+1);
          //   mip_thre_histo[i][j][k]->Draw();
          //   double thre_abso = channel_thre[i][j][k] + HG_ped[i][j][k];
          //   TLine *v_line= new TLine(thre_abso,0,thre_abso,100); //declare the vertical line 
          //   v_line->SetLineColor(kRed);
          //   v_line->SetLineWidth(2);
          //   v_line->Draw("same");
          //   mipthrecanvas->Update();
          // }

          double* position = EBUdecode_itr(i,j,k);
          double int_x = i%2==0 ? _yInterval : _xInterval;
          double int_y = i%2==0 ? _xInterval : _yInterval;
          double xMax_ = i%2==0 ? _yMax : _xMax;
          double yMax_ = i%2==0 ? _xMax : _yMax;
          int ix_ = (position[0]+xMax_+0.5)/int_x;
          int iy_ = (position[1]+yMax_+0.5)/int_y;
          int x_ = i%2==0 ? ix_ : iy_;
          int y_ = i%2==0 ? iy_ : ix_;

          gainMap[i]->SetBinContent(ix_+1,iy_+1, SPSgain[i][j][k]);
          mipMap[i]->SetBinContent(ix_+1,iy_+1, MIP[i][j][k]);

          double efficiency = channel_efficiency[i][x_][y_];
          if(channel_efficiency_entry[i][x_][y_]==0) continue;
          relation_gain->Fill(efficiency,SPSgain[i][j][k]);
          if(MIP[i][j][k]>0){
            relation_mip->Fill(efficiency,MIP[i][j][k]);
            if(i<=3 || i>=28) relation_mip_15->Fill(efficiency,MIP[i][j][k]);
            else relation_mip_10->Fill(efficiency,MIP[i][j][k]);
          }
          // if(efficiency>1 || efficiency<0 || SPSgain[i][j][k]<0) continue;
          // relation->SetPoint(relationCount++,efficiency,SPSgain[i][j][k]);
          // if(efficiency>1 || efficiency<0 || MIP[i][j][k]<0) continue;
          relation->SetPoint(relationCount++,efficiency,channel_hit[i][x_][y_]);
          relation_eff_fit->Fill(efficiency,channel_efficiency_entry[i][x_][y_]);
          relation_eff_fit_hit->Fill(efficiency,channel_fit_hit[i][x_][y_]);
          relation_eff_hit->Fill(efficiency,channel_hit[i][x_][y_]);

          // relation_eff_thre->Fill(efficiency, MIP[i][j][k] - channel_thre[i][j][k]);
          // if(efficiency<0.05) continue;
          // if(i<=3 || i>=28) relation_eff_thre_15->Fill(efficiency, MIP[i][j][k] - channel_thre[i][j][k]);
          // else if(MIP_entry[i][j][j]>500) relation_eff_thre_10->Fill(efficiency, MIP[i][j][k] - channel_thre[i][j][k]);
          //////////////////////////////////////////////////////////////////////////////////////
          // relation_eff_thre->Fill(efficiency, channel_thre[i][j][k]);
          // if(efficiency<0.05) continue;
          // if(i<=3 || i>=28) relation_eff_thre_15->Fill(efficiency, channel_thre[i][j][k]);
          // else if(MIP_entry[i][j][j]>500) relation_eff_thre_10->Fill(efficiency, channel_thre[i][j][k]);
          relation_eff_thre->Fill(efficiency, channel_thre[i][j][k]/MIP[i][j][k]);
          if(efficiency<0.05) continue;
          if(i<=3 || i>=28) relation_eff_thre_15->Fill(efficiency, channel_thre[i][j][k]/MIP[i][j][k]);
          else if(MIP_entry[i][j][j]>500) relation_eff_thre_10->Fill(efficiency, channel_thre[i][j][k]/MIP[i][j][k]);
        }
      }
    }
    // mipthrecanvas->Update();
    check->cd();
    // relation->SetMarkerStyle(4);
    // relation->SetMarkerSize(0.5);
    // relation->Draw("AP");
    // threshold_check->Draw();
    // MIP_entry_check->Draw();

    check_eff_vs_fit->cd();
    relation_eff_fit->Draw("colz");

    check_eff_vs_hit_fit->cd();
    relation_eff_fit_hit->Draw("colz");

    check_eff_vs_hit->cd();
    relation_eff_hit->Draw("colz");

    gStyle->SetOptStat(0);
    check_gain->cd();
    // gPad->SetGrid(0,0);
    relation_gain->Draw("colz");
    gStyle->SetOptStat(0);
    check_mip->cd();
    // gPad->SetGrid(0,0);
    relation_mip->Draw("colz");
    check_mip_15->cd();
    relation_mip_15->Draw("colz");
    check_mip_10->cd();
    relation_mip_10->Draw("colz");

    check_eff_vs_thre->cd();
    relation_eff_thre->Draw("colz");
    check_eff_vs_thre_15->cd();
    relation_eff_thre_15->Draw("colz");
    check_eff_vs_thre_10->cd();
    relation_eff_thre_10->Draw("colz");
  }


  // TCanvas *hitCanvas[30];
  // TCanvas *effCanvas[30];
  // for(int i=0; i<LayerNo; i++){
  //   hitCanvas[i] = new TCanvas(Form("hitCanvas_%d",i),Form("hitCanvas_%d",i),1200,600);
  //   hitCanvas[i]->Divide(2, 1);
  //   hitCanvas[i]->cd(1);
  //   hit_2d[i]->Draw("colz");
  //   hitCanvas[i]->cd(2);
  //   fit_2d[i]->Draw("colz");
  // 
  //   effCanvas[i] = new TCanvas(Form("effCanvas_%d",i),Form("effCanvas_%d",i),1200,600);
  //   effCanvas[i]->Divide(2, 1);
  //   effCanvas[i]->cd(1);
  //   fit_hit_2d[i]->Draw("colz");
  //   effCanvas[i]->cd(2);
  //   eff_2d[i]->GetZaxis()->SetRangeUser(0, 1);
  //   eff_2d[i]->Draw("colz");
  // }
  
  bool draw2dmap = true;
  bool draw2dmap_gain_mip = false;
  bool writeCanvas = true;
  if(draw2dmap){
    TCanvas *hitCanvas[2];
    TCanvas *fitCanvas[2];
    TCanvas *fithitCanvas[2];
    TCanvas *effCanvas[2];
    for(int i=0;i<2;i++){
      string eo = i%2==0 ? "even" : "odd";
      hitCanvas[i] = new TCanvas(Form("hitCanvas_%s_%d",eo.c_str(),i),Form("hitCanvas_%s_%d",eo.c_str(),i), 1440,810);
      hitCanvas[i]->Divide(4, 4);

      fitCanvas[i] = new TCanvas(Form("fitCanvas_%s_%d",eo.c_str(),i),Form("fitCanvas_%s_%d",eo.c_str(),i), 1440,810);
      fitCanvas[i]->Divide(4, 4);

      fithitCanvas[i] = new TCanvas(Form("fithitCanvas_%s_%d",eo.c_str(),i),Form("fithitCanvas_%s_%d",eo.c_str(),i), 1440,810);
      fithitCanvas[i]->Divide(4, 4);
  
      effCanvas[i] = new TCanvas(Form("effCanvas_%s_%d",eo.c_str(),i),Form("effCanvas_%s_%d",eo.c_str(),i), 1440,810);
      effCanvas[i]->Divide(4, 4);
    }
    for(int i=0; i<LayerNo; i++){
      int even_odd = i%2;
      int ilayer_eo = (i-i%2)/2;
  
      hitCanvas[even_odd]->cd();
      hitCanvas[even_odd]->cd(ilayer_eo+1);
      gStyle->SetOptStat(0);
      gPad->SetGrid(0,0);
      hit_2d[i]->Draw("colz");

      fitCanvas[even_odd]->cd();
      fitCanvas[even_odd]->cd(ilayer_eo+1);
      gStyle->SetOptStat(0);
      gPad->SetGrid(0,0);
      fit_2d[i]->Draw("colz");

      fithitCanvas[even_odd]->cd();
      fithitCanvas[even_odd]->cd(ilayer_eo+1);
      gStyle->SetOptStat(0);
      gPad->SetGrid(0,0);
      fit_hit_2d[i]->Draw("colz");
  
      effCanvas[even_odd]->cd();
      effCanvas[even_odd]->cd(ilayer_eo+1);
      gStyle->SetOptStat(0);
      eff_2d[i]->GetZaxis()->SetRangeUser(0, 1);
      gPad->SetGrid(0,0);
      eff_2d[i]->Draw("colz");
    }
    if(writeCanvas){
      hitCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/hitCanvas_even.png"),beamType.c_str(),runName.c_str());
      hitCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/hitCanvas_odd.png"),beamType.c_str(),runName.c_str());
      fitCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/fitCanvas_even.png"),beamType.c_str(),runName.c_str());
      fitCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/fitCanvas_odd.png"),beamType.c_str(),runName.c_str());
      fithitCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/fithitCanvas_even.png"),beamType.c_str(),runName.c_str());
      fithitCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/fithitCanvas_odd.png"),beamType.c_str(),runName.c_str());
      effCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/effCanvas_even.png"),beamType.c_str(),runName.c_str());
      effCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/effCanvas_odd.png"),beamType.c_str(),runName.c_str());
    }

    /*
    if(draw2dmap_gain_mip){
      TCanvas *mipCanvas[2];
      TCanvas *gainCanvas[2];
      for(int i=0;i<2;i++){
        string eo = i%2==0 ? "even" : "odd";
        mipCanvas[i] = new TCanvas(Form("mipCanvas_%s_%d",eo.c_str(),i),Form("mipCanvas_%s_%d",eo.c_str(),i), 1440,810);
        mipCanvas[i]->Divide(4, 4);

        gainCanvas[i] = new TCanvas(Form("gainCanvas_%s_%d",eo.c_str(),i),Form("gainCanvas_%s_%d",eo.c_str(),i), 1440,810);
        gainCanvas[i]->Divide(4, 4);
      }
      for(int i=0; i<LayerNo; i++){
        int even_odd = i%2;
        int ilayer_eo = (i-i%2)/2;

        mipCanvas[even_odd]->cd();
        mipCanvas[even_odd]->cd(ilayer_eo+1);
        mipMap[i]->Draw("colz");

        gainCanvas[even_odd]->cd();
        gainCanvas[even_odd]->cd(ilayer_eo+1);
        // gainMap[i]->GetZaxis()->SetRangeUser(0, 1);
        gainMap[i]->Draw("colz");
      }
      if(writeCanvas){
        // check_mip->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/efficiency_mip.png"),beamType.c_str());
        // check_gain->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/efficiency_gain.png)"));
        gainCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/gainCanvas_even.png"),beamType.c_str());
        gainCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/gainCanvas_odd.png"),beamType.c_str());
        mipCanvas[0]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/mipCanvas_even.png"),beamType.c_str());
        mipCanvas[1]->SaveAs(Form("../../result/trackFit/efficiency/figures/%s/%s/mipCanvas_odd.png"),beamType.c_str());
      }
    }
    */
  }
}