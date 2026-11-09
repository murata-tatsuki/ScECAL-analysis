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
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include "TLegend.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "EBUdecode.h"
#include "EBUdecode.cxx"
#include "langaus.C"        // langaus(chargeH, &fitFunc, &peakP, &peakPError);
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// エネルギー分布を求めるマクロ
// calibを使う

const bool pedestal_figure_saving = false;
const bool mip_figure_saving = true;
const bool gain_figure_saving = true;

const int layerNu = 32;
const int chipNu = 6;
const int channelNu = 36;
const int doublelayerNu = 2;
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

    const double _xInterval = 5.3;    // 300 um gap in width direction
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
    //     _position[0] = -y0;
    //     _position[1] = -x0;
    // }
    // if(LayerIDs%4==1) {
    //     _position[0] = x0;
    //     _position[1] = y0;
    // }
    // if(LayerIDs%4==2) {
    //     _position[0] = y0;
    //     _position[1] = x0;
    // }
    // if(LayerIDs%4==3) {
    //     _position[0] = -x0;
    //     _position[1] = -y0;
    // }
    // _position[2] = layerZ[LayerIDs];

    return _position;
}


int main(int argc, char* argv[])
{
    if(argc != 2){                                                     //エラー処理
        cout << "usage: ./between_files  output.root " << endl;
        return 1;
    }
  gROOT->SetBatch(kTRUE);

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;



  double CycleID, TriggerID;
  vector<int> *CellID = nullptr;
  vector<int> *BCID = nullptr;
  vector<int> *HitTag = nullptr;
  vector<int> *GainTag = nullptr;
  vector<double> *HG_Charge = nullptr;
  vector<double> *LG_Charge = nullptr;
  vector<vector<double>> *Temperature = nullptr;

  // double cycleID_ssa, triggerID_ssa;
  // vector<int> *cellIDs_ssa = nullptr;
  // vector<int> *BCIDs_ssa = nullptr;
  // vector<int> *hitTags_ssa = nullptr;
  // vector<int> *gainTags_ssa = nullptr;
  // vector<double> *charges_ssa = nullptr;
  // vector<double> *times_ssa = nullptr;
  // vector<double> *temp_ssa = nullptr;
  // vector<double> *posX_ssa = nullptr;
  // vector<double> *posY_ssa = nullptr;
  // vector<double> *posZ_ssa = nullptr;
  // vector<int> *ssaTag_ssa = nullptr;

  // double EnergyDep;
  // vector<int> *cellID = nullptr;
  // vector<double> *Hit_Energy = nullptr;
  // vector<double> *Hit_X = nullptr;
  // vector<double> *Hit_Y = nullptr;
  // vector<double> *Hit_Z = nullptr;

  double accurateBinY[10] = {-113.3,-68.3, -67.9,-22.9, -22.5,22.5, 22.9,67.9, 68.3,113.3};
  double accurateBinX[84] = {0};
  for(int i=0;i<rowNu;i++){
    accurateBinX[2*i] = -111.15 + i*_xInterval;
    accurateBinX[2*i+1] = accurateBinX[2*i] + 5.0;
  }
  const int range_maximum = 12000;
  const int range_maximum_PS = 500;
  const double binWidth_PS = 0.1;


  TFile fileout(argv[1],"RECREATE");
  // gain_histo->Write();


  // pedestalは昔のコードかも *要確認
  /*
  cout<<"=========== Pedestal Extract ============"<<endl;
    TFile* pedFile = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/pedestal2023_SPS.root");
    TTree* pedTree = (TTree*)pedFile->Get("ChnLevel"); 
    if(!pedTree)  {cout<<"!!! GET PED FILE FAILED !!!"<<endl;}
    vector<int>* _pedCellID = nullptr;
    vector<double>* _pedMeanHigh = nullptr;
    vector<double>* _pedSigmaHigh = nullptr;
    vector<double>* _pedMeanLow = nullptr;
    vector<double>* _pedSigmaLow = nullptr;
    pedTree ->SetBranchAddress("CellID",&_pedCellID);
    pedTree ->SetBranchAddress("PedHighMean",&_pedMeanHigh);
    pedTree ->SetBranchAddress("PedHighSig",&_pedSigmaHigh);
    pedTree ->SetBranchAddress("PedLowMean",&_pedMeanLow);
    pedTree ->SetBranchAddress("PedLowSig",&_pedSigmaLow);

    double meanPedHigh15  = 0 ; double meanPedHigh10  = 0;
    double sigmaPedHigh15 = 0 ; double sigmaPedHigh10 = 0;
    double chnNoPedHigh15 = 0 ; double chnNoPedHigh10 = 0;
    double meanPedLow15   = 0 ; double meanPedLow10   = 0;
    double sigmaPedLow15  = 0 ; double sigmaPedLow10  = 0;
    double chnNoPedLow15  = 0 ; double chnNoPedLow10  = 0;
    double pedMeanHigh[layerNu][chipNu][channelNu];
    double pedSigmaHigh[layerNu][chipNu][channelNu];
    double pedMeanLow[layerNu][chipNu][channelNu];
    double pedSigmaLow[layerNu][chipNu][channelNu];
    TH2F *ped_map[layerNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          pedMeanHigh[i_layer][i_chip][i_channel] = 1;
          pedSigmaHigh[i_layer][i_chip][i_channel] = 1;
          pedMeanLow[i_layer][i_chip][i_channel] = 1;
          pedSigmaLow[i_layer][i_chip][i_channel] = 1;
        }
      }
      double xNu_ = i_layer%2==0 ? columnNu : rowNu;
      double yNu_ = i_layer%2==0 ? rowNu : columnNu;
      double xMax_ = i_layer%2==0 ? _yMax : _xMax;
      double yMax_ = i_layer%2==0 ? _xMax : _yMax;
      // double accuratexNu_ = xNu_ * 2 - 1;
      // double accurateyNu_ = yNu_ * 2 - 1;
      // double* accurateXBin = i_layer%2==0 ? accurateBinY : accurateBinX;
      // double* accurateYBin = i_layer%2==0 ? accurateBinX : accurateBinY;
      ped_map[i_layer] = new TH2F(Form("ped_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
    }
    TH1F *pedestals[2][2];
    TH1F *pedestals_sigma[2][2];
    for(int pix=0;pix<2;pix++){
      for(int hl=0;hl<2;hl++){
        int pixel = pix == 0 ? 10 : 15;
        string HL = hl == 0 ? "high" : "low";
        pedestals[pix][hl] = new TH1F(Form("ped_%s_%d",HL.c_str(),pixel), Form("%s gain pedestal (%d um)",HL.c_str(),pixel), 300,300,600);
        pedestals_sigma[pix][hl] = new TH1F(Form("ped_sigma_%s_%d",HL.c_str(),pixel), Form("%s gain pedestal sigma (%d um)",HL.c_str(),pixel), 1000,0,10);
      }
    }
    for(int entry=0; entry!=pedTree->GetEntries(); ++entry){
      pedTree->GetEntry(entry);
      //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
      for(size_t ihit=0; ihit!=_pedCellID->size(); ++ihit)
      {
        int layerID = _pedCellID->at(ihit)/1e5;
        int chipID = _pedCellID->at(ihit)%100000/1e4;
        int channelID = _pedCellID->at(ihit)%100;        

        if (layerID>3 && layerID<28) {
          // if(_ChiSqr/(double)_NDF>2) continue;
          chnNoPedHigh10++;
          chnNoPedLow10++;
          meanPedHigh10 += _pedMeanHigh->at(ihit);
          meanPedLow10  += _pedMeanLow ->at(ihit);
          sigmaPedHigh10 += _pedSigmaHigh->at(ihit);
          sigmaPedLow10  += _pedSigmaLow ->at(ihit);
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
          chnNoPedHigh15++;
          chnNoPedLow15++;
          meanPedHigh15 += _pedMeanHigh->at(ihit);
          meanPedLow15  += _pedMeanLow ->at(ihit);
          sigmaPedHigh15 += _pedSigmaHigh->at(ihit);
          sigmaPedLow15  += _pedSigmaLow ->at(ihit);
        }
        pedMeanHigh[layerID][chipID][channelID] = _pedMeanHigh->at(ihit);
        pedSigmaHigh[layerID][chipID][channelID] = _pedSigmaHigh->at(ihit);
        pedMeanLow[layerID][chipID][channelID] = _pedMeanLow->at(ihit);
        pedSigmaLow[layerID][chipID][channelID] = _pedSigmaLow->at(ihit);
      }
    }
    meanPedHigh15 /= chnNoPedHigh15;
    meanPedLow15  /= chnNoPedLow15;
    meanPedHigh10 /= chnNoPedHigh10;
    meanPedLow10  /= chnNoPedLow10;
    sigmaPedHigh15 /= chnNoPedHigh15;
    sigmaPedLow15  /= chnNoPedLow15;
    sigmaPedHigh10 /= chnNoPedHigh10;
    sigmaPedLow10  /= chnNoPedLow10;
    cout << " HIGH GAIN PEDESTAL ... ..."  << endl;
    cout << " channel count(15um:10um) : " << chnNoPedHigh15 << " , " << chnNoPedHigh10 << endl;
    cout << " mean of pedestal High : "    << meanPedHigh15  << " , " << meanPedHigh10  << endl;
    cout << " sigma of pedestal High: "    << sigmaPedHigh15 << " , " << sigmaPedHigh10 << endl;
    cout << " LOW GAIN PEDESTAL ... ..."   << endl;
    cout << " channel count(15um:10um) : " << chnNoPedLow15  << " , " << chnNoPedLow10  << endl;
    cout << " mean of pedestal Low : "     << meanPedLow15   << " , " << meanPedLow10   << endl;
    cout << " sigma of pedestal Low: "     << sigmaPedLow15  << " , " << sigmaPedLow10  << endl;

    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if(pedMeanHigh[il][ip][ic]==1) cout << il << ", " << ip << ", " << ic << endl;
          if (il>3 && il<28) {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh10 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow10 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh10 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow10 : pedSigmaLow[il][ip][ic];

            pedestals[0][0]->Fill(pedMeanHigh[il][ip][ic]);
            pedestals[0][1]->Fill(pedMeanLow[il][ip][ic]);
            pedestals_sigma[0][0]->Fill(pedSigmaHigh[il][ip][ic]);
            pedestals_sigma[0][1]->Fill(pedSigmaLow[il][ip][ic]);
          } else {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh15 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow15 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh15 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow15 : pedSigmaLow[il][ip][ic];

            if(il<30){
              pedestals[1][0]->Fill(pedMeanHigh[il][ip][ic]);
              pedestals[1][1]->Fill(pedMeanLow[il][ip][ic]);
              pedestals_sigma[1][0]->Fill(pedSigmaHigh[il][ip][ic]);
              pedestals_sigma[1][1]->Fill(pedSigmaLow[il][ip][ic]);
            }
          }

          if(ip==5 && ic>29) continue;
          double* position = EBUdecode_itr(il,ip,ic);
          double int_x = il%2==0 ? _yInterval : _xInterval;
          double int_y = il%2==0 ? _xInterval : _yInterval;
          double xMax_ = il%2==0 ? _yMax : _xMax;
          double yMax_ = il%2==0 ? _xMax : _yMax;
          int ix_ = (position[0]+xMax_+0.5)/int_x;
          int iy_ = (position[1]+yMax_+0.5)/int_y;
          // int x_ = il%2==0 ? ix_ : iy_;
          // int y_ = il%2==0 ? iy_ : ix_;

          ped_map[il]->SetBinContent(ix_+1,iy_+1, pedMeanHigh[il][ip][ic]);
        }
      }
    }
    if(pedestal_figure_saving){
      TFile fileout_pedestal("../result/calibrations/pedestal.root","RECREATE");
      fileout_pedestal.cd();
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_ped_map_per_layer%d", i_layer), Form("ped map Layer%d", i_layer), 2560, 2560);
        ped_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/pedestal_map_Layer%d.png",i_layer));
        ped_map[i_layer]->Write(Form("pedestal_map_Layer%d", i_layer));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_ped_map_layer%d", i), Form("ped map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          ped_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/ped_map%d.png",i));
        c1->Write(Form("ped_map%d",i));
        delete c1;
      }
      for(int i=0;i<2;i++){
        int pixel = i == 0 ? 10 : 15;
        TCanvas* c1 = new TCanvas(Form("c_ped_%d", i), Form("ped %d", i), 2560, 1440);
        c1->Divide(2,2);
        c1->cd();
        for(int ii=0;ii<2;ii++){
          c1->cd(ii+1);
          gStyle->SetOptStat(1);
          pedestals[i][ii]->Draw();
          c1->cd(ii+3);
          pedestals_sigma[i][ii]->Draw();
        }
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/ped_%dum.png",pixel));
        c1->Write(Form("ped_%dum.png",pixel));
        delete c1;
      }
    }
  pedFile->Close();
  fileout.cd();
  */




  /*
  cout<<"=========== MIP Extract ============"<<endl;
    TFile *fileinMIP;
    fileinMIP = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/all_auto_muon_v4_trackfit.root");
    TTree *treeMIP = (TTree*) fileinMIP->Get("MIP_Fit");
    if(!treeMIP)  {cout<<"!!! GET MIP FILE FAILED !!!"<<endl;}
    int entry_max_MIP = treeMIP->GetEntries();
    Int_t _mipCellID, _NDF;
    Double_t _landauMPV, _LandauWidth, _GauSigma, _ChiSqr;
    treeMIP->SetBranchAddress("CellID",&_mipCellID);
    treeMIP->SetBranchAddress("LandauMPV",&_landauMPV);
    treeMIP->SetBranchAddress("LandauWidth",&_LandauWidth);
    treeMIP->SetBranchAddress("GauSigma",&_GauSigma);
    treeMIP->SetBranchAddress("ChiSquare",&_ChiSqr);
    treeMIP->SetBranchAddress("NDF",&_NDF);
    double MIPPeakValue[layerNu][chipNu][channelNu];
    double MIPWidth[layerNu][chipNu][channelNu];
    double MIPGauSigma[layerNu][chipNu][channelNu];
    TH2F *mip_map[layerNu];
    TH2F *mipWidth_map[layerNu];
    TH2F *mipGauSigma_map[layerNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          MIPPeakValue[i_layer][i_chip][i_channel] = 1;
          MIPWidth[i_layer][i_chip][i_channel] = 0;
          MIPGauSigma[i_layer][i_chip][i_channel] = 0;
        }
      }
      double xNu_ = i_layer%2==0 ? columnNu : rowNu;
      double yNu_ = i_layer%2==0 ? rowNu : columnNu;
      double xMax_ = i_layer%2==0 ? _yMax : _xMax;
      double yMax_ = i_layer%2==0 ? _xMax : _yMax;
      // double accuratexNu_ = xNu_ * 2 - 1;
      // double accurateyNu_ = yNu_ * 2 - 1;
      // double* accurateXBin = i_layer%2==0 ? accurateBinY : accurateBinX;
      // double* accurateYBin = i_layer%2==0 ? accurateBinX : accurateBinY;
      mip_map[i_layer] = new TH2F(Form("mip_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
      mipWidth_map[i_layer] = new TH2F(Form("mipWidth_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
      mipGauSigma_map[i_layer] = new TH2F(Form("mipGauSigma_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
    }
    TH1F *MIPs[2];
    TH1F *MIPs_width[2];
    TH1F *MIPs_GauSgima[2];
    for(int pix=0;pix<2;pix++){
      int pixel = pix == 0 ? 10 : 15;
      MIPs[pix] = new TH1F(Form("mip_%d",pixel), Form("MIP (%d um)",pixel), 1000,0,1000);
      MIPs_width[pix] = new TH1F(Form("mip_width_%d",pixel), Form("MIP landau width (%d um)",pixel), 1000,0.1,200);
      MIPs_GauSgima[pix] = new TH1F(Form("MIPs_GauSgima_%d",pixel), Form("MIP gaussian sigma (%d um)",pixel), 1000,0.1,200);
    }
    double meanMIP15  = 0 ; double meanMIP10  = 0;
    double chnNoMIP15 = 0 ; double chnNoMIP10 = 0;
    for(int entry=0; entry!=treeMIP->GetEntries(); entry++){
      treeMIP->GetEntry(entry);
      int layerID = _mipCellID/1e5;
      if(layerID>=30)  continue;  
      int chipID  = (_mipCellID%100000) /1e4;
      int chanID  = _mipCellID%100;
      if(_landauMPV<0) continue;
      if (layerID>3 && layerID<28) {
        // if(_ChiSqr/(double)_NDF>2) continue;
        chnNoMIP10++;
        meanMIP10 +=_landauMPV;
      } else {
        // if(_ChiSqr/(double)_NDF>1.7) continue;
        chnNoMIP15++;
        meanMIP15 +=_landauMPV;
      }
      MIPPeakValue[layerID][chipID][chanID] = _landauMPV;
      MIPWidth[layerID][chipID][chanID] = _LandauWidth;
      MIPGauSigma[layerID][chipID][chanID] = _GauSigma;
    }
    meanMIP10 /= chnNoMIP10;
    meanMIP15 /= chnNoMIP15;
    cout<<" channel count(15um:10um) : "<<chnNoMIP15<<" , "<<chnNoMIP10<<endl;
    cout<<" mean of MIPPeakValue : "<<meanMIP15<<" , "<<meanMIP10<<endl;

    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if(ip==5 && ic>29) continue;
          double means = (il>3 && il<28) ? meanMIP10 : meanMIP15;
          int pix = (il>3 && il<28) ? 0 : 1;

          double* position = EBUdecode_itr(il,ip,ic);
          double int_x = il%2==0 ? _yInterval : _xInterval;
          double int_y = il%2==0 ? _xInterval : _yInterval;
          double xMax_ = il%2==0 ? _yMax : _xMax;
          double yMax_ = il%2==0 ? _xMax : _yMax;
          int ix_ = (position[0]+xMax_+0.5)/int_x;
          int iy_ = (position[1]+yMax_+0.5)/int_y;
          // int x_ = il%2==0 ? ix_ : iy_;
          // int y_ = il%2==0 ? iy_ : ix_;
          if(MIPPeakValue[il][ip][ic]!=1) {
            MIPs[pix]->Fill(MIPPeakValue[il][ip][ic]);
            mip_map[il]->SetBinContent(ix_+1,iy_+1, MIPPeakValue[il][ip][ic]);
          }
          if(MIPWidth[il][ip][ic]>0.1) {
            MIPs_width[pix]->Fill(MIPWidth[il][ip][ic]);
            mipWidth_map[il]->SetBinContent(ix_+1,iy_+1, MIPWidth[il][ip][ic]);
          }
          if(MIPGauSigma[il][ip][ic]>0.1) {
            MIPs_GauSgima[pix]->Fill(MIPGauSigma[il][ip][ic]);
            mipGauSigma_map[il]->SetBinContent(ix_+1,iy_+1, MIPGauSigma[il][ip][ic]);
          }
          MIPPeakValue[il][ip][ic] = (MIPPeakValue[il][ip][ic]==1) ? means : MIPPeakValue[il][ip][ic];
        }
      }
    }
    if(mip_figure_saving){
      TFile fileout_mip("../result/calibrations/mip.root","RECREATE");
      fileout_mip.cd();
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_mip_map_per_layer%d", i_layer), Form("MIP map Layer%d", i_layer), 2560, 2560);
        mip_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/raw/mip_map_Layer%d.png",i_layer));
        mip_map[i_layer]->Write(Form("mip_map_Layer%d", i_layer));
        delete c1;
      }
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_mipWidth_map_per_layer%d", i_layer), Form("MIP width map Layer%d", i_layer), 2560, 2560);
        mipWidth_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/raw/mipWidth_map_Layer%d.png",i_layer));
        mipWidth_map[i_layer]->Write(Form("mipWidth_map_Layer%d", i_layer));
        delete c1;
      }
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_mipGauSig_map_per_layer%d", i_layer), Form("MIP gaussian sigma map Layer%d", i_layer), 2560, 2560);
        mipGauSigma_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/raw/mipGauSigma_map_Layer%d.png",i_layer));
        mipGauSigma_map[i_layer]->Write(Form("mipGauSigma_map_Layer%d", i_layer));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_mip_map_layer%d", i), Form("MIP map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          mip_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/mip_map%d.png",i));
        c1->Write(Form("mip_map%d",i));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_mipWidth_map_layer%d", i), Form("MIP width map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          mipWidth_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/mipWidth_map%d.png",i));
        c1->Write(Form("mipWidth_map%d",i));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_mipGauSigma_map_layer%d", i), Form("MIP gaussian sigma map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          mipGauSigma_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/mip/mipGauSigma_map%d.png",i));
        c1->Write(Form("mipGauSigma_map%d",i));
        delete c1;
      }
      TCanvas* c1_ = new TCanvas(Form("c_mip"), Form("MIP"), 2560, 1440);
      c1_->Divide(3,2);
      c1_->cd();
      for(int i=0;i<2;i++){
        int pixel = (i == 0) ? 10 : 15;
        int id = i*3;
        c1_->cd(id+1);
        gStyle->SetOptStat(1);
        MIPs[i]->Draw();
        c1_->cd(id+2);
        MIPs_width[i]->Draw();
        c1_->cd(id+3);
        MIPs_GauSgima[i]->Draw();
      }
      c1_->Update();
      c1_->SaveAs(Form("../result/calibrations/figures/mip/mip.png"));
      c1_->Write(Form("mip"));
      delete c1_;
    }
  fileinMIP->Close();
  fileout.cd();
  */


  /*
  cout<<"=========== Threshold Extract ============"<<endl;
    TFile *ThresholdFile = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root","READ");
    TTree* ThresholdTree = (TTree*)ThresholdFile->Get("threshold");
    Int_t _thresholdCellID;
    Double_t _threshold, _thresholdSigma, _thoresold_with_pedestal;
    ThresholdTree->SetBranchAddress("CellID",&_thresholdCellID);
    ThresholdTree->SetBranchAddress("threshold",&_threshold);
    ThresholdTree->SetBranchAddress("sigma",&_thresholdSigma);

    double thresholdValue[layerNu][chipNu][channelNu];
    double thresholdSigmaValue[layerNu][chipNu][channelNu];
    double thresholdValue_chip[layerNu][chipNu];
    double thresholdSigmaValue_chip[layerNu][chipNu];
    double thresholdValue_chipCount[layerNu][chipNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          thresholdValue[i_layer][i_chip][i_channel] = 1;
          thresholdSigmaValue[i_layer][i_chip][i_channel] = 0;
        }
        thresholdValue_chip[i_layer][i_chip] = 1;
        thresholdSigmaValue_chip[i_layer][i_chip] = 0;
        thresholdValue_chipCount[i_layer][i_chip] = 0;
      }
    }

    double meanThreshold15=0; double meanThreshold10=0;
    double sigmaThreshold15=0; double sigmaThreshold10=0;
    double chnNoThreshold15=0; double chnNoThreshold10=0;
    for(int entry=0; entry!=ThresholdTree->GetEntries(); entry++){
      ThresholdTree->GetEntry(entry);
      int layerID = _thresholdCellID/1e5;
      if(layerID>=30)  continue;  
      int chipID  = (_thresholdCellID%100000) /1e4;
      int chanID   = _thresholdCellID%100;
      if (layerID>3 && layerID<28) {
        // if(_ChiSqr/(double)_NDF>2) continue;
        chnNoThreshold10++;
        meanThreshold10  += _threshold;
        sigmaThreshold10 += _thresholdSigma;
        // meanThreshold10 += _thoresold_with_pedestal - thres;
        // meanThreshold10 += _thoresold_with_pedestal -_threshold + thres;
      } else {
        // if(_ChiSqr/(double)_NDF>1.7) continue;
        chnNoThreshold15++;
        meanThreshold15  += _threshold;
        sigmaThreshold15 += _thresholdSigma;
        // meanThreshold15 += _thoresold_with_pedestal - thres;
        // meanThreshold15 += _thoresold_with_pedestal -_threshold + thres;
      }
        // if (_threshold>50) {
        thresholdValue_chip[layerID][chipID] += _threshold;
        thresholdSigmaValue_chip[layerID][chipID] += _thresholdSigma;
        thresholdValue_chipCount[layerID][chipID]++;
        // }
        thresholdValue[layerID][chipID][chanID] = _threshold;
        thresholdSigmaValue[layerID][chipID][chanID] = _thresholdSigma;
        // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal - thres;
        // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal -_threshold + thres;
    }
    meanThreshold10 /= chnNoThreshold10;
    meanThreshold15 /= chnNoThreshold15;
    sigmaThreshold10 /= chnNoThreshold10;
    sigmaThreshold15 /= chnNoThreshold15;
    cout << " channel count(15um:10um) : " << chnNoThreshold15 << " , " << chnNoThreshold10 << endl;
    cout << " mean of Threshold : " << meanThreshold15 << " , " << meanThreshold10 << endl;
    cout << " sigma of Threshold : " << sigmaThreshold15 << " , " << sigmaThreshold10 << endl;

    for (int il=0; il<30; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if (il>3 && il<28) {
            thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? thresholdValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdValue[il][ip][ic];
            thresholdSigmaValue[il][ip][ic] = (thresholdSigmaValue[il][ip][ic]==0) ? thresholdSigmaValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdSigmaValue[il][ip][ic];
          } 
        }
      }
    }
  ThresholdFile->Close();
  fileout.cd();
  */



  cout<<"=========== SiPMGain Extract ============"<<endl;
    TFile* SiPMGainFile = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/20230504_1425_Calib_LED_all.root");
    TTree* SiPMGainTree = (TTree*)SiPMGainFile->Get("T_Event"); 
    if(!SiPMGainTree)  {cout<<"!!! GET SiPM Gain FILE FAILED !!!"<<endl;}
    Int_t _gainCellID;
    Double_t _gain, _gainError;
    SiPMGainTree->SetBranchAddress("cellIDs",&_gainCellID);
    SiPMGainTree->SetBranchAddress("gain",&_gain);
    SiPMGainTree->SetBranchAddress("gainError",&_gainError);

    double meanGain15=0; double meanGain10=0;
    double meanGainError15=0; double meanGainError10=0;
    double chnNoGain15=0; double chnNoGain10=0;
    double SiPMGain[layerNu][chipNu][channelNu];
    double SiPMGainError[layerNu][chipNu][channelNu];
    TH2F *gain_map[layerNu];
    TH2F *gainError_map[layerNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          SiPMGain[i_layer][i_chip][i_channel] = 0;
          SiPMGainError[i_layer][i_chip][i_channel] = 0;
        }
      }
      double xNu_ = i_layer%2==0 ? columnNu : rowNu;
      double yNu_ = i_layer%2==0 ? rowNu : columnNu;
      double xMax_ = i_layer%2==0 ? _yMax : _xMax;
      double yMax_ = i_layer%2==0 ? _xMax : _yMax;
      // double accuratexNu_ = xNu_ * 2 - 1;
      // double accurateyNu_ = yNu_ * 2 - 1;
      // double* accurateXBin = i_layer%2==0 ? accurateBinY : accurateBinX;
      // double* accurateYBin = i_layer%2==0 ? accurateBinX : accurateBinY;
      gain_map[i_layer] = new TH2F(Form("gain_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
      gainError_map[i_layer] = new TH2F(Form("gainError_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
    }
    TH1F *gains[2];
    TH1F *gains_error[2];
    for(int pix=0;pix<2;pix++){
      int pixel = pix == 0 ? 10 : 15;
      gains[pix] = new TH1F(Form("gain_%d",pixel), Form("gain (%d um)",pixel), 350,5,40);
      gains_error[pix] = new TH1F(Form("gain_error%d",pixel), Form("gain error (%d um)",pixel), 100,0,5);
    }
    for(int entry=0; entry!=SiPMGainTree->GetEntries(); ++entry){
      SiPMGainTree->GetEntry(entry);
      //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
      // for(size_t ihit=0; ihit!=_pedCellID->size(); ++ihit){
        // int layerID = _pedCellID->at(ihit)/1e5;
        // int chipID = _pedCellID->at(ihit)%100000/1e4;
        // int channelID = _pedCellID->at(ihit)%100;        
        int layerID = _gainCellID/1e6;
        if(layerID>=30)    continue;    
        int chipID  = (_gainCellID%1000000) /1e4;
        int chanID  = _gainCellID%100;

        if (layerID>3 && layerID<28) {
          // if(_ChiSqr/(double)_NDF>2) continue;
          chnNoGain10++;
          meanGain10 += _gain;
          meanGainError10 += _gainError;
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
          chnNoGain15++;
          meanGain15 += _gain;
          meanGainError15 += _gainError;
        }
        SiPMGain[layerID][chipID][chanID] = _gain;
        SiPMGainError[layerID][chipID][chanID] = _gainError;
      // }
    }
    meanGain10 /= chnNoGain10;
    meanGain15 /= chnNoGain15;
    meanGainError10 /= chnNoGain10;
    meanGainError15 /= chnNoGain15;
    cout << " channel count(15um:10um) : " << chnNoGain15 << " , " << chnNoGain10 << endl;
    cout << " mean of SiPM Gain : " << meanGain15 << " , " << meanGain10 << endl;
    cout << " mean of SiPM Gain Error: " << meanGainError15 << " , " << meanGainError10 << endl;
    
    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if(ip==5 && ic>29) continue;
          double means = (il>3 && il<28) ? meanGain10 : meanGain15;
          double errors = (il>3 && il<28) ? meanGainError10 : meanGainError15;
          int pix = (il>3 && il<28) ? 0 : 1;

          double* position = EBUdecode_itr(il,ip,ic);
          double int_x = il%2==0 ? _yInterval : _xInterval;
          double int_y = il%2==0 ? _xInterval : _yInterval;
          double xMax_ = il%2==0 ? _yMax : _xMax;
          double yMax_ = il%2==0 ? _xMax : _yMax;
          int ix_ = (position[0]+xMax_+0.5)/int_x;
          int iy_ = (position[1]+yMax_+0.5)/int_y;
          // int x_ = il%2==0 ? ix_ : iy_;
          // int y_ = il%2==0 ? iy_ : ix_;


          if(SiPMGain[il][ip][ic]>6.01 && SiPMGainError[il][ip][ic]!=0) {
            gains[pix]->Fill(SiPMGain[il][ip][ic]);
            gain_map[il]->SetBinContent(ix_+1,iy_+1, SiPMGain[il][ip][ic]);
            
            gains_error[pix]->Fill(SiPMGainError[il][ip][ic]);
            gainError_map[il]->SetBinContent(ix_+1,iy_+1, SiPMGainError[il][ip][ic]);
          }
          SiPMGain[il][ip][ic] = (SiPMGain[il][ip][ic]==0) ? means : SiPMGain[il][ip][ic];
          SiPMGainError[il][ip][ic] = (SiPMGainError[il][ip][ic]==0) ? errors : SiPMGainError[il][ip][ic];
        }
      }
    }
    if(gain_figure_saving){
      TFile fileout_gain("../result/calibrations/gain.root","RECREATE");
      fileout_gain.cd();
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_gain_map_per_layer%d", i_layer), Form("gain map Layer%d", i_layer), 2560, 2560);
        gain_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/gain/raw/gain_map_Layer%d.png",i_layer));
        gain_map[i_layer]->Write(Form("gain_map_Layer%d", i_layer));
        delete c1;
      }
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_gainError_map_per_layer%d", i_layer), Form("gain error map Layer%d", i_layer), 2560, 2560);
        gainError_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/gain/raw/gainError_map_Layer%d.png",i_layer));
        gainError_map[i_layer]->Write(Form("gainError_map_Layer%d", i_layer));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_gain_map_layer%d", i), Form("gain map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          gain_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/gain/gain_map%d.png",i));
        c1->Write(Form("gain_map%d",i));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_gainError_map_layer%d", i), Form("gain error map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          gainError_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/calibrations/figures/gain/gainError_map%d.png",i));
        c1->Write(Form("gainError_map%d",i));
        delete c1;
      }
      TCanvas* c1_ = new TCanvas(Form("c_gain"), Form("gain"), 2560, 1440);
      c1_->Divide(2,2);
      c1_->cd();
      for(int i=0;i<2;i++){
        int pixel = (i == 0) ? 10 : 15;
        int id = i*2;
        c1_->cd(id+1);
        gStyle->SetOptStat(1);
        gains[i]->Draw();
        c1_->cd(id+2);
        gains_error[i]->Draw();
      }
      c1_->Update();
      c1_->SaveAs(Form("../result/calibrations/figures/gain/gain.png"));
      c1_->Write(Form("gain"));
      delete c1_;
    }
  SiPMGainFile->Close();
  fileout.cd();




    
    


  






  /*
void Extract::HLGainRatioExtract(string inHLGainRatioFileName)
{
    cout<<"=========== HLGainRatioExtract ============"<<endl;
    fill_n(&HLRatio[0][0][0],         LAYERNO*CHIPNO*CHANNO, 1);
    fill_n(&HighADCSwitch[0][0][0], LAYERNO*CHIPNO*CHANNO, 1);
    TFile* HLRatioFile = new TFile(TString(inHLGainRatioFileName),"READ");
    TTree* HLRatioTree = (TTree*)HLRatioFile->Get("InterCalib");
    if(!HLRatioTree)    cout<<"!!! GET HLRatio FILE FAILED !!!"<<endl;
    Int_t _HLRatioCellID;
    Double_t _slope, _slopeError, _XMax;
    HLRatioTree->SetBranchAddress("CellID",&_HLRatioCellID);
    HLRatioTree->SetBranchAddress("Slope",&_slope);
    HLRatioTree->SetBranchAddress("SlopeError",&_slopeError);
    HLRatioTree->SetBranchAddress("XMax",&_XMax);
    
    double meanSlope       = 0; 
    double chnNoSlope      = 0;
    double minHighSwitch = 4000;
    for(int entry=0; entry!=HLRatioTree->GetEntries(); entry++)
    {
        HLRatioTree->GetEntry(entry);
        int layerID = _HLRatioCellID/1e5;
        if(layerID>=LAYERNO)    continue;    
        int chipID  = (_HLRatioCellID%100000) /1e4;
        int chanID   = _HLRatioCellID%100;
        
        meanSlope += 1. / _slope;
        chnNoSlope++;
        minHighSwitch = ((_XMax-600) < minHighSwitch) ? (_XMax-600) : minHighSwitch;

        HLRatio[layerID][chipID][chanID] = 1. / _slope;
        HighADCSwitch[layerID][chipID][chanID] = _XMax-600;
    }
    meanSlope /= chnNoSlope;
    cout << " channel count : "      << chnNoSlope << endl;
    cout << " mean of HLRatio : " << meanSlope  << endl;
    
    for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
                    HLRatio[il][ip][ic] = (HLRatio[il][ip][ic]==1) ? meanSlope : HLRatio[il][ip][ic];
                    HighADCSwitch[il][ip][ic] = (HighADCSwitch[il][ip][ic]==1) ? minHighSwitch : HighADCSwitch[il][ip][ic];
            }
        }
    }
    HLRatioFile->Close(); // jiaxuan
}
  */



}