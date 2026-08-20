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


// track fitを用いて角channelのdetection efficitencyを求めるコード
// fitのparameterの分布の表示も行っている


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

int main(int argc, char* argv[])
{ 
  if(argc < 3){                           //エラー処理
    cout << "usage: ./o  output.root input_1.root ... input_n.root  figure_path " << endl;
    return 1;
  }
  gROOT->SetBatch(kTRUE);

  //cout << argv[1] << endl;
  // cout << "=====>  " << argv[1] << endl;

  cout << "=====>  " << argv[1] << endl;
	int rawfilenum = argc - 3;

  TFile fileout(argv[1],"RECREATE");




  const int LayerNo=30;
  // hit chn=0
  // double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,
  //               2,1,4,2,1,2,3,3,3,1,
  //               2,4,1,1,0,0,0,0,9,3};
  // hit chn<=10
  double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,2,1,4,2,1,2,3,3,3,1,2,4,1,2,0,1,0,0,9,3};

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

  TH1D* h_intercept[2];
  TH1D* h_intercept_error[2];
  TH1D* h_slope[2];
  TH1D* h_slope_error[2];
  for(int i=0;i<2;i++){
    string xy = i==0 ? "x" : "y";
    h_intercept[i] = new TH1D(Form("h_intercept_%s",xy.c_str()), Form(";%sz intercept [mm];Entries",xy.c_str()), 600, -150, 150);
    h_intercept_error[i] = new TH1D(Form("h_intercept_error_%s",xy.c_str()), Form(";%sz intercept error [mm];Entries",xy.c_str()), 200, 0, 2);
    h_slope[i]     = new TH1D(Form("h_slope_%s",xy.c_str()),     Form(";%sz slope;Entries",xy.c_str()),     500, -0.5, 0.5);
    h_slope_error[i]     = new TH1D(Form("h_slope_error_%s",xy.c_str()),     Form(";%sz slope error;Entries",xy.c_str()),     1000, 0, 0.01);
  }



  int deadCells[65][3] = { {0, 1, 3}, {0, 1, 33}, {0, 1, 34}, {1, 27, 0}, {2, 1, 1}, {2, 1, 28}, {2, 1, 30}, {2, 4, 37}, {3, 36, 4}, {4, 4, 11}, {5, 29, 2}, {5, 0, 3}, {5, 26, 3}, {6, 0, 7}, {6, 0, 8}, 
          {7, 4, 0}, {8, 2, 31}, {9, 29, 4}, {10, 1, 20}, {10, 2, 11}, {11, 3, 0}, {12, 1, 31}, {12, 1, 32}, {12, 2, 9}, {12, 2, 11}, {13, 8, 0}, {13, 9, 0}, {14, 0, 28}, {15, 25, 0}, 
          {15, 17, 1}, {16, 2, 40}, {16, 3, 15}, {16, 3, 19}, {17, 13, 3}, {17, 38, 3}, {17, 40, 3}, {18, 2, 40}, {18, 3, 10}, {18, 3, 12}, {19, 2, 4}, {20, 2, 26}, {20, 3, 0}, {21, 16, 1}, 
          {21, 4, 3}, {21, 1, 4}, {21, 40, 4}, {22, 0, 24}, {23, 3, 0}, {23, 4, 0}, {24, 3, 20}, {25, 30, 4}, {26, 4, 15}, {27, 35, 4}, {28, 0, 2}, {28, 0, 27}, {28, 0, 32}, {28, 2, 9}, 
          {28, 2, 26}, {28, 3, 3}, {28, 3, 5}, {28, 3, 21}, {28, 4, 26}, {29, 8, 0}, {29, 12, 2}, {29, 38, 2} };
          

  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    cout << irawfile << "/" << rawfilenum << endl;

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
    
    TFile *file = new TFile(argv[irawfile+2]);
    TTree *fNtuple1 = (TTree*)file->Get("T_Event");
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

    int totalHitLayer[LayerNo];
    for(int i=0; i<LayerNo; i++){
      char hname[50];
      sprintf(hname,"layer_%d",i);
      totalHitLayer[i] = 0;
    }

    cout<<irawfile<<" : "<<argv[irawfile+2]<<" : "<<_totalEntries<<endl;
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
      for(int ixy=0;ixy<2;ixy++){
        h_intercept[ixy]->Fill(_trackFitPars->at(ixy*4+2));
        h_intercept_error[ixy]->Fill(_trackFitPars->at(ixy*4+2+1));
        h_slope[ixy]->Fill(_trackFitPars->at(ixy*4));
        h_slope_error[ixy]->Fill(_trackFitPars->at(ixy*4+1));
      }
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

    delete file;
  }


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
        if(efficiency==0) cout << " ";
        if(efficiency<0.05) cout << i << ", " << x_-1 << ", " << y_-1 << endl;
      }
    }
    cout << " layer " << i << " efficiency " << (denominator==0 ? 0 : numerator / denominator) << endl;
  }

  
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
      hitCanvas[i] = new TCanvas(Form("hitCanvas_%s_%d",eo.c_str(),i),Form("hitCanvas_%s_%d",eo.c_str(),i), 2560, 1440);
      hitCanvas[i]->Divide(4, 4);

      fitCanvas[i] = new TCanvas(Form("fitCanvas_%s_%d",eo.c_str(),i),Form("fitCanvas_%s_%d",eo.c_str(),i), 2560, 1440);
      fitCanvas[i]->Divide(4, 4);

      fithitCanvas[i] = new TCanvas(Form("fithitCanvas_%s_%d",eo.c_str(),i),Form("fithitCanvas_%s_%d",eo.c_str(),i), 2560, 1440);
      fithitCanvas[i]->Divide(4, 4);
  
      effCanvas[i] = new TCanvas(Form("effCanvas_%s_%d",eo.c_str(),i),Form("effCanvas_%s_%d",eo.c_str(),i), 2560, 1440);
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
    TCanvas *fitParameterCanvas = new TCanvas(Form("fitParameterCanvas"),Form("fit parameters"), 2560, 1440);
    fitParameterCanvas->Divide(4, 2);
    fitParameterCanvas->cd();
    for(int ixy=0;ixy<2;ixy++){
      fitParameterCanvas->cd(ixy*4 + 1);
      gPad->SetLogy(1);
      gPad->SetGrid(0,0);
      h_slope[ixy]->Draw();

      fitParameterCanvas->cd(ixy*4 + 2);
      gPad->SetLogy(1);
      gPad->SetGrid(0,0);
      h_slope_error[ixy]->Draw();

      fitParameterCanvas->cd(ixy*4 + 3);
      gPad->SetLogy(1);
      gPad->SetGrid(0,0);
      h_intercept[ixy]->Draw();

      fitParameterCanvas->cd(ixy*4 + 4);
      gPad->SetLogy(1);
      gPad->SetGrid(0,0);
      h_intercept_error[ixy]->Draw();
    }
    if(writeCanvas){
      hitCanvas[0]->SaveAs(Form("%s/hitCanvas_even.png",argv[argc-1]));
      hitCanvas[1]->SaveAs(Form("%s/hitCanvas_odd.png",argv[argc-1]));
      fitCanvas[0]->SaveAs(Form("%s/fitCanvas_even.png",argv[argc-1]));
      fitCanvas[1]->SaveAs(Form("%s/fitCanvas_odd.png",argv[argc-1]));
      fithitCanvas[0]->SaveAs(Form("%s/fithitCanvas_even.png",argv[argc-1]));
      fithitCanvas[1]->SaveAs(Form("%s/fithitCanvas_odd.png",argv[argc-1]));
      effCanvas[0]->SaveAs(Form("%s/effCanvas_even.png",argv[argc-1]));
      effCanvas[1]->SaveAs(Form("%s/effCanvas_odd.png",argv[argc-1]));
      fitParameterCanvas->SaveAs(Form("%s/trackFit_parameters.png",argv[argc-1]));
    }
  }
}