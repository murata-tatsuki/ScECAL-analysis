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
#include <TDatime.h>
#include <map>
#include <cstdio>
#include "TLegend.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "EBUdecode.h"
#include "EBUdecode.cxx"
#include "langaus.C"
// #include "/home/jap/analyseCode/Diagnose/src/EBUdecode.cxx"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// LEDデータをまとめるコード
// 1runのLEDのデータをまとめるコード


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

const double timeWindow = 600.0; // 10分 = 600秒
const int total_cells = layerNu * chipNu * channelNu; // 6480 + 210*2



double ctap_calc(double p0, double p1){
  if(p0==0 || p1==0) return -10;

  double lambda = -log(p0);
  double true_p1 = lambda * p0;
  return 1 - p1/true_p1;
}

double ctap_eqf(double p0, double gain, double mean){
  double lambda = -log(p0);
  double mu = mean / gain;
  return mu / lambda;
}

double calc_intercept1(double a, double b, double c, int plmi=1){
  double result = c + plmi * b * sqrt(2*log(abs(a)));
  return result;
}

struct BinStat {
    long long count = 0;
    double sum = 0.0;
    double sum2 = 0.0;

    inline void add(double val) {
        count++;
        sum += val;
        sum2 += val * val;
    }

    // 平均と標準誤差 (SEM) を算出
    bool get_results(double &mean, double &sem) const {
        if (count == 0) return false;
        mean = sum / count;
        if (count > 1) {
            double variance = (sum2 - (sum * sum) / count) / (count - 1);
            sem = std::sqrt(std::max(0.0, variance) / count);
        } else {
            sem = 0.0;
        }
        return true;
    }
};

int main(int argc, char* argv[])
{ 
    if(argc < 5){                                                     //エラー処理
        cout << "usage: ./between_files output.root outputDirectory input1.root input2.root  ..." << endl;
        return 1;
    }
  struct FileInfo {
    int runNo;
    string start_str;
    int iarg;
  };

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = argc - 3;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  vector<FileInfo> filename;
  int total_entry_max=0;

	string picDirName(argv[2]);
	cout << picDirName << endl;

  double rawFileDate[rawfilenum];

  for(int i=0; i<rawfilenum; i++){
    string txtname(argv[i+3]);
    txtname.erase(0,txtname.find_last_of("Run")+1);
    int runNo = stoi(txtname);
    txtname.erase(0,txtname.find_first_of("_")+1);
    txtname.erase(txtname.find_last_of("."),-1);
    filename.push_back({runNo, txtname, i+3});
    // cout << argv[i+2] << ",  " << stoi(txtname) << endl;
  }
  std::sort(filename.begin(), filename.end(), [](const FileInfo& a, const FileInfo& b) { return a.runNo < b.runNo; });
  // sort(filename.begin(),filename.end());
  for(int i=0; i<rawfilenum; i++){
    filein[i] = new TFile(argv[filename[i].iarg]);
    tree[i] = (TTree*) filein[i]->Get("Raw_Hit");
    if(tree[i]!=nullptr) entry_max[i] = tree[i]->GetEntries();
    else entry_max[i] = 0;
  }

    // /megraid01/users/tsuji/scecal/led/data/analyseCode/LEDCalibration/LEDgroupID.txt にある
  int group_channel[14][15] = {{2 ,      6 ,      8 ,      1  ,     5  ,     202,     190,     192,     194,     179,     175,     187,     176,     173,     171},
    {0 ,      4 ,      7 ,      3  ,     201,     203,     189,     191,     193,     177,     188,     186,     178,     174,     172},
    {10,      12,      14,      204,     206,     208,     196,     198,     200,     185,     183,     181,     169,     167,     165},
    {9 ,      11,      13,      205,     207,     209,     195,     197,     199,     184,     182,     180,     170,     168,     166},
    {17,      21,      25,      16 ,     20 ,     24 ,     76 ,     78 ,     80 ,     164,     160,     156,     161,     157,     153},
    {15,      19,      23,      18 ,     22 ,     26 ,     75 ,     77 ,     79 ,     162,     158,     154,     163,     159,     155},
    {29,      33,      35,      28 ,     32 ,     88 ,     82 ,     84 ,     86 ,     152,     148,     73 ,     149,     146,     144},
    {27,      31,      34,      30 ,     87 ,     89 ,     81 ,     83 ,     85 ,     150,     74 ,     72 ,     151,     147,     145},
    {38,      42,      46,      37 ,     41 ,     45 ,     91 ,     93 ,     95 ,     108,     112,     116,     111,     115,     119},
    {36,      40,      44,      39 ,     43 ,     47 ,     90 ,     92 ,     94 ,     110,     114,     118,     109,     113,     117},
    {50,      54,      58,      49 ,     53 ,     57 ,     97 ,     99 ,     101,     120,     124,     128,     123,     127,     131},
    {48,      52,      56,      51 ,     55 ,     59 ,     96 ,     98 ,     100,     122,     126,     130,     121,     125,     129},
    {62,      66,      70,      61 ,     65 ,     69 ,     103,     105,     107,     132,     136,     140,     135,     139,     143},
    {60,      64,      68,      63 ,     67 ,     71 ,     102,     104,     106,     134,     138,     142,     133,     137,     141}};





  int cycleID, triggerID, Event_Time;
  vector<int> *CellID = nullptr;
  vector<int> *BCID = nullptr;
  vector<int> *HitTag = nullptr;
  vector<int> *GainTag = nullptr;
  vector<double> *HG_Charge = nullptr;
  vector<double> *LG_Charge = nullptr;
  vector<vector<double>> *temp = nullptr;


  // const int fit_sigma_maximum = 9;
  const int fit_sigma_maximum = 6;



  TFile fileout(argv[1],"RECREATE");
  // gain_histo->Write();



  cout<<"=========== PedestalExtract ============"<<endl;
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
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          pedMeanHigh[i_layer][i_chip][i_channel] = 1;
          pedSigmaHigh[i_layer][i_chip][i_channel] = 1;
          pedMeanLow[i_layer][i_chip][i_channel] = 1;
          pedSigmaLow[i_layer][i_chip][i_channel] = 1;
        }
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
          if (il>3 && il<28) {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh10 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow10 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh10 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow10 : pedSigmaLow[il][ip][ic];
          } else {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh15 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow15 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh15 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow15 : pedSigmaLow[il][ip][ic];
          }
        }
      }
    }
  pedFile->Close();
  fileout.cd();


  cout<<"=========== Threshold Extract ============"<<endl;
    TFile *ThresholdFile = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root","READ");
    // TTree* ThresholdTree = (TTree*)ThresholdFile->Get("threshold");
    TTree* ThresholdTree = (TTree*)ThresholdFile->Get("expErfThre");
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



  // dataを入れるもの
  TH1F *HG_histo[layerNu][chipNu][channelNu];
  TH1F *LG_histo[layerNu][chipNu][channelNu];
  TH1F *channel_temperature[layerNu][chipNu][channelNu];
  int temp_count[layerNu][chipNu][channelNu];
  TH1F *sensor_temperature[layerNu][16];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        HG_histo[i][j][k] = new TH1F(Form("HG_histo_%d_%d_%d",i,j,k),Form("high gain layer %d chip %d channel %d",i,j,k),4096,0,4095);
        LG_histo[i][j][k] = new TH1F(Form("LG_histo_%d_%d_%d",i,j,k),Form("low gain layer %d chip %d channel %d",i,j,k),4096,0,4095);
        channel_temperature[i][j][k] = new TH1F(Form("channel_temperature_%d_%d_%d",i,j,k),Form("temperature layer %d chip %d channel %d",i,j,k),200,15,35);
      }
    }
    for(int j=0; j<16;j++){
      sensor_temperature[i][j] = new TH1F(Form("sensor_temperature_%d_%d",i,j),Form("temperature layer%d sensor%d",i,j),1000,20,30);
      sensor_temperature[i][j]->SetXTitle("temperature [#circC]");
    }
  }




  // 保存形式
  double _gain, _gainError, _gausSigma, _gausSigmaError, _chiSquare, _NDF, _temp, _tempError, _ctap, _eqf;
  int _cellID_fit, _fit_entry;
  TF1 *_fitting;
  TTree* gain_tree = new TTree("T_Event", "gain calibration");
  gain_tree->Branch("CellID", &_cellID_fit);
  gain_tree->Branch("gain", &_gain);
  gain_tree->Branch("gainError", &_gainError);
  gain_tree->Branch("gausSigma", &_gausSigma);
  gain_tree->Branch("gausSigmaError", &_gausSigmaError);
  gain_tree->Branch("chiSquare", &_chiSquare);
  gain_tree->Branch("NDF", &_NDF);
  gain_tree->Branch("temp", &_temp);
  gain_tree->Branch("tempError", &_tempError);
  gain_tree->Branch("fitting", &_fitting);
  gain_tree->Branch("ctap", &_ctap);
  gain_tree->Branch("eqf", &_eqf);
  gain_tree->Branch("fit_entry", &_fit_entry);

  double channel_temperature_mean[32][6][36] = {0};
  double channel_temperature_sigma[32][6][36] = {0};




  cout << "made save TTree" << endl;

  vector<vector<BinStat>> cell_data(total_cells);
  double global_start_time = 0.0;

  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;
    if(entry_max[irawfile]==0) continue;

    // ★ ディスクI/O高速化: 一旦全ブランチの読み込みを無効化し、必要なものだけ有効化
    // tree[irawfile]->SetBranchStatus("*", 0);
    // tree[irawfile]->SetBranchStatus("Event_Time", 1);
    // tree[irawfile]->SetBranchStatus("CellID", 1);
    // tree[irawfile]->SetBranchStatus("HG_Charge", 1);

    tree[irawfile]->SetBranchAddress("Event_Time", &Event_Time);
    tree[irawfile]->SetBranchAddress("CellID", &CellID);
    tree[irawfile]->SetBranchAddress("HG_Charge", &HG_Charge);
    tree[irawfile]->SetBranchAddress("LG_Charge", &LG_Charge);
    tree[irawfile]->SetBranchAddress("HitTag", &HitTag);
    tree[irawfile]->SetBranchAddress("Temperature", &temp);

    string start_str = filename[irawfile].start_str;
    int year, month, day, hour, min, sec;
    sscanf(start_str.c_str(), "%4d%2d%2d_%2d%2d%2d", &year, &month, &day, &hour, &min, &sec);
    TDatime da(year, month, day, hour, min, sec);
    double current_run_start = (double)da.Convert();

    if (irawfile == 0) global_start_time = current_run_start;
    double run_offset_sec = current_run_start - global_start_time;

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);

      double total_elapsed_sec = run_offset_sec + Event_Time;
      size_t bin_idx = static_cast<size_t>(total_elapsed_sec / timeWindow);

      for(int i=0; i<(int)CellID->size(); i++){
        // cout << "   " << i << "/" << cellIDs->size() << endl;
        if(HitTag->at(i)==0) continue;
        int layerID = CellID->at(i)/1e5;
        int chipID = CellID->at(i)%100000/1e4;
        int channelID = CellID->at(i)%100;

        // int index = chipID * 36 + channelID;
        // if(find(group_channel[groupNum], group_channel[groupNum] + 15, index) == group_channel[groupNum] + 15) continue;

        double _hg_charge = HG_Charge->at(i);
        double _lg_charge = LG_Charge->at(i);
        // int _hitTags = hitTags->at(i);
        // double _x = posX->at(i);
        // double _y = posY->at(i);
        // double _z = posZ->at(i);

        HG_histo[layerID][chipID][channelID]->Fill(_hg_charge);
        LG_histo[layerID][chipID][channelID]->Fill(_lg_charge);

        double *_position = EBUdecode(layerID,chipID,channelID);
        // cout << _position[0] << ", " << _position[1] << ", " << _position[2] << endl;
        double SiPMtemp = 20;
        if(temp->at(layerID).size()!=0) {
          SiPMtemp = tempReconstruction(layerID, _position, temp->at(layerID));
          // cout << SiPMtemp << endl;
        }
        else SiPMtemp=20;
        channel_temperature[layerID][chipID][channelID]->Fill(SiPMtemp);
        

        int c_idx = layerID * 216 + chipID * 36 + channelID;
        if (c_idx < 0 || c_idx >= total_cells) continue;
        if (bin_idx >= cell_data[c_idx].size()) {
          cell_data[c_idx].resize(bin_idx + 1);
        }
        // その場で加算処理（メモリに電荷配列を保存しない）
        cell_data[c_idx][bin_idx].add(SiPMtemp);
      }

      for(int i=0;i<(int)temp->size();i++){
        for(int j=0;j<(int)temp->at(i).size();j++){
          sensor_temperature[i][j]->Fill(temp->at(i).at(j));
        }
      }
    }
    CellID = nullptr;
    BCID = nullptr;
    HitTag = nullptr;
    GainTag = nullptr;
    HG_Charge = nullptr;
    LG_Charge = nullptr;
    temp = nullptr;
  }

  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        channel_temperature_mean[i][j][k] = channel_temperature[i][j][k]->GetMean();
        channel_temperature_sigma[i][j][k] = channel_temperature[i][j][k]->GetRMS();
      }
    }
  }

  TH1F *gain = new TH1F("gain","gain;gain;",200/25*29,6,35);
  TH1F *gain_15 = new TH1F("gain_15","15um gain;gain;",200/25*29,6,35);
  TH1F *gain_10 = new TH1F("gain_10","10um gain;gain;",200/25*29,6,35);
  TH2F *gain_layer = new TH2F("gain_layer","gain vs layer;gain;layer",200/25*29,6,35,32,0,32);
  TH2F *sigma_layer = new TH2F("sigma_layer","sigma_layer",100,0,fit_sigma_maximum,32,0,32);
  TH1F *gain_15_ratio = new TH1F("gain_15_ratio","15um gain ratio;ratio (beam test gain / cosmic ray gain);",200,0,3);
  TH1F *gain_10_ratio = new TH1F("gain_10_ratio","10um gain ratio;ratio (beam test gain / cosmic ray gain);",200,0,3);

  TF1 *highfun = new TF1("highfun","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1]-[4])/[5])^2) + [6]*exp(-0.5*((x-[1]-[4]*2)/[7])^2)",300,900);
  highfun->SetNpx(4000);
  // TF1 *lowfun = new TF1("lowfun","gaus(0)+gaus(3)+gaus(6)",range_l,range_h);
  highfun->SetParLimits(4,6,34);
  highfun->SetParLimits(2,0,fit_sigma_maximum);
  highfun->SetParLimits(5,0,fit_sigma_maximum);
  highfun->SetParLimits(7,0,fit_sigma_maximum);

  fileout.mkdir("chip_canvas");
  fileout.cd("chip_canvas");
  TCanvas *C_MIP_high[layerNu][chipNu];
  // TCanvas *C_MIP_low[layerNu][chipNu];
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      // sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_MIP_high[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_high", i_layer, i_chip), Form("layer%d_chip%d_high", i_layer, i_chip), 2560, 1440);
      C_MIP_high[i_layer][i_chip]->Divide(6, 6);
      // C_MIP_low[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_low", i_layer, i_chip), Form("layer%d_chip%d_low", i_layer, i_chip), 2560, 1440);
      // C_MIP_low[i_layer][i_chip]->Divide(6, 6);
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if(i_chip==5&&i_channel>29) continue;
        int range_l = thresholdValue[i_layer][i_chip][i_channel]+pedMeanHigh[i_layer][i_chip][i_channel] - 50;
        int range_h = thresholdValue[i_layer][i_chip][i_channel]+pedMeanHigh[i_layer][i_chip][i_channel] + 300;
        C_MIP_high[i_layer][i_chip]->cd(i_channel + 1);
        int fit_entry = HG_histo[i_layer][i_chip][i_channel]->Integral(range_l,range_h);
        // HG_histo[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(range_l,range_h);
        // HG_histo[i_layer][i_chip][i_channel]->Rebin(2);
        // HG_histo[i_layer][i_chip][i_channel]->Draw();
        gPad-> SetLogy();
        // int bin0pe = HG_histo[i_layer][i_chip][i_channel]->GetMaximumBin();
        // int amp0pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe);
        // int amp1pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+25);
        // int amp2pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+50);
        int init_gain = i_layer<4 || i_layer>27 ? 25 : 15;
        // int bin0pe = pedMeanHigh[i_layer][i_chip][i_channel];
        int bin0pe = thresholdValue[i_layer][i_chip][i_channel] + pedMeanHigh[i_layer][i_chip][i_channel];
        int amp0pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe);
        int amp1pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+init_gain);
        int amp2pe = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+init_gain);
        int maximum_entry = HG_histo[i_layer][i_chip][i_channel]->GetBinContent(HG_histo[i_layer][i_chip][i_channel]->GetMaximumBin());
        // highfun->SetParLimits(1,bin0pe-10,bin0pe+10);
        highfun->SetParLimits(1,bin0pe-20,bin0pe+20);
        highfun->SetParLimits(0,0,maximum_entry*1.5);
        highfun->SetParLimits(3,0,maximum_entry*1.5);
        highfun->SetParLimits(6,0,maximum_entry*1.5);
        if(i_layer<4 || i_layer>27) highfun->SetParameters(amp0pe,bin0pe,4,amp1pe,25,4,amp2pe,4);
        else highfun->SetParameters(amp0pe,bin0pe,4,amp1pe,15,4,amp2pe,4);
        HG_histo[i_layer][i_chip][i_channel]->Fit("highfun","QR","",bin0pe-init_gain,bin0pe+4*init_gain);
        // double fit_range_low = calc_intercept1(highfun->GetParameter(0), highfun->GetParameter(1), highfun->GetParameter(2), -1);
        // double fit_range_high = calc_intercept1(highfun->GetParameter(6), highfun->GetParameter(1)+2*highfun->GetParameter(4), highfun->GetParameter(7), 1);
        double fit_range_low = highfun->GetParameter(1) - highfun->GetParameter(2) * 2;
        double fit_range_high = highfun->GetParameter(1)+2*highfun->GetParameter(4) + highfun->GetParameter(7) * 2;
        HG_histo[i_layer][i_chip][i_channel]->Fit("highfun","QR","",fit_range_low,fit_range_high);
        if(i_layer<30) gain->Fill(highfun->GetParameter(4));
        if(i_layer<4 || (i_layer>27 && i_layer<30)) {
          gain_15->Fill(highfun->GetParameter(4));
          // gain_15_ratio->Fill(highfun->GetParameter(4) / SPSgain[i_layer][i_chip][i_channel]);
        }
        if(i_layer>=4 && i_layer<=27) {
          gain_10->Fill(highfun->GetParameter(4));
          // gain_10_ratio->Fill(highfun->GetParameter(4) / SPSgain[i_layer][i_chip][i_channel]);
        }
        gain_layer->Fill(highfun->GetParameter(4), i_layer);
        sigma_layer->Fill(highfun->GetParameter(2), i_layer);
        sigma_layer->Fill(highfun->GetParameter(5), i_layer);
        sigma_layer->Fill(highfun->GetParameter(7), i_layer);

        double p0 = highfun->GetParameter(0) /2.50662827463/ highfun->GetParameter(2) / (double)HG_histo[i_layer][i_chip][i_channel]->GetEntries();
        double p1 = highfun->GetParameter(3) /2.50662827463/ highfun->GetParameter(5) / (double)HG_histo[i_layer][i_chip][i_channel]->GetEntries();
        double mean = HG_histo[i_layer][i_chip][i_channel]->GetMean();

        gStyle->SetOptFit(1111);
        HG_histo[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(range_l,range_h);
        HG_histo[i_layer][i_chip][i_channel]->Draw();

        _cellID_fit = i_layer*1e5 + i_chip*1e4 + i_channel;
        _gain = highfun->GetParameter(4);
        _gainError = highfun->GetParError(4);
        _gausSigma = highfun->GetParameter(2);
        _gausSigmaError = highfun->GetParError(2);
        _chiSquare = highfun->GetChisquare();
        _NDF = highfun->GetNDF();
        _temp = channel_temperature_mean[i_layer][i_chip][i_channel];
        _tempError = channel_temperature_sigma[i_layer][i_chip][i_channel];
        _fitting = highfun;
        _ctap = ctap_calc(p0,p1);
        _eqf = ctap_eqf(p0, highfun->GetParameter(4), mean);
        _fit_entry = fit_entry;
        gain_tree->Fill();
     
        // C_MIP_low[i_layer][i_chip]->cd(i_channel + 1);
        // gStyle->SetOptFit(1111);
        // LG_histo[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(range_l,range_h);
        // LG_histo[i_layer][i_chip][i_channel]->Draw();
        // gPad-> SetLogy();
      }
      C_MIP_high[i_layer][i_chip]->Update();
      C_MIP_high[i_layer][i_chip]->Write();
      C_MIP_high[i_layer][i_chip]->SaveAs(Form("%s/histo/high_Layer%d_Chip%d.png",picDirName.c_str(),i_layer,i_chip));
      // C_MIP_low[i_layer][i_chip]->Write();
      // C_MIP_low[i_layer][i_chip]->SaveAs(Form("%s/low_Layer%d_Chip%d.png",picDirName.c_str(),i_layer,i_chip));

      delete C_MIP_high[i_layer][i_chip];
      // delete C_MIP_low[i_layer][i_chip];
    }
  }
  fileout.cd();
  fileout.mkdir("raw_histogram");
  fileout.cd("raw_histogram");
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        LG_histo[i][j][k]->Write();
        HG_histo[i][j][k]->Write();
        delete LG_histo[i][j][k];
        delete HG_histo[i][j][k];
      }
    }
  }
  fileout.cd();

  gain->Write();
  gain_layer->Write();
  sigma_layer->Write();
  gain_10->Write();
  gain_15->Write();
  gain_10_ratio->Write();
  gain_15_ratio->Write();



  cout << "saved gain histograms" << endl;



  fileout.cd();
  fileout.mkdir("temperature_histo");
  fileout.cd("temperature_histo");
  TCanvas *C_temperature[layerNu];
  TF1 *temp_fun = new TF1("temp_fun","gaus",20,30);
  TH2F *temp_sensor = new TH2F("temp_sensor","temperature;layer;sensor;temperature [#circC]",32,0,32,16,0,16);
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    C_temperature[i_layer] = new TCanvas(Form("layer%d_temperature", i_layer), Form("layer%d_temperature", i_layer), 2560, 1440);
    C_temperature[i_layer]->Divide(4, 4);
    for(int i_sensor=0;i_sensor<16;i_sensor++){
      C_temperature[i_layer]->cd(i_sensor + 1);
      sensor_temperature[i_layer][i_sensor]->Draw();
      sensor_temperature[i_layer][i_sensor]->Fit("temp_fun","QR","",20,30);
      sensor_temperature[i_layer][i_sensor]->Write();
      gPad-> SetLogy();

      temp_sensor->SetBinContent(i_layer+1,i_sensor+1, temp_fun->GetParameter(1));
    }
    C_temperature[i_layer]->Write();
    for(int i_sensor=0;i_sensor<16;i_sensor++) delete sensor_temperature[i_layer][i_sensor];
  }
  temp_sensor->Write();
  for(int i_layer=0;i_layer<layerNu;i_layer++) delete C_temperature[i_layer];
  delete temp_sensor;


  TCanvas *C_temperature_channel[layerNu][chipNu];
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      // sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_temperature_channel[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_temperature", i_layer, i_chip), Form("layer%d chip%d temperature", i_layer, i_chip), 2560, 1440);
      C_temperature_channel[i_layer][i_chip]->Divide(6, 6);
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if(i_chip==5&&i_channel>29) continue;
        C_temperature_channel[i_layer][i_chip]->cd(i_channel + 1);
        gPad-> SetLogy();
        // channel_temperature[i_layer][i_chip][i_channel]->Fit("temp_channel_fun","QR","",20,30);
        channel_temperature[i_layer][i_chip][i_channel]->Draw();
      }
      C_temperature_channel[i_layer][i_chip]->Write();
      C_temperature_channel[i_layer][i_chip]->SaveAs(Form("%s/temp/histo/Layer%d_Chip%d.png",picDirName.c_str(),i_layer,i_chip));
    }
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      delete C_temperature_channel[i_layer][i_chip];
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) delete channel_temperature[i_layer][i_chip][i_channel];
    }
  }
  TCanvas *C_temperature_graph_channel[layerNu][chipNu];
  TGraphErrors *channel_temperature_graph[layerNu][chipNu][channelNu];
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      // sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_temperature_graph_channel[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_temperature_graph", i_layer, i_chip), Form("layer%d chip%d temperature", i_layer, i_chip), 2560, 1440);
      gStyle->SetPadLeftMargin(0.15);   // 左余白（Y軸タイトル用）
      gStyle->SetPadBottomMargin(0.15); // 下余白（X軸タイトル・時刻用）
      gStyle->SetPadRightMargin(0.05);  // 右余白
      gStyle->SetPadTopMargin(0.10);    // 上余白
      gStyle->SetTitleY(0.98);
      gStyle->SetTitleFontSize(0.07);
      C_temperature_graph_channel[i_layer][i_chip]->Divide(6, 6);
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if(i_chip==5&&i_channel>29) continue;

        int c_idx = i_layer * 216 + i_chip * 36 + i_channel;
        vector<double> x_pts, y_pts, x_errs, y_errs;
        for (size_t bin_idx = 0; bin_idx < cell_data[c_idx].size(); ++bin_idx) {
          double mean = 0.0, sem = 0.0;
          if (!cell_data[c_idx][bin_idx].get_results(mean, sem)) continue;
          double bin_center_sec = (bin_idx + 0.5) * timeWindow;
          x_pts.push_back(bin_center_sec);
          y_pts.push_back(mean);
          x_errs.push_back(timeWindow / 2.0);
          y_errs.push_back(sem);
        }
        if (x_pts.empty()) continue;
        channel_temperature_graph[i_layer][i_chip][i_channel] = new TGraphErrors(x_pts.size(), &x_pts[0], &y_pts[0], &x_errs[0], &y_errs[0]);

        C_temperature_graph_channel[i_layer][i_chip]->cd(i_channel + 1);
        // gPad-> SetLogy();
        // channel_temperature[i_layer][i_chip][i_channel]->Fit("temp_channel_fun","QR","",20,30);
        channel_temperature_graph[i_layer][i_chip][i_channel]->SetTitle(Form("Layer %d Chip %d Channel %d;;temperature [#circC]", i_layer, i_chip, i_channel));
        channel_temperature_graph[i_layer][i_chip][i_channel]->SetMarkerSize(0.5);
        channel_temperature_graph[i_layer][i_chip][i_channel]->SetMarkerStyle(20);
        TAxis *xaxis = channel_temperature_graph[i_layer][i_chip][i_channel]->GetXaxis();
        xaxis->SetTimeDisplay(1);                         // 時刻描画を有効化
        xaxis->SetTimeFormat("#splitline{%m/%d}{%H:%M}"); // フォーマット（上段: 月/日, 下段: 時:分）
        xaxis->SetTimeOffset(global_start_time, "gmt");   // X=0秒の位置を基準時刻に指定
        xaxis->SetTitleSize(0.07);   // 軸タイトルの文字サイズ (デフォルトは 0.04 前後)
        xaxis->SetLabelSize(0.055);  // 目盛り数字（または時刻）の文字サイズ
        xaxis->SetTitleOffset(1.1);  // 軸タイトルと目盛りの距離（重なり防止
        xaxis->SetLabelOffset(0.02);
        TAxis *yaxis = channel_temperature_graph[i_layer][i_chip][i_channel]->GetYaxis();
        yaxis->SetTitleSize(0.07);   // 軸タイトルの文字サイズ
        yaxis->SetLabelSize(0.055);  // 目盛り数字の文字サイズ
        yaxis->SetTitleOffset(1.2);  // Y軸タイトルが数値と被らないよう少し離す
        channel_temperature_graph[i_layer][i_chip][i_channel]->Draw("AP");
      }
      C_temperature_graph_channel[i_layer][i_chip]->Write();
      C_temperature_graph_channel[i_layer][i_chip]->SaveAs(Form("%s/temp/graph/Layer%d_Chip%d.png",picDirName.c_str(),i_layer,i_chip));
    }
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if(channel_temperature_graph[i_layer][i_chip][i_channel] != nullptr) delete channel_temperature_graph[i_layer][i_chip][i_channel];
      }
      delete C_temperature_graph_channel[i_layer][i_chip];
    }
  }

  cout << "saved temperature" << endl;

  fileout.cd();
  gain_tree->Write();


}