#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
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
#include <cmath>
#include <limits>
#include <numeric>
#include "TCanvas.h"
#include "langaus.C"        // langaus(chargeH, &fitFunc, &peakP, &peakPError);
#include "EBUdecode.h"
#include "EBUdecode.cxx"
// #include "/home/jap/analyseCode/Diagnose/src/EBUdecode.cxx"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// LEDデータをまとめるコード
// 1runのLEDのデータをまとめるコード

int dacList[9]={2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500};


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

bool best_fitting_amp(double n0, double n1, double n2){
  double amp_range[2] = {0.4, 2};
  bool res = amp_range[0] < n0/n1 && n0/n1 < amp_range[1];
  res = res && (amp_range[0] < n2/n1 && n2/n1 < amp_range[1]);
  return res;
}
bool best_fitting_gausSigma(double n0, double n1, double n2){
  double sigma_range[2] = {0.5, 1.5};
  bool res = sigma_range[0] < n1/n0 && n1/n0 < sigma_range[1];
  res = res && (sigma_range[0] < n2/n1 && n2/n1 < sigma_range[1]);
  return res;
}
bool best_fitting_chi2(double chi2_ndf){
  double chi_range[2] = {0.5, 1.5};
  bool res = chi_range[0]< chi2_ndf && chi2_ndf < chi_range[1];
  return res;
}

bool best_fitting_amp1(double n0, double n1, double n2){
  double amp_range[2] = {0.2, 3};
  bool res = amp_range[0] < n0/n1 && n0/n1 < amp_range[1];
  res = res && (amp_range[0] < n2/n1 && n2/n1 < amp_range[1]);
  return res;
}
bool best_fitting_gausSigma1(double n0, double n1, double n2){
  double sigma_range[2] = {0.5, 1.5};
  bool res = sigma_range[0] < n1/n0 && n1/n0 < sigma_range[1];
  res = res && (sigma_range[0] < n2/n1 && n2/n1 < sigma_range[1]);
  return res;
}
bool best_fitting_chi21(double chi2_ndf){
  double chi_range[2] = {0.5, 1.5};
  bool res = chi_range[0]< chi2_ndf && chi2_ndf < chi_range[1];
  return res;
}



int main(int argc, char* argv[])
{ 
    if(argc < 4){                                                     //エラー処理
        cout << "usage: ./between_files output.root input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  for(int i=0; i<rawfilenum; i++){
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("T_Event");
    if(tree[i]!=nullptr) entry_max[i] = tree[i]->GetEntries();
    else entry_max[i] = 0;
  }

  int dac[rawfilenum];
  for(int i=0; i<rawfilenum; i++){
    string txtname(argv[i+2]);
    txtname.erase(0,txtname.find_last_of("_")+1);
    // txtname.erase(txtname.find_last_of("_"),txtname.end()-1);
    dac[i] = stoi(txtname);
    // cout << argv[i+2] << ",  " << stoi(txtname) << endl;
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





  double cycleID, triggerID;
  vector<int> *cellIDs = nullptr;
  vector<int> *BCIDs = nullptr;
  vector<int> *hitTags = nullptr;
  vector<int> *gainTags = nullptr;
  vector<double> *charges = nullptr;
  vector<double> *times = nullptr;
  vector<vector<double>> *temp = nullptr;


  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;


  // const int fit_sigma_maximum = 9;
  const int fit_sigma_maximum = 6;



  TFile fileout(argv[1],"RECREATE");
  // gain_histo->Write();



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



  // dataを入れるもの
  TCanvas *canvas_highGain[layerNu][chipNu][channelNu];
  TCanvas *canvas_highGain_chip[layerNu][chipNu][rawfilenum];
  TCanvas *canvas_chip[layerNu][chipNu];
  TH1F *raw_histos[layerNu][chipNu][channelNu][rawfilenum];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        canvas_highGain[i][j][k] = new TCanvas(Form("canvas_highGain_%d_%d_%d",i,j,k),Form("canvas_highGain_%d_%d_%d",i,j,k),1);
        canvas_highGain[i][j][k]->Divide(3,3);
      }
      for(int irawfile=0; irawfile<rawfilenum; irawfile++){
        canvas_highGain_chip[i][j][irawfile] = new TCanvas(Form("canvas_highGain_chip_%d_%d_%d",i,j,irawfile),Form("canvas_highGain_chip_%d_%d_%d",i,j,irawfile),1);
        canvas_highGain_chip[i][j][irawfile]->Divide(6,6);
      }
      canvas_chip[i][j] = new TCanvas(Form("canvas_chip_%d_%d",i,j),Form("canvas_chip_%d_%d",i,j),1);
      canvas_chip[i][j]->Divide(6,6);
    }
  }
  cout << "made save TTree" << endl;

  double chi2_ndf_scan[layerNu][chipNu][channelNu][rawfilenum];
  double gain_scan[layerNu][chipNu][channelNu][rawfilenum];
  // double gainError_scan[layerNu][chipNu][channelNu][rawfilenum];
  double gausSigma_scan_n0[layerNu][chipNu][channelNu][rawfilenum];
  double gausSigma_scan_n1[layerNu][chipNu][channelNu][rawfilenum];
  double gausSigma_scan_n2[layerNu][chipNu][channelNu][rawfilenum];
  // double chiSquare_scan[layerNu][chipNu][channelNu][rawfilenum];
  // double NDF_scan[layerNu][chipNu][channelNu][rawfilenum];
  // double temp_scan[layerNu][chipNu][channelNu][rawfilenum];
  // double tempError_scan[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n0[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n1[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n2[layerNu][chipNu][channelNu][rawfilenum];
  // double ctap_scan[layerNu][chipNu][channelNu][rawfilenum];
  // double eqf_scan[layerNu][chipNu][channelNu][rawfilenum];
  double fit_entry_scan[layerNu][chipNu][channelNu][rawfilenum];
  double max_chi2_ndf[layerNu][chipNu][channelNu];
  int max_chi2_ndf_dac[layerNu][chipNu][channelNu];
  int max_chi2_ndf_dac1[layerNu][chipNu][channelNu];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        for(int irawfile=0; irawfile<rawfilenum; irawfile++){
          chi2_ndf_scan[i][j][k][irawfile] = 0;
          gain_scan[i][j][k][irawfile] = 0;
          // gainError_scan[i][j][k][irawfile] = 0;
          gausSigma_scan_n0[i][j][k][irawfile] = 0;
          gausSigma_scan_n1[i][j][k][irawfile] = 0;
          gausSigma_scan_n2[i][j][k][irawfile] = 0;
          // chiSquare_scan[i][j][k][irawfile] = 0;
          // NDF_scan[i][j][k][irawfile] = 0;
          // temp_scan[i][j][k][irawfile] = 0;
          // tempError_scan[i][j][k][irawfile] = 0;
          fitting_scan_n0[i][j][k][irawfile] = 0;
          fitting_scan_n1[i][j][k][irawfile] = 0;
          fitting_scan_n2[i][j][k][irawfile] = 0;
          // ctap_scan[i][j][k][irawfile] = 0;
        }
        max_chi2_ndf[i][j][k] = 1000;
        max_chi2_ndf_dac[i][j][k] = -1;
        max_chi2_ndf_dac1[i][j][k] = -1;
      }
    }
  }
  cout << "made save TTree" << endl;

  double _gain, _gainError, _gausSigma, _chiSquare, _NDF, _temp, _tempError, _ctap, _eqf;
  int _cellIDs, _fit_entry;
  TF1 *_fitting = nullptr;

  // // // // TF1 *highfun = new TF1("highfun","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1]-[4])/[5])^2) + [6]*exp(-0.5*((x-[1]-[4]*2)/[7])^2)",300,600);


  double channel_temperature_mean[32][6][36] = {0};
  double channel_temperature_sigma[32][6][36] = {0};





  // double out_gain, out_gainError, out_chiSquare, out_NDF, out_temp, out_tempError, out_ctap, out_eqf, out_sigma0, out_sigma1, out_sigma2, out_amp0, out_amp1, out_amp2;
  // int out_cellIDs, out_fit_entry, out_dac;

  // TTree* gain_tree = new TTree("T_Event", "gain calibration");
  // gain_tree->Branch("dac", &out_dac);
  // gain_tree->Branch("cellIDs", &out_cellIDs);
  // gain_tree->Branch("gain", &out_gain);
  // // gain_tree->Branch("gainError", &out_gainError);
  // gain_tree->Branch("chiSquare", &out_chiSquare);
  // gain_tree->Branch("NDF", &out_NDF);
  // gain_tree->Branch("temp", &out_temp);
  // gain_tree->Branch("sigma0", &out_sigma0);
  // gain_tree->Branch("sigma1", &out_sigma1);
  // gain_tree->Branch("sigma2", &out_sigma2);
  // gain_tree->Branch("amp0", &out_amp0);
  // gain_tree->Branch("amp1", &out_amp1);
  // gain_tree->Branch("amp2", &out_amp2);
  // gain_tree->Branch("tempError", &tempError);
  // gain_tree->Branch("fitting", &fitting);
  // gain_tree->Branch("ctap", &ctap);
  // gain_tree->Branch("eqf", &eqf);
  // gain_tree->Branch("fit_entry", &out_fit_entry);

  cout << "made save TTree" << endl;





  double SPSgain[layerNu][chipNu][channelNu];
  double SPSgain_error[layerNu][chipNu][channelNu];
  TH1F *gain_histo = new TH1F("gain_histo","gain",120,10,40);
  // TGraphErrors *gain_graph = new TGraphErrors();
  // int count_gain_graph = 0;

  int _SPScellID;
  double _SPSgain, _SPSgainError;
  // TFile *gainfile = new TFile("/megraid01/users/tsuji/scecal/ScECAL_CR/lightCalib/analyseCode/CaliSPS/macro/SPSgain.root");
  TFile *gainfile = new TFile("SPSgain.root");
  TTree *t_gainfile = (TTree*) gainfile->Get("T_Event");
  t_gainfile->SetBranchAddress("cellID",&_SPScellID);
  t_gainfile->SetBranchAddress("gain",&_SPSgain);
  t_gainfile->SetBranchAddress("gainError",&_SPSgainError);
  int SPSentry = 0;
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        if(j==5&&k>29) continue;
        t_gainfile->GetEntry(SPSentry);
        int cellid_now = (30+i)*1000000 + j*10000 + k;
        //if(cellid_now!=_SPScellID){
        //  cout << i << ", " << j << ", " << k << "   " << _SPScellID << " " << cellid_now;
        //  cout << "    something is wrong at gain search" << endl;
        //}
        SPSgain[i][j][k] = _SPSgain;
        SPSgain_error[i][j][k] = _SPSgainError;
        SPSentry++;
        //cout << i << ", " << j << ", " << k << "   gain : "  << _SPSgain << "  error : " << _SPSgainError << endl;
        gain_histo->Fill(_SPSgain);
        // gain_graph->SetPoint(count_gain_graph,count_gain_graph,_SPSgain);
        // gain_graph->SetPointError(count_gain_graph,0,_SPSgainError);
        // count_gain_graph++;
      }
    }
  }
  gainfile->Close();
  cout << "CR gain reading finished" << endl;
  fileout.cd();






  // // // // TF1 *highfun = new TF1("highfun","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1]-[4])/[5])^2) + [6]*exp(-0.5*((x-[1]-[4]*2)/[7])^2)",300,600);
  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

    tree[irawfile]->SetBranchAddress("cellIDs", &_cellIDs);
    tree[irawfile]->SetBranchAddress("gain", &_gain);
    tree[irawfile]->SetBranchAddress("gainError", &_gainError);
    tree[irawfile]->SetBranchAddress("chiSquare", &_chiSquare);
    tree[irawfile]->SetBranchAddress("NDF", &_NDF);
    tree[irawfile]->SetBranchAddress("temp", &_temp);
    tree[irawfile]->SetBranchAddress("tempError", &_tempError);
    tree[irawfile]->SetBranchAddress("fitting", &_fitting);
    tree[irawfile]->SetBranchAddress("ctap", &_ctap);
    tree[irawfile]->SetBranchAddress("eqf", &_eqf);
    tree[irawfile]->SetBranchAddress("fit_entry", &_fit_entry);

    // out_dac = dacList[irawfile];

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);

      int _layer = _cellIDs/1e6;
      int _chip = (_cellIDs/10000)%100;
      int _channel = _cellIDs%100;

      // _fitting->GetParameter();

      // out_cellIDs = _cellIDs;
      // out_gain = _gain;
      // out_chiSquare = _chiSquare;
      // out_NDF = _NDF;
      // out_temp = _temp;
      // out_sigma0 = _fitting->GetParameter(2);
      // out_sigma1 = _fitting->GetParameter(5);
      // out_sigma2 = _fitting->GetParameter(7);
      // out_amp0 = _fitting->GetParameter(0);
      // out_amp1 = _fitting->GetParameter(3);
      // out_amp2 = _fitting->GetParameter(6);
      // gain_tree->Fill();

      chi2_ndf_scan[_layer][_chip][_channel][irawfile] = _NDF!=0 ? _chiSquare / _NDF : 1000;
      gain_scan[_layer][_chip][_channel][irawfile] = _gain;
      // gainError_scan[_layer][_chip][_channel][irawfile] = _gainError;
      gausSigma_scan_n0[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(2);
      gausSigma_scan_n1[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(5);
      gausSigma_scan_n2[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(7);
      // chiSquare_scan[_layer][_chip][_channel][irawfile] = _chiSquare;
      // NDF_scan[_layer][_chip][_channel][irawfile] = _NDF;
      // temp_scan[_layer][_chip][_channel][irawfile] = _temp;
      // tempError_scan[_layer][_chip][_channel][irawfile] = _tempError;
      fitting_scan_n0[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(0);
      fitting_scan_n1[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(3);
      fitting_scan_n2[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(6);
      // ctap_scan[_layer][_chip][_channel][irawfile] = _ctap;
      // eqf_scan[_layer][_chip][_channel][irawfile] = _eqf;
      fit_entry_scan[_layer][_chip][_channel][irawfile] = _fit_entry;


    }

    for(int i=0; i<layerNu;i++){
      for(int j=0; j<chipNu;j++){
        for(int k=0; k<channelNu;k++){
          raw_histos[i][j][k][irawfile] = (TH1F*)filein[irawfile]->Get(Form("raw_histogram/time_%d_%d_%d",i,j,k));
          raw_histos[i][j][k][irawfile]->GetXaxis()->SetRangeUser(350, 550);
          TF1 *ff = raw_histos[i][j][k][irawfile]->GetFunction("highfun");
          if(ff!=nullptr) ff->SetNpx(1000);
          raw_histos[i][j][k][irawfile]->SetTitle(Form("%d_%d_%d dac %d", i,j,k, dac[irawfile]));
          // canvas_highGain[i][j][k]->cd();
          // canvas_highGain[i][j][k]->cd(irawfile+1);
          // raw_histos[i][j][k][irawfile]->Draw();

          // canvas_highGain_chip[i][j][irawfile]->cd();
          // canvas_highGain_chip[i][j][irawfile]->cd(k+1);
          // raw_histos[i][j][k][irawfile]->Draw();
        }
      }
    }

    // delete filein[irawfile];
  }

  // gain_tree->Write();



  // best fit の探索
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        if(j==5&&k>29) continue;
        for(int irawfile=0; irawfile<rawfilenum; irawfile++){
          double chi2_ndf = chi2_ndf_scan[i][j][k][irawfile]>1 ? chi2_ndf_scan[i][j][k][irawfile] : 1.0/chi2_ndf_scan[i][j][k][irawfile];
          if(!best_fitting_chi2(chi2_ndf_scan[i][j][k][irawfile])) continue;
          if(fit_entry_scan[i][j][k][irawfile] < 10) continue;
          double n0 =  fitting_scan_n0[i][j][k][irawfile];
          double n1 =  fitting_scan_n1[i][j][k][irawfile];
          double n2 =  fitting_scan_n2[i][j][k][irawfile];
          // if(n0<2||n2<2) continue;
          if(!best_fitting_amp(n0,n1,n2)) continue;
          if(!best_fitting_gausSigma(gausSigma_scan_n0[i][j][k][irawfile],gausSigma_scan_n1[i][j][k][irawfile],gausSigma_scan_n2[i][j][k][irawfile])) continue;
          bool misfit = (i<4 || i>27) ? (gain_scan[i][j][k][irawfile] == 25) : (gain_scan[i][j][k][irawfile] == 15);
          if(misfit) continue;
          if(max_chi2_ndf_dac[i][j][k]==-1){          /////////////////////////// ここをchi2殻だけじゃなくて、tripleになってるかを判断するべき
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac[i][j][k] = irawfile;
            // cout << max_chi2_ndf[i][j][k] << " " << chi2_ndf << " " << chi2_ndf_scan[i][j][k][irawfile] << "  " << chi2_ndf_scan[i][j][k][irawfile] << "  " << 1.0/chi2_ndf_scan[i][j][k][irawfile] << endl;
            continue;
          }

          if(chi2_ndf < max_chi2_ndf[i][j][k]){
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac[i][j][k] = irawfile;
          }
        }
      }
    }
  }
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        if(j==5&&k>29) continue;
        for(int irawfile=0; irawfile<rawfilenum; irawfile++){
          if(max_chi2_ndf_dac[i][j][k]!=-1) max_chi2_ndf_dac1[i][j][k] = max_chi2_ndf_dac[i][j][k]; continue;

          double chi2_ndf = chi2_ndf_scan[i][j][k][irawfile]>1 ? chi2_ndf_scan[i][j][k][irawfile] : 1.0/chi2_ndf_scan[i][j][k][irawfile];
          if(!best_fitting_chi21(chi2_ndf_scan[i][j][k][irawfile])) continue;
          if(fit_entry_scan[i][j][k][irawfile] < 10) continue;
          double n0 =  fitting_scan_n0[i][j][k][irawfile];
          double n1 =  fitting_scan_n1[i][j][k][irawfile];
          double n2 =  fitting_scan_n2[i][j][k][irawfile];
          if(!best_fitting_amp1(n0,n1,n2)) continue;
          if(!best_fitting_gausSigma1(gausSigma_scan_n0[i][j][k][irawfile],gausSigma_scan_n1[i][j][k][irawfile],gausSigma_scan_n2[i][j][k][irawfile])) continue;
          bool misfit = (i<4 || i>27) ? (gain_scan[i][j][k][irawfile] == 25) : (gain_scan[i][j][k][irawfile] == 15);
          if(misfit) continue;
          if(max_chi2_ndf_dac1[i][j][k]==-1){          /////////////////////////// ここをchi2殻だけじゃなくて、tripleになってるかを判断するべき
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac1[i][j][k] = irawfile;
            // cout << max_chi2_ndf[i][j][k] << " " << chi2_ndf << " " << chi2_ndf_scan[i][j][k][irawfile] << "  " << chi2_ndf_scan[i][j][k][irawfile] << "  " << 1.0/chi2_ndf_scan[i][j][k][irawfile] << endl;
            continue;
          }

          if(chi2_ndf < max_chi2_ndf[i][j][k]){
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac1[i][j][k] = irawfile;
          }
        }
      }
    }
  }



  // cout << "saving raw histograms" << endl;
  // fileout.mkdir("raw_histogram");
  // fileout.cd("raw_histogram");
  // for(int i=0; i<layerNu;i++){
  //   cout << i << "/" << layerNu << endl;
  //   for(int j=0; j<chipNu;j++){
  //     for(int k=0; k<channelNu;k++){
  //       canvas_highGain[i][j][k]->Write();
  //       // delete canvas_highGain[i][j][k];
  //     }
  //   }
  // }
  // cout << "saved raw histograms" << endl;
  // fileout.cd();

  cout << "saving raw histograms chips" << endl;
  fileout.mkdir("raw_histogram_chip");
  fileout.cd("raw_histogram_chip");
  for(int i=0; i<layerNu;i++){
    cout << i << "/" << layerNu << endl;
    for(int j=0; j<chipNu;j++){
      canvas_chip[i][j]->cd();
      for(int k=0; k<channelNu;k++){
        canvas_chip[i][j]->cd(k+1);
        if(max_chi2_ndf_dac1[i][j][k]==-1) continue;
        raw_histos[i][j][k][max_chi2_ndf_dac1[i][j][k]]->Draw();
      }
      // for(int irawfile=0; irawfile<rawfilenum; irawfile++){
      //   canvas_highGain_chip[i][j][irawfile]->Write();
      //   // delete canvas_highGain[i][j][k];
      // }
      canvas_chip[i][j]->Write();
    }
  }
  cout << "saved raw histograms chips" << endl;
  // fileout.cd();
  
  // return 1;



}