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
#include <algorithm>
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
// 全部のrunのデータからgainの温度依存等を求めるコード



int main(int argc, char* argv[])
{
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files output.root input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;

  int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  vector<pair<int,int>> filename;
  int total_entry_max=0;

  double rawFileDate[rawfilenum];

  // for(int i=0; i<rawfilenum; i++){
  //   string txtname(argv[i+2]);
  //   txtname.erase(0,txtname.find_last_of("p")+1);
  //   // txtname.erase(txtname.find_last_of("_"),txtname.end()-1);
  //   filename.push_back({stoi(txtname), i+2});
  //   // cout << argv[i+2] << ",  " << stoi(txtname) << endl;
  // }
  // sort(filename.begin(),filename.end());
  for(int i=0; i<rawfilenum; i++){
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("T_Event");
    if(tree[i]!=nullptr) entry_max[i] = tree[i]->GetEntries();
    else entry_max[i] = 0;
  }

    // /megraid01/users/tsuji/scecal/led/data/analyseCode/LEDCalibration/LEDgroupID.txt にある
  int group_channel[14][15] = {
    {2 ,      6 ,      8 ,      1  ,     5  ,     202,     190,     192,     194,     179,     175,     187,     176,     173,     171},
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




  double _gain, _gainError, _gausSigma, _chiSquare, _NDF, _temp, _tempError;
  int _cellIDs;




  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;


  const double temp_max=30, temp_min=21, gain_max_10=20, gain_min_10=10, gain_max_15=30, gain_min_15=15;


  const int fit_sigma_maximum = 9;

  double channel_temperature[32][6][36] = {0};


  vector<vector<int>> dead_channel{{0, 3, 28}, {0, 5, 3}, {0, 5, 4}, {1, 4, 15}, {2, 0, 7}, {2, 3, 32}, {2, 4, 18}, {2, 5, 0}, {3, 0, 8}, {4, 1, 12}, {5, 0, 22}, {5, 1, 35}, {5, 2, 3}, {6, 3, 19}, {6, 3, 21}, 
    {7, 3, 27}, {8, 5, 19}, {9, 0, 15}, {10, 2, 2}, {10, 2, 24}, {11, 3, 29}, {12, 2, 24}, {12, 2, 26}, {12, 5, 1}, {12, 5, 2}, {13, 3, 17}, {13, 3, 19}, {14, 4, 17}, {15, 3, 0}, {15, 4, 11}, {16, 1, 5}, 
    {16, 2, 16}, {16, 5, 10}, {17, 0, 3}, {17, 1, 9}, {17, 5, 21}, {18, 1, 11}, {18, 1, 15}, {18, 5, 10}, {19, 1, 30}, {20, 1, 35}, {20, 2, 6}, {21, 0, 2}, {21, 1, 27}, {21, 1, 32}, {22, 4, 9}, {23, 3, 27}, 
    {23, 3, 29}, {28, 0, 21}, {28, 0, 32}, {28, 1, 25}, {28, 1, 29}, {28, 2, 6}, {28, 2, 26}, {28, 3, 31}, {28, 4, 15}, {28, 4, 23}, {29, 2, 23}, {29, 3, 19}, {29, 5, 12}, {30, 1, 23}, {30, 3, 27}, 
    {30, 5, 15}, {31, 3, 18} };


  TFile fileout(argv[1],"RECREATE");
  // gain_histo->Write();



  // dataを入れるもの
  TGraphErrors *gain_temp[layerNu][chipNu][channelNu];
  int gain_temp_count[layerNu][chipNu][channelNu];
  TGraphErrors *gain_temp_cut[layerNu][chipNu][channelNu];
  int gain_temp_cut_count[layerNu][chipNu][channelNu];
  // TH1F *time_histo[layerNu][chipNu][channelNu];
  // int charge_count[layerNu][chipNu][channelNu];

  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        gain_temp[i][j][k] = new TGraphErrors();
        gain_temp[i][j][k]->SetTitle(Form("gain vs temp channel%d",k));
        gain_temp[i][j][k]->GetXaxis()->SetTitle("temperature [C]");
        gain_temp[i][j][k]->GetYaxis()->SetTitle("gain [ADC]");
        gain_temp_count[i][j][k] = 0;

        gain_temp_cut[i][j][k] = new TGraphErrors();
        gain_temp_cut[i][j][k]->SetTitle(Form("gain vs temp channel%d",k));
        gain_temp_cut[i][j][k]->GetXaxis()->SetTitle("temperature [C]");
        gain_temp_cut[i][j][k]->GetYaxis()->SetTitle("gain [ADC]");
        gain_temp_cut_count[i][j][k] = 0;
        // time_histo[i][j][k] = new TH1F(Form("time_%d_%d_%d",i,j,k),Form("time_%d_%d_%d",i,j,k),4096,0,4095);
        // charge_count[i][j][k] = 0;
      }
    }
  }

  TH1F *sensor_temperature[layerNu][16];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<16;j++){
      sensor_temperature[i][j] = new TH1F(Form("sensor_temperature_%d_%d",i,j),Form("temperature layer%d sensor%d",i,j),1000,20,30);
      sensor_temperature[i][j]->SetXTitle("temperature [C]");
    }
  }




  // runごとにつくるもの







  // 保存形式
  double gain, gainError, gausSigma, chiSquare, NDF, temp, tempError;
  int cellIDs, runNo;
  TTree* gain_tree = new TTree("T_Event", "gain calibration");
  gain_tree->Branch("runNo", &runNo);
  gain_tree->Branch("cellIDs", &cellIDs);
  gain_tree->Branch("gain", &gain);
  gain_tree->Branch("gainError", &gainError);
  gain_tree->Branch("chiSquare", &chiSquare);
  gain_tree->Branch("NDF", &NDF);
  gain_tree->Branch("temp", &temp);
  gain_tree->Branch("tempError", &tempError);






  cout << "made save TTree" << endl;





  double SPSgain[layerNu][chipNu][channelNu];
  double SPSgain_error[layerNu][chipNu][channelNu];
  TH1F *SPSgain_histo = new TH1F("SPSgain_histo","gain",120,10,40);
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
        SPSgain_histo->Fill(_SPSgain);
        // gain_graph->SetPoint(count_gain_graph,count_gain_graph,_SPSgain);
        // gain_graph->SetPointError(count_gain_graph,0,_SPSgainError);
        // count_gain_graph++;
      }
    }
  }
  gainfile->Close();
  cout << "CR gain reading finished" << endl;
  fileout.cd();


  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

    if(entry_max[irawfile]==0) continue;

    tree[irawfile]->SetBranchAddress("cellIDs", &_cellIDs);
    tree[irawfile]->SetBranchAddress("gain", &_gain);
    tree[irawfile]->SetBranchAddress("gainError", &_gainError);
    tree[irawfile]->SetBranchAddress("chiSquare", &_chiSquare);
    tree[irawfile]->SetBranchAddress("NDF", &_NDF);
    tree[irawfile]->SetBranchAddress("temp", &_temp);
    tree[irawfile]->SetBranchAddress("tempError", &_tempError);

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);

      int _layer = _cellIDs/1000000;
      int _chip = (_cellIDs/10000)%100;
      int _channel = _cellIDs%100;
      if(( ((_layer>3&&_layer<28) && (_gain>=gain_min_10 && _gain<=gain_max_10)) || (!(_layer>3&&_layer<28) && (_gain>=gain_min_15 && _gain<=gain_max_15)) ) && _gainError>=0.01) {
        gain_temp[_layer][_chip][_channel]->SetPoint(irawfile, _temp, _gain);
        gain_temp[_layer][_chip][_channel]->SetPointError(irawfile, _tempError, _gainError);
      }

      // if(_chiSquare/_NDF <1.5){
      if(_gainError <1.5){
        gain_temp_cut[_layer][_chip][_channel]->SetPoint(irawfile, _temp, _gain);
        gain_temp_cut[_layer][_chip][_channel]->SetPointError(irawfile, _tempError, _gainError);
      }

      cellIDs = _cellIDs;
      runNo = irawfile;
      gain = _gain;
      gainError = _gainError;
      chiSquare = _chiSquare;
      NDF = _NDF;
      temp = _temp;
      tempError = _tempError;
      gain_tree->Fill();
    }


    // _cellIDs = nullptr;
    // _gain = nullptr;
    // _gainError = nullptr;
    // _gausSigma = nullptr;
    // _chiSquare = nullptr;
    // _NDF = nullptr;
    // _temp = nullptr;
  }
  gain_tree->Write();


  fileout.mkdir("raw_histogram");
  fileout.cd("raw_histogram");
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        // gain_temp[i][j][k]->SetMarkerStyle(4);
        gain_temp[i][j][k]->Write(Form("gain_temp_%d_%d_%d",i,j,k));
        gain_temp_cut[i][j][k]->Write(Form("gain_temp_cut_%d_%d_%d",i,j,k));
        // time_histo[i][j][k]->Write();
        // charge_count[i][j][k] = 0;
        // cout << channel_temperature[i][j][k] << ", " << time_histo[i][j][k]->Integral() << endl;
      }
    }
  }
  fileout.cd();
  cout << "saved raw histograms" << endl;



  double slope, intercept, slope_error, intercept_error, NDF_temp, chiSquare_temp;
  int cellIDs_temp, layer_temp, chip_temp, channel_temp;
  TTree* temperatureDep = new TTree("temperatureDep", "gain calibration");
  temperatureDep->Branch("cellIDs", &cellIDs_temp);
  temperatureDep->Branch("layer", &layer_temp);
  temperatureDep->Branch("chip", &chip_temp);
  temperatureDep->Branch("channel", &channel_temp);
  temperatureDep->Branch("slope", &slope);
  temperatureDep->Branch("intercept", &intercept);
  temperatureDep->Branch("slope_error", &slope_error);
  temperatureDep->Branch("intercept_error", &intercept_error);
  temperatureDep->Branch("NDF", &NDF_temp);
  temperatureDep->Branch("chiSquare", &chiSquare_temp);

  double slope_cut, intercept_cut, slope_error_cut, intercept_error_cut, NDF_temp_cut, chiSquare_temp_cut;
  int cellIDs_temp_cut, layer_temp_cut, chip_temp_cut, channel_temp_cut;
  TTree* temperatureDep_cut = new TTree("temperatureDep_cut", "gain calibration");
  temperatureDep_cut->Branch("cellIDs", &cellIDs_temp_cut);
  temperatureDep_cut->Branch("layer", &layer_temp_cut);
  temperatureDep_cut->Branch("chip", &chip_temp_cut);
  temperatureDep_cut->Branch("channel", &channel_temp_cut);
  temperatureDep_cut->Branch("slope", &slope_cut);
  temperatureDep_cut->Branch("intercept", &intercept_cut);
  temperatureDep_cut->Branch("slope_error", &slope_error_cut);
  temperatureDep_cut->Branch("intercept_error", &intercept_error_cut);
  temperatureDep_cut->Branch("NDF", &NDF_temp_cut);
  temperatureDep_cut->Branch("chiSquare", &chiSquare_temp_cut);

  TF1 *fit = new TF1("fit","pol1",15,35);
  // fit->SetParLimits(1,-2,-0.005);
  // fit->SetParLimits(1,-2,0);
  TH2F *slope_channel = new TH2F("slope_channel","slope_channel", 192,0,192, 36,0,36);
  TH2F *slope_channel_normalized = new TH2F("slope_channel_normalized","slope_channel_normalized", 192,0,192, 36,0,36);
  TH2F *slope_channel_normalized_cut = new TH2F("slope_channel_normalized_cut","slope_channel_normalized_cut", 192,0,192, 36,0,36);
  slope_channel->SetZTitle("[ADC/C]");
  slope_channel_normalized->SetZTitle("[%/C]");
  slope_channel_normalized_cut->SetZTitle("[%/C]");

  TGraphErrors *Tdep_layer = new TGraphErrors();
  Tdep_layer->SetTitle(Form("temperature dependence of the gain at each layer"));
  Tdep_layer->GetXaxis()->SetTitle("layer");
  Tdep_layer->GetYaxis()->SetTitle("temperature dependence [%/C]");
  TGraphErrors *Tdep_chip = new TGraphErrors();
  Tdep_chip->SetTitle(Form("temperature dependence of the gain at each chip"));
  Tdep_chip->GetXaxis()->SetTitle("layer*6 + chip");
  Tdep_chip->GetYaxis()->SetTitle("temperature dependence [%/C]");
  double Tdep_layer_dep[32]={0}, Tdep_layer_dep_error[32]={0}, Tdep_layer_count[32]={0};
  double Tdep_chip_dep[192]={0}, Tdep_chip_dep_error[192]={0}, Tdep_hip_count[192]={0};
  
  
  fileout.mkdir("chip_canvas");
  fileout.cd("chip_canvas");
  TCanvas *C_MIP_high[layerNu][chipNu];
  TCanvas *C_MIP_high_cut[layerNu][chipNu];
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      // sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_MIP_high[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_gain", i_layer, i_chip), Form("layer%d_chip%d_gain", i_layer, i_chip), 1);
      C_MIP_high[i_layer][i_chip]->Divide(6, 6);
      C_MIP_high_cut[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_gain_cut", i_layer, i_chip), Form("layer%d_chip%d_gain_cut", i_layer, i_chip), 1);
      C_MIP_high_cut[i_layer][i_chip]->Divide(6, 6);
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        fit->SetParameters(30,-0.3);
        if(i_layer>3 && i_layer<28) fit->SetParameters(20,-0.02);
        if(i_chip==5&&i_channel>29) continue;
        vector<int> v{i_layer,i_chip,i_channel};
        if(find(dead_channel.begin(),dead_channel.end(),v) != dead_channel.end()) continue;
        C_MIP_high[i_layer][i_chip]->cd(i_channel + 1);
        gain_temp[i_layer][i_chip][i_channel]->GetXaxis()->SetLimits(temp_min,temp_max);
        if(i_layer>3 && i_layer<28) {
          gain_temp[i_layer][i_chip][i_channel]->SetMaximum(gain_max_10);
          gain_temp[i_layer][i_chip][i_channel]->SetMinimum(gain_min_10);
        } else {
          gain_temp[i_layer][i_chip][i_channel]->SetMaximum(gain_max_15);
          gain_temp[i_layer][i_chip][i_channel]->SetMinimum(gain_min_15);
        }
        gain_temp[i_layer][i_chip][i_channel]->Draw("AP");
        if(gain_temp[i_layer][i_chip][i_channel]->GetN()<5) continue;
        gain_temp[i_layer][i_chip][i_channel]->Fit("fit","QR","",temp_min,temp_max);
        cellIDs_temp = i_layer*1000000 + i_chip*10000 + i_channel;
        layer_temp = i_layer;
        chip_temp = i_chip;
        channel_temp = i_channel;
        slope = fit->GetParameter(1);
        intercept = fit->GetParameter(0);
        slope_error = fit->GetParError(1);
        intercept_error = fit->GetParError(0);
        NDF_temp = fit->GetNDF();
        chiSquare_temp = fit->GetChisquare();
        temperatureDep->Fill();

        double gain_20c = slope * 20.0 + intercept;

        slope_channel->SetBinContent(i_layer*6 + i_chip+1, i_channel+1, -slope);
        slope_channel_normalized->SetBinContent(i_layer*6 + i_chip+1, i_channel+1, -slope/gain_20c*100);
        if(-slope/gain_20c>0.0005) {
          slope_channel_normalized_cut->SetBinContent(i_layer*6 + i_chip+1, i_channel+1, -slope/gain_20c*100);

          double dependence_error = sqrt(pow(intercept*slope_error,2) + pow(slope*intercept_error,2)) * 100/gain_20c/gain_20c;

          Tdep_layer_dep[i_layer] += slope/gain_20c*100;
          Tdep_layer_dep_error[i_layer] += pow(dependence_error,2);
          Tdep_layer_count[i_layer]++;

          Tdep_chip_dep[i_layer*6+i_chip] += slope/gain_20c*100;
          Tdep_chip_dep_error[i_layer*6+i_chip] += pow(dependence_error,2);
          Tdep_hip_count[i_layer*6+i_chip]++;
        }


        C_MIP_high_cut[i_layer][i_chip]->cd(i_channel + 1);
        gain_temp_cut[i_layer][i_chip][i_channel]->Draw("AP");
        if(gain_temp_cut[i_layer][i_chip][i_channel]->GetN()>0){
          fit->SetParameters(40,-0.5);
          gain_temp_cut[i_layer][i_chip][i_channel]->Fit("fit","QR","",20,30);
          cellIDs_temp_cut = i_layer*1000000 + i_chip*10000 + i_channel;
          layer_temp_cut = i_layer;
          chip_temp_cut = i_chip;
          channel_temp_cut = i_channel;
          slope_cut = fit->GetParameter(1);
          intercept_cut = fit->GetParameter(0);
          slope_error_cut = fit->GetParError(1);
          intercept_error_cut = fit->GetParError(0);
          NDF_temp_cut = fit->GetNDF();
          chiSquare_temp_cut = fit->GetChisquare();
          temperatureDep_cut->Fill();
        }
        // gPad-> SetLogy();

      }
      C_MIP_high[i_layer][i_chip]->Write();
      C_MIP_high_cut[i_layer][i_chip]->Write();

      Tdep_chip->SetPoint(i_layer*6+i_chip,i_layer*6+i_chip, Tdep_chip_dep[i_layer*6+i_chip]/Tdep_hip_count[i_layer*6+i_chip]);
      Tdep_chip->SetPointError(i_layer*6+i_chip,0, sqrt(Tdep_chip_dep_error[i_layer*6+i_chip])/Tdep_hip_count[i_layer*6+i_chip]);
      // cout << Tdep_chip_dep[i_layer*6+i_chip]/Tdep_hip_count[i_layer*6+i_chip] << "," << sqrt(Tdep_chip_dep_error[i_layer*6+i_chip])/Tdep_hip_count[i_layer*6+i_chip] << endl;
    }
    Tdep_layer->SetPoint(i_layer,i_layer, Tdep_layer_dep[i_layer]/Tdep_layer_count[i_layer]);
    Tdep_layer->SetPointError(i_layer,0, sqrt(Tdep_layer_dep_error[i_layer])/Tdep_layer_count[i_layer]);
  }

  fileout.cd();
  temperatureDep->Write();
  temperatureDep_cut->Write();
  slope_channel->Write();
  slope_channel_normalized->Write();
  slope_channel_normalized_cut->Write();
  Tdep_layer->SetMarkerStyle(8);
  Tdep_layer->Write("Tdep_layer");
  Tdep_chip->SetMarkerStyle(8);
  Tdep_chip->Write("Tdep_chip");



  // cout << "saved gain histograms" << endl;






}