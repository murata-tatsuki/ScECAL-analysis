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
// dac scanの中から一番いいfitを見積もるマクロ


bool best_fitting(double n0, double n1, double n2){
  bool res = n0/n1 < 2.0 && n0/n1 > 0.5;
  res = res && (n2/n1 < 2.0 && n2/n1 > 0.5);
  return res;
}



int main(int argc, char* argv[])
{
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files output.root input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;

  int rawfilenum = argc - 2;
  if(rawfilenum != 9) cout << "something is wrong  file number is not 9" << endl;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  vector<pair<int,int>> filename;
  int total_entry_max=0;

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




  double _gain, _gainError, _gausSigma, _chiSquare, _NDF, _temp, _tempError, _ctap, _eqf;
  int _cellIDs, _fit_entry;
  TF1 *_fitting = nullptr;




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




  // runごとにつくるもの
  double chi2_ndf_scan[layerNu][chipNu][channelNu][rawfilenum];
  double gain_scan[layerNu][chipNu][channelNu][rawfilenum];
  double gainError_scan[layerNu][chipNu][channelNu][rawfilenum];
  double gausSigma_scan[layerNu][chipNu][channelNu][rawfilenum];
  double chiSquare_scan[layerNu][chipNu][channelNu][rawfilenum];
  double NDF_scan[layerNu][chipNu][channelNu][rawfilenum];
  double temp_scan[layerNu][chipNu][channelNu][rawfilenum];
  double tempError_scan[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n0[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n1[layerNu][chipNu][channelNu][rawfilenum];
  double fitting_scan_n2[layerNu][chipNu][channelNu][rawfilenum];
  double ctap_scan[layerNu][chipNu][channelNu][rawfilenum];
  double eqf_scan[layerNu][chipNu][channelNu][rawfilenum];
  double fit_entry_scan[layerNu][chipNu][channelNu][rawfilenum];
  double max_chi2_ndf[layerNu][chipNu][channelNu];
  int max_chi2_ndf_dac[layerNu][chipNu][channelNu];
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        for(int irawfile=0; irawfile<rawfilenum; irawfile++){
          chi2_ndf_scan[i][j][k][irawfile] = 0;
          gain_scan[i][j][k][irawfile] = 0;
          gainError_scan[i][j][k][irawfile] = 0;
          gausSigma_scan[i][j][k][irawfile] = 0;
          chiSquare_scan[i][j][k][irawfile] = 0;
          NDF_scan[i][j][k][irawfile] = 0;
          temp_scan[i][j][k][irawfile] = 0;
          tempError_scan[i][j][k][irawfile] = 0;
          fitting_scan_n0[i][j][k][irawfile] = 0;
          fitting_scan_n1[i][j][k][irawfile] = 0;
          fitting_scan_n2[i][j][k][irawfile] = 0;
          ctap_scan[i][j][k][irawfile] = 0;
        }
        max_chi2_ndf[i][j][k] = 1000;
        max_chi2_ndf_dac[i][j][k] = -1;
      }
    }
  }





  // 保存形式
  double gain, gainError, gausSigma, chiSquare, NDF, temp, tempError, ctap, eqf;
  int cellIDs, fit_entry, dac;  //, runNo;
  // TF1 *fitting = new TF1();
  TTree* gain_tree = new TTree("T_Event", "gain calibration");
  gain_tree->Branch("dac", &dac);
  gain_tree->Branch("cellIDs", &cellIDs);
  gain_tree->Branch("gain", &gain);
  gain_tree->Branch("gainError", &gainError);
  gain_tree->Branch("chiSquare", &chiSquare);
  gain_tree->Branch("NDF", &NDF);
  gain_tree->Branch("temp", &temp);
  gain_tree->Branch("tempError", &tempError);
  // gain_tree->Branch("fitting", &fitting);
  gain_tree->Branch("ctap", &ctap);
  gain_tree->Branch("eqf", &eqf);
  gain_tree->Branch("fit_entry", &fit_entry);






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
        int cellid_now = i*1000000 + j*10000 + k;
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


  TH1F *spSpectrum[layerNu][chipNu][channelNu];





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
    tree[irawfile]->SetBranchAddress("fitting", &_fitting);
    tree[irawfile]->SetBranchAddress("ctap", &_ctap);
    tree[irawfile]->SetBranchAddress("eqf", &_eqf);
    tree[irawfile]->SetBranchAddress("fit_entry", &_fit_entry);

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);

      int _layer = _cellIDs/1e6;
      int _chip = (_cellIDs/10000)%100;
      int _channel = _cellIDs%100;

      chi2_ndf_scan[_layer][_chip][_channel][irawfile] = _NDF!=0 ? _chiSquare / _NDF : 1000;
      gain_scan[_layer][_chip][_channel][irawfile] = _gain;
      gainError_scan[_layer][_chip][_channel][irawfile] = _gainError;
      // gausSigma_scan[_layer][_chip][_channel][irawfile] = _chiSquare;
      chiSquare_scan[_layer][_chip][_channel][irawfile] = _chiSquare;
      NDF_scan[_layer][_chip][_channel][irawfile] = _NDF;
      temp_scan[_layer][_chip][_channel][irawfile] = _temp;
      tempError_scan[_layer][_chip][_channel][irawfile] = _tempError;
      fitting_scan_n0[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(0);
      fitting_scan_n1[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(3);
      fitting_scan_n2[_layer][_chip][_channel][irawfile] = _fitting->GetParameter(6);
      ctap_scan[_layer][_chip][_channel][irawfile] = _ctap;
      eqf_scan[_layer][_chip][_channel][irawfile] = _eqf;
      fit_entry_scan[_layer][_chip][_channel][irawfile] = _fit_entry;
    }

    _fitting = nullptr;
  }

  cout << "best fitting finding" << endl;

  // best fit の探索
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        if(j==5&&k>29) continue;
        for(int irawfile=0; irawfile<rawfilenum; irawfile++){
          double chi2_ndf = chi2_ndf_scan[i][j][k][irawfile]>1 ? chi2_ndf_scan[i][j][k][irawfile] : 1.0/chi2_ndf_scan[i][j][k][irawfile];
          if(fit_entry_scan[i][j][k][irawfile] < 10) continue;
          double n0 =  fitting_scan_n0[i][j][k][irawfile];
          double n1 =  fitting_scan_n1[i][j][k][irawfile];
          double n2 =  fitting_scan_n2[i][j][k][irawfile];
          if(!best_fitting(n0,n1,n2)) continue;
          if(max_chi2_ndf_dac[i][j][k]==-1){          /////////////////////////// ここをchi2殻だけじゃなくて、tripleになってるかを判断するべき
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac[i][j][k] = irawfile;
            // cout << max_chi2_ndf[i][j][k] << " " << chi2_ndf << " " << chi2_ndf_scan[i][j][k][irawfile] << "  " << chi2_ndf_scan[i][j][k][irawfile] << "  " << 1.0/chi2_ndf_scan[i][j][k][irawfile] << endl;
            continue;
          }
          bool misfit = (layerNu<4 || layerNu>27) ? (gain_scan[i][j][k][irawfile] == 25) : (gain_scan[i][j][k][irawfile] == 15);
          if(misfit) continue;

          if(chi2_ndf < max_chi2_ndf[i][j][k]){
            max_chi2_ndf[i][j][k] = chi2_ndf;
            max_chi2_ndf_dac[i][j][k] = irawfile;
          }
        }
      }
    }
  }

  cout << "best fitting finding finished" << endl;

  int coutcheck = 0;

  cout << coutcheck++ << endl;

  // treeに詰めている
  fileout.cd();
  fileout.mkdir("chip_canvas");
  fileout.cd("chip_canvas");
  TCanvas *C_gain_high[layerNu][chipNu];
  for(int i=0; i<layerNu;i++){
    coutcheck = 0;
    for(int j=0; j<chipNu;j++){
      C_gain_high[i][j] = new TCanvas(Form("layer%d_chip%d_gain", i, j), Form("layer%d_chip%d_gain", i, j), 1);
      C_gain_high[i][j]->Divide(6, 6);

      for(int k=0; k<channelNu;k++){
        if(j==5 && k>29) continue;
cout << coutcheck++ << endl;
        int max_file = max_chi2_ndf_dac[i][j][k];
cout << coutcheck++ << endl;
        if(max_chi2_ndf_dac[i][j][k] == -1) continue;
cout << coutcheck++ << endl;
        

        // if(( ((_layer>3&&_layer<28) && (_gain>=gain_min_10 && _gain<=gain_max_10)) || (!(_layer>3&&_layer<28) && (_gain>=gain_min_15 && _gain<=gain_max_15)) ) && _gainError>=0.01) {
        //   gain_temp[_layer][_chip][_channel]->SetPoint(irawfile, _temp, _gain);
        //   gain_temp[_layer][_chip][_channel]->SetPointError(irawfile, _tempError, _gainError);
        // }
        // if(_chiSquare/_NDF <1.5){
        //   gain_temp_cut[_layer][_chip][_channel]->SetPoint(irawfile, _temp, _gain);
        //   gain_temp_cut[_layer][_chip][_channel]->SetPointError(irawfile, _tempError, _gainError);
        // }
        // cout << i << "," << j << ","  << k << endl;
        // if(i*1e6 + j*1e4 + k > 17010003) continue;

        cellIDs = i*1e6 + j*1e4 + k;
// cout << coutcheck++ << endl;
        // runNo = irawfile;
        dac = max_file*50 + 2100;
cout << coutcheck++ << endl;
        gain = gain_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        gainError = gainError_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        chiSquare = chiSquare_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        NDF = NDF_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        temp = temp_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        tempError = tempError_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        // fitting = fitting_scan[i][j][k][max_file];
        ctap = ctap_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        eqf = eqf_scan[i][j][k][max_file];
cout << coutcheck++ << endl;
        fit_entry = fit_entry_scan[i][j][k][max_file];
cout << coutcheck++ << endl;

        gain_tree->Fill();


        
        spSpectrum[i][j][k] = (TH1F*)filein[max_file]->Get(Form("raw_histogram/time_%d_%d_%d",i,j,k));
cout << coutcheck++ << endl;
        // if(spSpectrum[i][j][k] == nullptr) continue;
        int entry = spSpectrum[i][j][k]->Integral(300,600);
        // cout << entry << ", " << spSpectrum[i][j][k] << endl;
        C_gain_high[i][j]->cd(k + 1);
        // spSpectrum[i][j][k]->GetXaxis()->SetRange(300,600);
        spSpectrum[i][j][k]->GetXaxis()->SetRange(150,300);
        spSpectrum[i][j][k]->Draw();
cout << coutcheck++ << endl;
      }
      C_gain_high[i][j]->Write();
cout << coutcheck++ << endl;
    }
  }
  fileout.cd();
  gain_tree->Write();







}