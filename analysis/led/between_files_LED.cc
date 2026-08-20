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



int main(int argc, char* argv[])
{ 
    if(argc < 5){                                                     //エラー処理
        cout << "usage: ./between_files output.root outputDirectory input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = argc - 3;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  vector<pair<int,int>> filename;
  int total_entry_max=0;

	string picDirName(argv[2]);
	cout << picDirName << endl;

  double rawFileDate[rawfilenum];

  for(int i=0; i<rawfilenum; i++){
    string txtname(argv[i+3]);
    txtname.erase(0,txtname.find_last_of("p")+1);
    // txtname.erase(txtname.find_last_of("_"),txtname.end()-1);
    filename.push_back({stoi(txtname), i+3});
    // cout << argv[i+2] << ",  " << stoi(txtname) << endl;
  }
  // sort(filename.begin(),filename.end());
  for(int i=0; i<rawfilenum; i++){
    filein[i] = new TFile(argv[filename[i].second]);
    tree[i] = (TTree*) filein[i]->Get("T_Event");
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
  TH1F *charge_histo[layerNu][chipNu][channelNu];
  TH1F *time_histo[layerNu][chipNu][channelNu];
  TH1F *channel_temperature[layerNu][chipNu][channelNu];
  // int charge_count[layerNu][chipNu][channelNu];

  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        charge_histo[i][j][k] = new TH1F(Form("charge_%d_%d_%d",i,j,k),Form("charge_%d_%d_%d",i,j,k),4096,0,4095);
        time_histo[i][j][k] = new TH1F(Form("time_%d_%d_%d",i,j,k),Form("time_%d_%d_%d",i,j,k),4096,0,4095);

        channel_temperature[i][j][k] = new TH1F(Form("channel_temperature_%d_%d_%d",i,j,k),Form("channel_temperature_%d_%d_%d",i,j,k),200,15,35);
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
  double _gain, _gainError, _gausSigma, _chiSquare, _NDF, _temp, _tempError, _ctap, _eqf;
  int _cellID_fit, _fit_entry;
  TF1 *_fitting;
  TTree* gain_tree = new TTree("T_Event", "gain calibration");
  gain_tree->Branch("cellIDs", &_cellID_fit);
  gain_tree->Branch("gain", &_gain);
  gain_tree->Branch("gainError", &_gainError);
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


  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

    if(entry_max[irawfile]==0) continue;


    tree[irawfile]->SetBranchAddress("cellIDs", &cellIDs);
    tree[irawfile]->SetBranchAddress("charges", &charges);
    tree[irawfile]->SetBranchAddress("times", &times);
    tree[irawfile]->SetBranchAddress("hitTags", &hitTags);
    tree[irawfile]->SetBranchAddress("temp", &temp);

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);

      int groupNum = filename[irawfile].first;

      for(int i=0; i<(int)cellIDs->size(); i++){
        // cout << "   " << i << "/" << cellIDs->size() << endl;
        if(hitTags->at(i)==0) continue;
        int _layer = cellIDs->at(i)/1e6;
        int _chip = (cellIDs->at(i)/10000)%100;
        int _channel = cellIDs->at(i)%100;

        int index = _chip * 36 + _channel;
        if(find(group_channel[groupNum], group_channel[groupNum] + 15, index) == group_channel[groupNum] + 15) continue;

        double _charges = charges->at(i);
        double _times = times->at(i);
        // int _hitTags = hitTags->at(i);
        // double _x = posX->at(i);
        // double _y = posY->at(i);
        // double _z = posZ->at(i);

        charge_histo[_layer][_chip][_channel]->Fill(_charges);
        time_histo[_layer][_chip][_channel]->Fill(_times);
        // charge_count[_layer][_chip][_channel]++;

        double *_position = EBUdecode(_layer,_chip,_channel);
        // cout << _position[0] << ", " << _position[1] << ", " << _position[2] << endl;
        double SiPMtemp = 20;
        if(temp->at(_layer).size()!=0) {
          SiPMtemp = tempReconstruction(_layer, _position, temp->at(_layer));
          // cout << SiPMtemp << endl;
        }
        else SiPMtemp=20;
        channel_temperature[_layer][_chip][_channel]->Fill(SiPMtemp);
      }

      for(int i=0;i<temp->size();i++){
        for(int j=0;j<temp->at(i).size();j++){
          sensor_temperature[i][j]->Fill(temp->at(i).at(j));
        }
      }
    }


    cellIDs = nullptr;
    charges = nullptr;
    times = nullptr;
    hitTags = nullptr;
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


  TF1 *highfun = new TF1("highfun","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1]-[4])/[5])^2) + [6]*exp(-0.5*((x-[1]-[4]*2)/[7])^2)",300,600);
  highfun->SetNpx(4000);
  // TF1 *lowfun = new TF1("lowfun","gaus(0)+gaus(3)+gaus(6)",300,600);
  highfun->SetParLimits(4,6,34);
  highfun->SetParLimits(2,0,fit_sigma_maximum);
  highfun->SetParLimits(5,0,fit_sigma_maximum);
  highfun->SetParLimits(7,0,fit_sigma_maximum);

  fileout.mkdir("chip_canvas");
  fileout.cd("chip_canvas");
  TCanvas *C_MIP_high[layerNu][chipNu];
  TCanvas *C_MIP_low[layerNu][chipNu];
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    for(int i_chip=0;i_chip<chipNu;i_chip++){
      // sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_MIP_high[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_high", i_layer, i_chip), Form("layer%d_chip%d_high", i_layer, i_chip), 1);
      C_MIP_high[i_layer][i_chip]->Divide(6, 6);
      C_MIP_low[i_layer][i_chip] = new TCanvas(Form("layer%d_chip%d_low", i_layer, i_chip), Form("layer%d_chip%d_low", i_layer, i_chip), 1);
      C_MIP_low[i_layer][i_chip]->Divide(6, 6);
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if(i_chip==5&&i_channel>29) continue;
        C_MIP_high[i_layer][i_chip]->cd(i_channel + 1);
        int fit_entry = time_histo[i_layer][i_chip][i_channel]->Integral(300,600);
        time_histo[i_layer][i_chip][i_channel]->GetXaxis()->SetRange(300,600);
        time_histo[i_layer][i_chip][i_channel]->Rebin(2);
        time_histo[i_layer][i_chip][i_channel]->Draw();
        gPad-> SetLogy();
        // int bin0pe = time_histo[i_layer][i_chip][i_channel]->GetMaximumBin();
        // int amp0pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe);
        // int amp1pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+25);
        // int amp2pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+50);
        int init_gain = i_layer<4 || i_layer>27 ? 25 : 15;
        int bin0pe = HG_ped[i_layer][i_chip][i_channel];
        int amp0pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe);
        int amp1pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+init_gain);
        int amp2pe = time_histo[i_layer][i_chip][i_channel]->GetBinContent(bin0pe+init_gain);
        int maximum_entry = time_histo[i_layer][i_chip][i_channel]->GetBinContent(time_histo[i_layer][i_chip][i_channel]->GetMaximumBin());
        // highfun->SetParLimits(1,bin0pe-10,bin0pe+10);
        highfun->SetParLimits(1,bin0pe-20,bin0pe+20);
        highfun->SetParLimits(0,0,maximum_entry*1.5);
        highfun->SetParLimits(3,0,maximum_entry*1.5);
        highfun->SetParLimits(6,0,maximum_entry*1.5);
        if(i_layer<4 || i_layer>27) highfun->SetParameters(amp0pe,bin0pe,4,amp1pe,25,4,amp2pe,4);
        else highfun->SetParameters(amp0pe,bin0pe,4,amp1pe,15,4,amp2pe,4);
        time_histo[i_layer][i_chip][i_channel]->Fit("highfun","QR","",bin0pe-init_gain,bin0pe+4*init_gain);
        // double fit_range_low = calc_intercept1(highfun->GetParameter(0), highfun->GetParameter(1), highfun->GetParameter(2), -1);
        // double fit_range_high = calc_intercept1(highfun->GetParameter(6), highfun->GetParameter(1)+2*highfun->GetParameter(4), highfun->GetParameter(7), 1);
        double fit_range_low = highfun->GetParameter(0) - highfun->GetParameter(2) * 2;
        double fit_range_high = highfun->GetParameter(1)+2*highfun->GetParameter(4) + highfun->GetParameter(7) * 2;
        time_histo[i_layer][i_chip][i_channel]->Fit("highfun","QR","",fit_range_low,fit_range_high);
        if(i_layer<30) gain->Fill(highfun->GetParameter(4));
        if(i_layer<4 || (i_layer>27 && i_layer<30)) {
          gain_15->Fill(highfun->GetParameter(4));
          gain_15_ratio->Fill(highfun->GetParameter(4) / SPSgain[i_layer][i_chip][i_channel]);
        }
        if(i_layer>=4 && i_layer<=27) {
          gain_10->Fill(highfun->GetParameter(4));
          gain_10_ratio->Fill(highfun->GetParameter(4) / SPSgain[i_layer][i_chip][i_channel]);
        }
        gain_layer->Fill(highfun->GetParameter(4), i_layer);
        sigma_layer->Fill(highfun->GetParameter(2), i_layer);
        sigma_layer->Fill(highfun->GetParameter(5), i_layer);
        sigma_layer->Fill(highfun->GetParameter(7), i_layer);

        double p0 = highfun->GetParameter(0) /2.50662827463/ highfun->GetParameter(2) / (double)time_histo[i_layer][i_chip][i_channel]->GetEntries();
        double p1 = highfun->GetParameter(3) /2.50662827463/ highfun->GetParameter(5) / (double)time_histo[i_layer][i_chip][i_channel]->GetEntries();
        double mean = time_histo[i_layer][i_chip][i_channel]->GetMean();

        _cellID_fit = i_layer*1000000 + i_chip*10000 + i_channel;
        _gain = highfun->GetParameter(4);
        _gainError = highfun->GetParError(4);
        _chiSquare = highfun->GetChisquare();
        _NDF = highfun->GetNDF();
        _temp = channel_temperature_mean[i_layer][i_chip][i_channel];
        _tempError = channel_temperature_sigma[i_layer][i_chip][i_channel];
        _fitting = highfun;
        _ctap = ctap_calc(p0,p1);
        _eqf = ctap_eqf(p0, highfun->GetParameter(4), mean);
        _fit_entry = fit_entry;
        gain_tree->Fill();
     
        C_MIP_low[i_layer][i_chip]->cd(i_channel + 1);
        gStyle->SetOptFit(1111);
        charge_histo[i_layer][i_chip][i_channel]->GetXaxis()->SetRange(300,600);
        charge_histo[i_layer][i_chip][i_channel]->Draw();
        gPad-> SetLogy();
      }
      C_MIP_high[i_layer][i_chip]->Write();
      C_MIP_low[i_layer][i_chip]->Write();
      // C_MIP_high[i_layer][i_chip]->SaveAs(Form("%s/high_Layer%d_Chip%d.pdf",picDirName.c_str(),i_layer,i_chip));
      // C_MIP_low[i_layer][i_chip]->SaveAs(Form("%s/Layer%d_Chip%d.pdf",picDirName.c_str(),i_layer,i_chip));
    }
  }

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
  TH2F *temp_sensor = new TH2F("temp_sensor","temperature;layer;sensor;temperature [C]",32,0,32,16,0,16);
  for(int i_layer=0;i_layer<layerNu;i_layer++){
    C_temperature[i_layer] = new TCanvas(Form("layer%d_temperature", i_layer), Form("layer%d_temperature", i_layer), 1);
    C_temperature[i_layer]->Divide(4, 4);
    for(int i_sensor=0;i_sensor<16;i_sensor++){
      C_temperature[i_layer]->cd(i_sensor + 1);
      sensor_temperature[i_layer][i_sensor]->Draw();
      sensor_temperature[i_layer][i_sensor]->Fit("temp_fun","QR","",20,30);
      gPad-> SetLogy();

      temp_sensor->SetBinContent(i_layer+1,i_sensor+1, temp_fun->GetParameter(1));
    }
    C_temperature[i_layer]->Write();
  }
  temp_sensor->Write();

  fileout.cd();
  gain_tree->Write();



  fileout.mkdir("raw_histogram");
  fileout.cd("raw_histogram");
  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        charge_histo[i][j][k]->Write();
        time_histo[i][j][k]->Write();
        // charge_count[i][j][k] = 0;
        // cout << channel_temperature[i][j][k] << ", " << time_histo[i][j][k]->Integral() << endl;
        // channel_temperature_mean[i][j][k] = channel_temperature[i][j][k]->GetMean();
        // channel_temperature_sigma[i][j][k] = channel_temperature[i][j][k]->GetRMS();
      }
    }
  }
  for(int i=0;i<32;i++){
    for(int j=0;j<16;j++){
      sensor_temperature[i][j]->Write();
    }
  }
  fileout.cd();
  cout << "saved raw histograms" << endl;
  



}