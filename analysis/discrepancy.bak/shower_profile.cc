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
#include <set>
#include <map>
#include <tuple>
#include "TStyle.h"
// #include "langaus_newConvoluted.C"
// #include "EBUdecode.h"
// #include "EBUdecode.cxx"
// #include "fft_fullrange_tempcorr_plusNphoto_v5.C"

using namespace std;

// ssaしたfileからshower profileなどを求めるマクロ
// 


const bool fitting = true;

double fixTemp = 25;



int main(int argc, char* argv[])
{ 
  if(argc < 3 && argc%2==1){                                                     //エラー処理
    cout << "usage: ./shower_profile output.root input1_ssa.root input1_decode.root  ..." << endl;
    return 1;
  }


  cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  int total_entry_max=0;

  for(int i=0; i<rawfilenum; i++){
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("SSA_Hit");
    // entry_max[i] = tree[i]->GetEntries();
  }


  int cycleID, triggerID;
  vector<int> *cellIDs = nullptr;
  vector<int> *BCIDs = nullptr;
  vector<int> *hitTags = nullptr;
  vector<int> *gainTags = nullptr;
  vector<double> *charges = nullptr;
  vector<double> *times = nullptr;
  vector<vector<double>> *temp = nullptr;
  // vector<double> *posX = nullptr;
  // vector<double> *posY = nullptr;
  // vector<double> *posZ = nullptr;

  int cycleID_ssa, triggerID_ssa;
  vector<int> *cellIDs_ssa = nullptr;
  vector<int> *BCIDs_ssa = nullptr;
  vector<int> *hitTags_ssa = nullptr;
  vector<int> *gainTags_ssa = nullptr;
  vector<double> *charges_ssa = nullptr;
  vector<double> *times_ssa = nullptr;
  vector<double> *temp_ssa = nullptr;
  vector<double> *posX_ssa = nullptr;
  vector<double> *posY_ssa = nullptr;
  vector<double> *posZ_ssa = nullptr;
  vector<int> *ssaTag_ssa = nullptr;
  


  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;



  /*
  // pedestal
  TFile *fileinPed;
  fileinPed = new TFile("/megraid01/users/data_beamtest/ECAL_data/backup/analysed_oldFormat/2023/ps/ped/fit/pedestal/ECAL_Run23_20230518_153022.root");
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
  */

  /*
    // gain w/ temperature correction
    TFile *fileinGain;
    fileinGain = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/ECAL/results/led/gain_all.root");
    TTree *treeGain = (TTree*) fileinGain->Get("temperatureDep");
    int entry_max_Gain = treeGain->GetEntries();
    int cellIDs_gain = 0;
    double slope = 0, intercept = 0, slope_error = 0, intercept_error = 0;
    double gain_slope[layerNu][chipNu][channelNu];
    double gain_intercept[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          gain_slope[i_layer][i_chip][i_channel] = 0;
          gain_intercept[i_layer][i_chip][i_channel] = 0;
        }
      }
    }
    treeGain->SetBranchAddress("cellIDs", &cellIDs_gain);
    treeGain->SetBranchAddress("slope", &slope);
    treeGain->SetBranchAddress("intercept", &intercept);
    treeGain->SetBranchAddress("slope_error", &slope_error);
    treeGain->SetBranchAddress("intercept_error", &intercept_error);
    for(int ientry=0; ientry<entry_max_Gain; ientry++){
      treeGain->GetEntry(ientry);
      // for(int i=0; i<(int)cellIDs_gain->size(); i++){
        int _layer = cellIDs_gain/1000000;
        int _chip = (cellIDs_gain/10000)%100;
        int _channel = cellIDs_gain%100;

        gain_slope[_layer][_chip][_channel] = slope;
        gain_intercept[_layer][_chip][_channel] = intercept;
      // }
    }
  fileinGain->Close();
  */

  /*
  // gain w/ temperature correction
    TFile *fileinGain;
    fileinGain = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/ECAL/results/led/20230518_1533_Calib_LED_all.root");
    TTree *treeGain = (TTree*) fileinGain->Get("T_Event");
    int entry_max_Gain = treeGain->GetEntries();
    int cellIDs_gain = 0;
    double gain_led = 0;
    double led_gain[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          led_gain[i_layer][i_chip][i_channel] = 0;
        }
      }
    }
    treeGain->SetBranchAddress("cellIDs", &cellIDs_gain);
    treeGain->SetBranchAddress("gain", &gain_led);
    for(int ientry=0; ientry<entry_max_Gain; ientry++){
      treeGain->GetEntry(ientry);
      // for(int i=0; i<(int)cellIDs_gain->size(); i++){
        int _layer = cellIDs_gain/1e6;
        int _chip = (cellIDs_gain/10000)%100;
        int _channel = cellIDs_gain%100;

        led_gain[_layer][_chip][_channel] = gain_led;
      // }
    }
  fileinGain->Close();
  
  // fft gain
    TFile *fileinGainFFT;
    fileinGainFFT = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/fft_gain.root");
    TTree *treeGainFFT = (TTree*) fileinGainFFT->Get("T_Event");
    int entry_max_GainFFT = treeGainFFT->GetEntries();
    int cellIDs_gainfft = 0;
    double gain_fft = 0;
    double fft_gain[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          fft_gain[i_layer][i_chip][i_channel] = 0;
        }
      }
    }
    treeGainFFT->SetBranchAddress("CellID", &cellIDs_gainfft);
    treeGainFFT->SetBranchAddress("gain_ssaCut", &gain_fft);
    for(int ientry=0; ientry<entry_max_GainFFT; ientry++){
      treeGainFFT->GetEntry(ientry);
      // for(int i=0; i<(int)cellIDs_gain->size(); i++){
        int _layer = cellIDs_gainfft/1e5;
        int _chip = cellIDs_gainfft%100000/1e4;
        int _channel = cellIDs_gainfft%100;

        fft_gain[_layer][_chip][_channel] = gain_fft;
        cout << gain_fft << endl;
      // }
    }
  fileinGainFFT->Close();
  */





  TFile fileout(argv[1],"RECREATE");
  fileout.mkdir("hitmap_canvas");
  // gain_histo->Write();
  


  /*
  // 保存形式
  double _landauMPV, _landauWidth, _gausSigma, _landauMPV_tempCorrected, _landauWidth_tempCorrected, _gausSigma_tempCorrected, _ChiSqr, _ChiSqr_tempCorrected;
  int _chn_Entries;
  int _cellID_fit;
  int _Ndf, _Ndf_tempCorrected;  
  double _gain, _gain_ssaCut, _gain_tempCorrected, _gain_led;
  TH1D *_MIP_histo = new TH1D();
  TH1D *_MIP_histo_quantized = new TH1D();
  TH1D *_MIP_histo_ssaCut = new TH1D();
  TH1D *_MIP_histo_tempCorrected = new TH1D();
  TGraphErrors *_temp_graph = new TGraphErrors();

  TTree* mip_tree = new TTree("T_Event", "after calibration");
  mip_tree->Branch("chnEntries", &_chn_Entries);
  mip_tree->Branch("CellID", &_cellID_fit);
  mip_tree->Branch("gain_led", &_gain_led);
  mip_tree->Branch("gain", &_gain);
  mip_tree->Branch("gain_ssaCut", &_gain_ssaCut);
  mip_tree->Branch("gain_tempCorrected", &_gain_tempCorrected);
  mip_tree->Branch("landauWidth", &_landauWidth);
  mip_tree->Branch("gausSigma", &_gausSigma);
  mip_tree->Branch("ChiSqr", &_ChiSqr);
  mip_tree->Branch("ndf", &_Ndf);
  mip_tree->Branch("histo", &_MIP_histo);
  mip_tree->Branch("histo_quantized", &_MIP_histo_quantized);
  mip_tree->Branch("histo_ssaCut", &_MIP_histo_ssaCut);
  mip_tree->Branch("histo_tempCorrected", &_MIP_histo_tempCorrected);
  mip_tree->Branch("landauWidth_tempCorrected", &_landauWidth_tempCorrected);
  mip_tree->Branch("gausSigma_tempCorrected", &_gausSigma_tempCorrected);
  mip_tree->Branch("landauMPV_tempCorrected", &_landauMPV_tempCorrected);
  mip_tree->Branch("ChiSqr_tempCorrected", &_ChiSqr_tempCorrected);
  mip_tree->Branch("ndf_tempCorrected", &_Ndf_tempCorrected);
  mip_tree->Branch("histo_tempCorrected", &_MIP_histo_tempCorrected);
  mip_tree->Branch("tempGraph", &_temp_graph);


  Double_t langaus_peakP, langaus_peakPError;
  TF1 *langaus_func;
  */



  cout << "made save TTree" << endl;


  // int bin_10=210, max_10=1000, min_10=-50;
  // int bin_15=620, max_15=3000, min_15=-100;
  int bin_10=1050, max_10=1000, min_10=-50;
  int binw_15 = 4;
  int bin_15=3100/binw_15, max_15=3000, min_15=-100;
  int bin_q=124, max_q=120, min_q=-4;


  char char_tmp[100];
  TH2D* hit_map[layerNu];
  // TH1D* h_MIP[layerNu][chipNu][channelNu];
  // TH1D* h_MIP_quantize[layerNu][chipNu][channelNu];
  // TH1D* h_MIP_ssaCut[layerNu][chipNu][channelNu];
  // TH1D* h_MIP_tempCorrected[layerNu][chipNu][channelNu];
  // TGraphErrors* g_temperature[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    sprintf(char_tmp, "Layer%d", i_layer);
    hit_map[i_layer] = new TH2D(char_tmp, char_tmp, 42,-111.3,111.3, 42,-111.3,111.3);
  }




  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;
    // vector<vector<int>> hitStrip;
    if(!tree[irawfile]){
      cout << "   error in root file conversion!!!!!" << endl;
      cout << "      " << argv[irawfile+2] << endl;
      continue;
    }
    entry_max[irawfile] = tree[irawfile]->GetEntries();
  
    // ssa
    tree[irawfile]->SetBranchAddress("Event_Num", &triggerID_ssa);
    // tree[irawfile]->SetBranchAddress("BCID", &BCIDs_ssa);
    tree[irawfile]->SetBranchAddress("CellID", &cellIDs_ssa);
    // tree[irawfile]->SetBranchAddress("LG_Charge", &charges_ssa);
    // tree[irawfile]->SetBranchAddress("HG_Charge", &times_ssa);
    // tree[irawfile]->SetBranchAddress("HitTag", &hitTags_ssa);
    // tree[irawfile]->SetBranchAddress("ssaTag", &gainTags_ssa);
    tree[irawfile]->SetBranchAddress("ssaTag", &ssaTag_ssa);
    // tree[irawfile]->SetBranchAddress("temp", &temp);
    tree[irawfile]->SetBranchAddress("Hit_X", &posX_ssa);
    tree[irawfile]->SetBranchAddress("Hit_Y", &posY_ssa);
    tree[irawfile]->SetBranchAddress("Hit_Z", &posZ_ssa);

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      if(ientry%10000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
      tree[irawfile]->GetEntry(ientry);
      // vector<int> v;
      // int layerHit[32] = {0};
      // int nlayerHit = 0;
      set<tuple<int, int, int>> unique_hits;

      if(cellIDs_ssa->size()==0) continue;
      // v.push_back((int)BCIDs_ssa->at(0));
      // v.push_back((int)triggerID_ssa);
      for(int i=0; i<(int)cellIDs_ssa->size(); i++){
        // cout << "   " << i << "/" << cellIDs_ssa->size() << endl;
        int _layer = cellIDs_ssa->at(i)/1e5;
        int _chip = (cellIDs_ssa->at(i)/10000)%10;
        int _channel = cellIDs_ssa->at(i)%100;
        // double _charges = charges_ssa->at(i);
        // double _times = times_ssa->at(i);
        // int _hitTags = hitTags_ssa->at(i);
        double ix = posX_ssa->at(i);
        double iy = posY_ssa->at(i);
        double iz = posZ_ssa->at(i);
        int _ssaTag = ssaTag_ssa->at(i);
        if(_ssaTag==0) continue;

        auto coord_key = std::make_tuple(ix, iy, iz);
        if (unique_hits.count(coord_key) > 0) continue;
        unique_hits.insert(coord_key); // 新しい座標の組み合わせを記録

        hit_map[_layer]->Fill(posX_ssa->at(i), posY_ssa->at(i));

        // auto result = find(v.begin()+1, v.end(), (int)cellIDs_ssa->at(i));
        // if (result == v.end()) {
        //   layerHit[_layer]++;
        //   v.push_back((int)cellIDs_ssa->at(i));
        // }

      }
      // for(int i=0;i<32;i++) if(layerHit[i]>0) nlayerHit++;
      // if(v.size()>=15 && v.size()<64 && nlayerHit>=15) hitStrip.push_back(v);
    }
    triggerID_ssa = 0;
    BCIDs_ssa = nullptr;
    cellIDs_ssa = nullptr;
    charges_ssa = nullptr;
    times_ssa = nullptr;
    hitTags_ssa = nullptr;
    gainTags_ssa = nullptr;
    ssaTag_ssa = nullptr;
    posX_ssa = nullptr;
    posY_ssa = nullptr;
    posZ_ssa = nullptr;
    // }
    filein[irawfile]->Close();

    /*
    // decode
    // for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    irawfile++;
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

    tree[irawfile]->SetBranchAddress("TriggerID", &triggerID);
    tree[irawfile]->SetBranchAddress("BCID", &BCIDs);
    tree[irawfile]->SetBranchAddress("CellID", &cellIDs);
    tree[irawfile]->SetBranchAddress("LG_Charge", &charges);
    tree[irawfile]->SetBranchAddress("HG_Charge", &times);
    tree[irawfile]->SetBranchAddress("HitTag", &hitTags);
    tree[irawfile]->SetBranchAddress("GainTag", &gainTags);
    tree[irawfile]->SetBranchAddress("Temperature", &temp);
    // t[irawfile]->SetBranchAddress("posX", &posX);
    // t[irawfile]->SetBranchAddress("posY", &posY);
    // t[irawfile]->SetBranchAddress("posZ", &posZ);

    TH1D* h_temperature[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          sprintf(char_tmp, "temperature_Layer%d_Chip%d_Chn%d_run%d", i_layer, i_chip, i_channel, (int)irawfile/2);
          h_temperature[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 200, 10, 30);
          h_temperature[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("temperature");
        }
      }
    }

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      if(ientry%10000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
      tree[irawfile]->GetEntry(ientry);
      if(cellIDs->size()==0) continue;

      int itr = 0;
      bool itr_ = false;
      for(itr=0;itr<(int)hitStrip.size();itr++){
        if(hitStrip[itr][0]==triggerID){
          itr_ = true;
          break;
        }
      }
      // if(cut) if(itr==hitStrip.size()) continue;

      if(!itr_) continue;

      // int nchip = cellIDs->size()/channelNu;

      for(int i=0; i<(int)cellIDs->size(); i++){
        auto result = itr_ ? find(hitStrip[itr].begin()+1, hitStrip[itr].end(), (int)cellIDs->at(i)) : hitStrip[itr].end();
        // if(cut) if(result == hitStrip[itr].end()) continue;
        // cout << "   " << i << "/" << cellIDs->size() << endl;
        int _layer = cellIDs->at(i)/1e5;
        int _chip = (cellIDs->at(i)/10000)%10;
        int _channel = cellIDs->at(i)%100;
        double _charges = charges->at(i);
        double _times = times->at(i);
        int _hitTags = hitTags->at(i);
        // int _gainTags = gainTags->at(i);
        // double _x = posX->at(i);
        // double _y = posY->at(i);
        // double _z = posZ->at(i);

        if(_hitTags==0) continue;

        double *_position = EBUdecode(_layer,_chip,_channel);
        double SiPMtemp;
        if(temp->at(_layer).size()!=0) SiPMtemp = tempReconstruction(_layer, _position, temp->at(_layer));
        else SiPMtemp=fixTemp;
        // double slope_now = gain_slope[_layer][_chip][_channel];
        // double intercept_now = gain_intercept[_layer][_chip][_channel];
        // double _charge_corrected = (_charges - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];
        // double _time_corrected = (_times - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];
        // double slope_now = _layer<4 || _layer>28 ? 3.5/2.3/100 : 1.6/1.35/100;
        // double intercept_now = 0;
        // double _charge_corrected = (_charges - HG_ped[_layer][_chip][_channel]) / (1-slope_now)  + HG_ped[_layer][_chip][_channel];
        // double _time_corrected = (_times - HG_ped[_layer][_chip][_channel]) / (1-slope_now)  + HG_ped[_layer][_chip][_channel];

        _times = _times - HG_ped[_layer][_chip][_channel];
        _times = (1+(SiPMtemp-fixTemp)*0.03)*_times;


        h_temperature[_layer][_chip][_channel]->Fill(SiPMtemp);

        // if(_layer==0 && _chip==2 && _channel==11) cout << slope_now << ", " << SiPMtemp << ", " << intercept_now << ", " << _times << ", " << HG_ped[_layer][_chip][_channel] << ", " << _time_corrected << ", " << slope_now*SiPMtemp+intercept_now << "," << slope_now*20.0+intercept_now << endl;


        h_MIP[_layer][_chip][_channel]->Fill(_times);
        if(result != hitStrip[itr].end()){
          if(fft_gain[_layer][_chip][_channel]!=0) h_MIP_quantize[_layer][_chip][_channel]->Fill(_times/fft_gain[_layer][_chip][_channel]);
          h_MIP_ssaCut[_layer][_chip][_channel]->Fill(_times);
          h_MIP_tempCorrected[_layer][_chip][_channel]->Fill(_times);
        }
      }
    }
    triggerID = 0;
    BCIDs = nullptr;
    cellIDs = nullptr;
    charges = nullptr;
    times = nullptr;
    hitTags = nullptr;
    gainTags = nullptr;
    filein[irawfile]->Close();

    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          g_temperature[i_layer][i_chip][i_channel]->SetPoint(irawfile/2, irawfile/2, h_temperature[i_layer][i_chip][i_channel]->GetMean());
          g_temperature[i_layer][i_chip][i_channel]->SetPointError(irawfile/2, 0, h_temperature[i_layer][i_chip][i_channel]->GetRMS());
          h_temperature[i_layer][i_chip][i_channel]->Delete();
        }
      }
    }
    */
    fileout.cd();
  }
  
  cout << "saving histograms ... " << endl;

  fileout.cd();
  fileout.cd("hitmap_canvas");
  cout << "hitmap_canvas" << endl;
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    cout << i_layer << "/" << layerNu << endl;
    hit_map[i_layer]->Write();
  }

  // data の保存をしている
  fileout.cd();
  TCanvas* C_hitmap[2];
  sprintf(char_tmp, "layer0_15");
  C_hitmap[0] = new TCanvas(char_tmp, char_tmp, 3000,3000);
  C_hitmap[0]->Divide(4,4);
  sprintf(char_tmp, "layer16_31");
  C_hitmap[1] = new TCanvas(char_tmp, char_tmp, 3000,3000);
  C_hitmap[1]->Divide(4,4);
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    int icanvas = i_layer/16;
    C_hitmap[icanvas]->cd(i_layer%16+1);
    gStyle->SetOptStat(0);
    hit_map[i_layer]->Draw("colz");
    gPad-> SetLogz();
  }
  C_hitmap[0]->Update();
  C_hitmap[0]->Write();
  C_hitmap[0]->SaveAs(Form("../result/discrepancy/simulation/layer0_15.png"));
  C_hitmap[1]->Update();
  C_hitmap[1]->Write();
  C_hitmap[1]->SaveAs(Form("../result/discrepancy/simulation/layer16_31.png"));


  // TCanvas* C_MIP_quantize[layerNu][chipNu];
  // TCanvas* C_MIP_ssaCut[layerNu][chipNu];
  // TCanvas* C_MIP_tempCorrected[layerNu][chipNu];
  // for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
  //   cout << i_layer << "/" << layerNu << endl;
  //     sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
  //     C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
  //     C_MIP[i_layer][i_chip]->Divide(6, 6);
  //     sprintf(char_tmp, "layer%d_chip%d_quantize", i_layer, i_chip);
  //     C_MIP_quantize[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
  //     C_MIP_quantize[i_layer][i_chip]->Divide(6, 6);
  //     sprintf(char_tmp, "layer%d_chip%d_ssaCut", i_layer, i_chip);
  //     C_MIP_ssaCut[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
  //     C_MIP_ssaCut[i_layer][i_chip]->Divide(6, 6);
  //     sprintf(char_tmp, "layer%d_chip%d_tempCorrected", i_layer, i_chip);
  //     C_MIP_tempCorrected[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
  //     C_MIP_tempCorrected[i_layer][i_chip]->Divide(6, 6);

  //     for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
  //       if (i_chip == 5 && i_channel > 29) continue;
  //       _chn_Entries = h_MIP[i_layer][i_chip][i_channel]->GetEntries();
  //       _cellID_fit = i_layer * 1e5 + i_chip * 1e4 + i_channel;
  //       _gain=0,  _gain_ssaCut=0, _gain_tempCorrected=0;
  //       _gain_led = led_gain[i_layer][i_chip][i_channel];


  //       C_MIP[i_layer][i_chip]->cd(i_channel + 1);
  //       h_MIP[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
  //       h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
  //       // h_MIP[i_layer][i_chip][i_channel]->Rebin(5);
  //       h_MIP[i_layer][i_chip][i_channel]->Draw("HIST");

  //       C_MIP_quantize[i_layer][i_chip]->cd(i_channel + 1);
  //       h_MIP_quantize[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
  //       h_MIP_quantize[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
  //       h_MIP_quantize[i_layer][i_chip][i_channel]->Draw("HIST");
        
  //       C_MIP_ssaCut[i_layer][i_chip]->cd(i_channel + 1);
  //       h_MIP_ssaCut[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
  //       h_MIP_ssaCut[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
  //       // h_MIP[i_layer][i_chip][i_channel]->Scale(1./h_MIP[i_layer][i_chip][i_channel]->Integral());
  //       // h_MIP[i_layer][i_chip][i_channel]->Draw("HIST");
  //       // h_MIP_ssaCut[i_layer][i_chip][i_channel]->SetLineColor(kBlue);
  //       // h_MIP_ssaCut[i_layer][i_chip][i_channel]->Scale(1./h_MIP_ssaCut[i_layer][i_chip][i_channel]->Integral());
  //       // h_MIP_ssaCut[i_layer][i_chip][i_channel]->Draw("same");
  //       h_MIP_ssaCut[i_layer][i_chip][i_channel]->Draw("HIST");

  //       C_MIP_tempCorrected[i_layer][i_chip]->cd(i_channel + 1);
  //       h_MIP_tempCorrected[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
  //       h_MIP_tempCorrected[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
  //       // h_MIP_tempCorrected[i_layer][i_chip][i_channel]->Rebin(5);
  //       h_MIP_tempCorrected[i_layer][i_chip][i_channel]->Draw("HIST");

  //       mip_tree->Fill();
  //     }
  //     C_MIP[i_layer][i_chip]->Write();
  //     C_MIP_quantize[i_layer][i_chip]->Write();
  //     C_MIP_ssaCut[i_layer][i_chip]->Write();
  //     C_MIP_tempCorrected[i_layer][i_chip]->Write();
  //   }
  // }
  // mip_tree->Write();

}