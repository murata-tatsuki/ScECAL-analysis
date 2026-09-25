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
#include "TStyle.h"
#include "TPad.h"
// #include "langaus_newConvoluted.C"
#include "FitMIP.C"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;

// ssa後にcutをかけてMPVを求めるマクロ
// pedestal を引いている

const bool cut = true;  // cutをするかどうかのbool

const bool threshold_fitting = true;
const bool fitting = true;



int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./MIP_calibration output.root input1.root input2.root  ..." << endl;
        return 1;
    }

  cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  int total_entry_max=0;
	
  for(int i=0; i<rawfilenum; i++){
    // cout << argv[i+2] << endl;
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("SSA_Hit");
    entry_max[i] = tree[i]->GetEntries();
    
    i++;
    // cout << argv[i+2] << endl;
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("Raw_Hit");
    entry_max[i] = tree[i]->GetEntries();
  }


  int cycleID, triggerID;
  vector<int> *cellIDs = nullptr;
  vector<int> *BCIDs = nullptr;
  vector<int> *hitTags = nullptr;
  vector<int> *gainTags = nullptr;
  vector<double> *charges = nullptr;
  vector<double> *times = nullptr;
  // vector<vector<double>> *temp = nullptr;
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


  // pedestal
  TFile *fileinPed;
  fileinPed = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/pedestal2023_SPS.root");
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
    // treePed->SetBranchAddress("cellIDs", &cellIDs_ped);
    // treePed->SetBranchAddress("pedMeanTimes", &pedMeanTimes);
    // treePed->SetBranchAddress("pedSigmaTimes", &pedSigmaTimes);
    // treePed->SetBranchAddress("pedMeanCharges", &pedMeanCharges);
    // treePed->SetBranchAddress("pedSigmaCharges", &pedSigmaCharges);
    treePed->SetBranchAddress("CellID", &cellIDs_ped);
    treePed->SetBranchAddress("PedHighMean", &pedMeanTimes);
    treePed->SetBranchAddress("PedHighSig", &pedSigmaTimes);
    treePed->SetBranchAddress("PedLowMean", &pedMeanCharges);
    treePed->SetBranchAddress("PedLowSig", &pedSigmaCharges);
    for(int ientry=0; ientry<entry_max_Ped; ientry++){
      treePed->GetEntry(ientry);
      for(int i=0; i<(int)cellIDs_ped->size(); i++){
        int _layer = cellIDs_ped->at(i)/1e5;
        int _chip = (cellIDs_ped->at(i)%100000)/1e4;
        int _channel = cellIDs_ped->at(i)%100;
        // int _layer = cellIDs_ped->at(i)/100000;
        // int _chip = (cellIDs_ped->at(i)/10000)%100;
        // int _channel = cellIDs_ped->at(i)%100;

        HG_ped[_layer][_chip][_channel] = pedMeanTimes->at(i);
        HG_pedSigma[_layer][_chip][_channel] = pedSigmaTimes->at(i);
        LG_ped[_layer][_chip][_channel] = pedMeanCharges->at(i);
        LG_pedSigma[_layer][_chip][_channel] = pedSigmaCharges->at(i);
      }
    }
  fileinPed->Close();

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


  cout<<"=========== Threshold Extract ============"<<endl;
    TFile *ThresholdFile = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root","READ");
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





  TFile fileout(argv[1],"RECREATE");
  fileout.mkdir("chip_canvas");
  // fileout.mkdir("chip_canvas_tempCorrected");
  // fileout.mkdir("chip_canvas_run");
  // fileout.mkdir("chip_canvas_tempCorrected_run");
  // gain_histo->Write();
  



  // 保存形式
  double _landauMPV, _landauWidth, _gausSigma, _landauMPV_error, _landauWidth_error, _gausSigma_error, _landauMPV_tempCorrected, _landauWidth_tempCorrected, _gausSigma_tempCorrected, _ChiSqr, _ChiSqr_tempCorrected;
  int _chn_Entries;
  int _cellID_fit;
  int _Ndf, _Ndf_tempCorrected;  
  TH1F *_MIP_histo = new TH1F();
  TH1F *_MIP_histo_tempCorrected = new TH1F();
  TGraphErrors *_temp_graph = new TGraphErrors();

  TTree* mip_tree = new TTree("MIP_Fit", "after calibration");
  mip_tree->Branch("chnEntries", &_chn_Entries);
  mip_tree->Branch("CellID", &_cellID_fit);
  mip_tree->Branch("LandauMPV", &_landauMPV);
  mip_tree->Branch("LandauWidth", &_landauWidth);
  mip_tree->Branch("GauSigma", &_gausSigma);
  mip_tree->Branch("LandauMPV_error", &_landauMPV_error);
  mip_tree->Branch("LandauWidth_error", &_landauWidth_error);
  mip_tree->Branch("GauSigma_error", &_gausSigma_error);
  mip_tree->Branch("ChiSquare", &_ChiSqr);
  mip_tree->Branch("NDF", &_Ndf);
  // mip_tree->Branch("histo", &_MIP_histo);
  // mip_tree->Branch("landauWidth_tempCorrected", &_landauWidth_tempCorrected);
  // mip_tree->Branch("gausSigma_tempCorrected", &_gausSigma_tempCorrected);
  // mip_tree->Branch("landauMPV_tempCorrected", &_landauMPV_tempCorrected);
  // mip_tree->Branch("ChiSqr_tempCorrected", &_ChiSqr_tempCorrected);
  // mip_tree->Branch("ndf_tempCorrected", &_Ndf_tempCorrected);
  // mip_tree->Branch("histo_tempCorrected", &_MIP_histo_tempCorrected);
  // mip_tree->Branch("tempGraph", &_temp_graph);


  Double_t langaus_peakP, langaus_peakPError;
  TF1 *langaus_func;



  cout << "made save TTree" << endl;


  char char_tmp[100];
  TH1F* h_MIP[layerNu][chipNu][channelNu];
  TH1F* h_MIP_tempCorrected[layerNu][chipNu][channelNu];
  TGraphErrors* g_temperature[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
        if (3 < i_layer && i_layer < 28)
          // h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -50, 1000);
          h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 210, -50, 1000);
        else
          // h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -100, 3000);
          h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 620, -100, 3000);
        h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
        h_MIP[i_layer][i_chip][i_channel]->SetLineWidth(2);
        h_MIP[i_layer][i_chip][i_channel]->SetLineColor(4);

        sprintf(char_tmp, "Layer%d_Chip%d_Chn%d_tempCorrected", i_layer, i_chip, i_channel);
        if (3 < i_layer && i_layer < 28)
          h_MIP_tempCorrected[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 210, -50, 1000);
          // h_MIP_tempCorrected[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 1050, -50, 1000);
        else
          h_MIP_tempCorrected[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 620, -100, 3000);
          // h_MIP_tempCorrected[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 3100, -100, 3000);
        h_MIP_tempCorrected[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
        h_MIP_tempCorrected[i_layer][i_chip][i_channel]->SetLineWidth(2);
        h_MIP_tempCorrected[i_layer][i_chip][i_channel]->SetLineColor(4);

        g_temperature[i_layer][i_chip][i_channel] = new TGraphErrors();
      }
    }
  }

  // runごとのnumber of hit の確認用
  // int Nhit_run[5][rawfilenum/2];    // 20011, 12,13,14,18 の5つのchannelをまずは確認する(中心のchannel)
  // for(int irawfile=0; irawfile<rawfilenum/2; irawfile++) {
  //   for(int iii=0; iii<5; iii++) {
  //     Nhit_run[iii][irawfile] = 0;
  //   }
  // }




  // data をとってきてる
  for(int irawfile=0; irawfile<rawfilenum; irawfile++){
    if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

    vector<vector<int>> hitStrip;
    // ssa
    tree[irawfile]->SetBranchAddress("Event_Num", &triggerID_ssa);
    // tree[irawfile]->SetBranchAddress("Event_Time", &BCIDs_ssa);
    tree[irawfile]->SetBranchAddress("CellID", &cellIDs_ssa);
    // tree[irawfile]->SetBranchAddress("charges", &charges_ssa);
    // tree[irawfile]->SetBranchAddress("times", &times_ssa);
    // tree[irawfile]->SetBranchAddress("hitTags", &hitTags_ssa);
    // tree[irawfile]->SetBranchAddress("gainTags", &gainTags_ssa);
    tree[irawfile]->SetBranchAddress("ssaTag", &ssaTag_ssa);
    // t[irawfile]->SetBranchAddress("temp", &temp);
    // t[irawfile]->SetBranchAddress("posX", &posX);
    // t[irawfile]->SetBranchAddress("posY", &posY);
    // t[irawfile]->SetBranchAddress("posZ", &posZ);
    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
      tree[irawfile]->GetEntry(ientry);
      vector<int> v;
      int layerHit[32] = {0};
      int nlayerHit = 0;

      if(cellIDs_ssa->size()==0) continue;
      // v.push_back((int)BCIDs_ssa->at(0));
      v.push_back((int)triggerID_ssa);
      for(int i=0; i<(int)cellIDs_ssa->size(); i++){
        // cout << "   " << i << "/" << cellIDs_ssa->size() << endl;
        // int _layer = cellIDs_ssa->at(i)/100000;
        // int _chip = (cellIDs_ssa->at(i)/10000)%10;
        // int _channel = cellIDs_ssa->at(i)%100;
        int _layer = cellIDs_ssa->at(i)/1e5;
        int _chip = (cellIDs_ssa->at(i)%100000)/1e4;
        int _channel = cellIDs_ssa->at(i)%100;
        // cout << _layer << " " << _chip << " " << _channel << endl;
        // double _charges = charges_ssa->at(i);
        // double _times = times_ssa->at(i);
        // int _hitTags = hitTags_ssa->at(i);
        // double _x = posX->at(i);
        // double _y = posY->at(i);
        // double _z = posZ->at(i);
        int _ssaTag = ssaTag_ssa->at(i);
        if(_ssaTag==0) continue;
        // if(_hitTags==0 || _ssaTag==0) continue;
        auto result = find(v.begin()+2, v.end(), (int)cellIDs_ssa->at(i));
        if (result == v.end()) {
          layerHit[_layer]++;
          v.push_back((int)cellIDs_ssa->at(i));
        }
      }
      for(int i=0;i<32;i++) if(layerHit[i]>0) nlayerHit++;
      if(v.size()>=22 && v.size()<64 && nlayerHit>=22) hitStrip.push_back(v);
    }
    triggerID_ssa = 0;
    BCIDs_ssa = nullptr;
    cellIDs_ssa = nullptr;
    charges_ssa = nullptr;
    times_ssa = nullptr;
    hitTags_ssa = nullptr;
    gainTags_ssa = nullptr;
    ssaTag_ssa = nullptr;
    // }
    filein[irawfile]->Close();


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
    // tree[irawfile]->SetBranchAddress("temp", &temp);
    // t[irawfile]->SetBranchAddress("posX", &posX);
    // t[irawfile]->SetBranchAddress("posY", &posY);
    // t[irawfile]->SetBranchAddress("posZ", &posZ);

    TH1F* h_temperature[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          sprintf(char_tmp, "temperature_Layer%d_Chip%d_Chn%d_run%d", i_layer, i_chip, i_channel, (int)irawfile);
          h_temperature[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 200, 10, 30);
          h_temperature[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("temperature");
        }
      }
    }

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);
      if(cellIDs->size()==0) continue;

      int itr = 0;
      for(itr=0;itr<(int)hitStrip.size();itr++){
        if(hitStrip[itr][0]==triggerID) break;
        // if(hitStrip[itr][0]==BCIDs->at(0) && hitStrip[itr][1]==triggerID) break;
      }
      if(cut) if(itr==hitStrip.size()) continue;

      // int nchip = cellIDs->size()/channelNu;

      for(int i=0; i<(int)cellIDs->size(); i++){
        auto result = find(hitStrip[itr].begin()+1, hitStrip[itr].end(), (int)cellIDs->at(i));
        if(cut) if(result == hitStrip[itr].end()) continue;
        // cout << "   " << i << "/" << cellIDs->size() << endl;
        int _layer = cellIDs->at(i)/1e5;
        int _chip = (cellIDs->at(i)%100000)/1e4;
        int _channel = cellIDs->at(i)%100;
        // cout << _layer << " " << _chip << " " << _channel << endl;
        double _charges = charges->at(i);
        double _times = times->at(i);
        int _hitTags = hitTags->at(i);
        int _gainTags = gainTags->at(i);
        // double _x = posX->at(i);
        // double _y = posY->at(i);
        // double _z = posZ->at(i);

        if(_hitTags==0) continue;

        double *_position = EBUdecode(_layer,_chip,_channel);
        // double SiPMtemp;
        // if(temp->at(_layer).size()!=0) SiPMtemp = tempReconstruction(_layer, _position, temp->at(_layer));
        // else SiPMtemp=20;
        // double slope_now = gain_slope[_layer][_chip][_channel];
        // double intercept_now = gain_intercept[_layer][_chip][_channel];
        // double _charge_corrected = (_charges - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];
        // double _time_corrected = (_times - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];
        // double slope_now = _layer<4 || _layer>28 ? 3.5/2.3/100 : 1.6/1.35/100;
        // double intercept_now = 0;
        // double _charge_corrected = (_charges - HG_ped[_layer][_chip][_channel]) / (1-slope_now)  + HG_ped[_layer][_chip][_channel];
        // double _time_corrected = (_times - HG_ped[_layer][_chip][_channel]) / (1-slope_now)  + HG_ped[_layer][_chip][_channel];


        // h_temperature[_layer][_chip][_channel]->Fill(SiPMtemp);

        // if(_layer==0 && _chip==2 && _channel==11) cout << slope_now << ", " << SiPMtemp << ", " << intercept_now << ", " << _times << ", " << HG_ped[_layer][_chip][_channel] << ", " << _time_corrected << ", " << slope_now*SiPMtemp+intercept_now << "," << slope_now*20.0+intercept_now << endl;

        // runごとのnumber of hitの確認用
        // if(cellIDs->at(i)==20011) Nhit_run[0][(int)irawfile/2]++;
        // if(cellIDs->at(i)==20012) Nhit_run[1][(int)irawfile/2]++;
        // if(cellIDs->at(i)==20013) Nhit_run[2][(int)irawfile/2]++;
        // if(cellIDs->at(i)==20014) Nhit_run[3][(int)irawfile/2]++;
        // if(cellIDs->at(i)==20018) Nhit_run[4][(int)irawfile/2]++;

        h_MIP[_layer][_chip][_channel]->Fill(_times-HG_ped[_layer][_chip][_channel]);
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

    // for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    //   for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
    //     for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
    //       g_temperature[i_layer][i_chip][i_channel]->SetPoint(irawfile/2, irawfile/2, h_temperature[i_layer][i_chip][i_channel]->GetMean());
    //       g_temperature[i_layer][i_chip][i_channel]->SetPointError(irawfile/2, 0, h_temperature[i_layer][i_chip][i_channel]->GetRMS());
    //       h_temperature[i_layer][i_chip][i_channel]->Delete();
    //     }
    //   }
    // }
    fileout.cd();
  }

  cout << "saving histograms ... " << endl;

  fileout.cd();
  fileout.cd("chip_canvas");
  cout << "chip_canvas" << endl;
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    cout << i_layer << "/" << layerNu << endl;
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        h_MIP[i_layer][i_chip][i_channel]->Write();
        // h_MIP[i_layer][i_chip][i_channel]->Delete();
      }
    }
  }
  fileout.cd();
  // fileout.cd("chip_canvas_tempCorrected");
  // cout << "chip_canvas_tempCorrected" << endl;
  // for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
  //   cout << i_layer << "/" << layerNu << endl;
  //   for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
  //     for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
  //       h_MIP_tempCorrected[i_layer][i_chip][i_channel]->Write();
  //       // h_MIP_tempCorrected[i_layer][i_chip][i_channel]->Delete();
  //     }
  //   }
  // }

  TH1F *chi2ndf[2];
  TH1F *MIPs[2];
  TH1F *MIPs_width[2];
  TH1F *MIPs_GauSigma[2];
  TH1F *MIPs_error[2];
  TH1F *MIPs_width_error[2];
  TH1F *MIPs_GauSigma_error[2];
  for(int pix=0;pix<2;pix++){
    int pixel = pix == 0 ? 10 : 15;
    MIPs[pix] = new TH1F(Form("mip_%d",pixel), Form("MIP (%d um)",pixel), 200,0,1000);
    MIPs_width[pix] = new TH1F(Form("mip_width_%d",pixel), Form("MIP landau width (%d um)",pixel), 1000,0.1,200);
    MIPs_GauSigma[pix] = new TH1F(Form("MIPs_GauSigma_%d",pixel), Form("MIP gaussian sigma (%d um)",pixel), 1000,0.1,200);

    MIPs_error[pix] = new TH1F(Form("mip_error_%d",pixel), Form("MIP error (%d um)",pixel), 1000,0,200);
    MIPs_width_error[pix] = new TH1F(Form("mip_width_error_%d",pixel), Form("MIP landau width error (%d um)",pixel), 201,0.1,20);
    MIPs_GauSigma_error[pix] = new TH1F(Form("MIPs_GauSigma_error_%d",pixel), Form("MIP gaussian sigma error (%d um)",pixel), 201,0.1,20);

    chi2ndf[pix] = new TH1F(Form("chi2ndf_%d",pixel), Form("chi2/ndf (%d um)",pixel), 1000,0,10);
  }


  // data の保存をしている
  fileout.cd();
  TCanvas* C_MIP[layerNu][chipNu];
  TCanvas* C_MIP_tempCorrected[layerNu][chipNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    cout << i_layer << "/" << layerNu << endl;
    int pix = (i_layer>3 && i_layer<28) ? 0 : 1;
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 2560, 1440);
      C_MIP[i_layer][i_chip]->Divide(6, 6);
      sprintf(char_tmp, "layer%d_chip%d_tempCorrected", i_layer, i_chip);
      C_MIP_tempCorrected[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 2560, 1440);
      C_MIP_tempCorrected[i_layer][i_chip]->Divide(6, 6);

      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if (i_chip == 5 && i_channel > 29) continue;
        _chn_Entries = h_MIP[i_layer][i_chip][i_channel]->GetEntries();
        _cellID_fit = i_layer * 1e5 + i_chip * 1e4 + i_channel;


        C_MIP[i_layer][i_chip]->cd(i_channel + 1);
        h_MIP[i_layer][i_chip][i_channel]->Draw();

        if(_chn_Entries<100) continue;

        // /////////////////////////// langaus setting 
        // double fr[2];
        // double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
        // fr[0]=i_layer<4 || i_layer>28 ? 0.3*h_MIP[i_layer][i_chip][i_channel]->GetMean() : 0;
        // fr[1]=3.0*h_MIP[i_layer][i_chip][i_channel]->GetMean();
        // pllo[0]=0.5; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
        // plhi[0]=5.0; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=5.0;
        // sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;
        // double chisqr;
        // int    ndf;
        // sv[1] = i_layer<4 || i_layer>28 ? 45 : 20;
        // sv[3] = i_layer<4 || i_layer>28 ? 100 : 30;
        // pllo[0] = 4;
        // pllo[1] = i_layer<4 || i_layer>28 ? 200 : 0;
        // pllo[3] = 0;
        // plhi[0] = 150;
        // plhi[1] = i_layer<4 || i_layer>28 ? 1000 : 400;
        // plhi[3] = 250;
        // ///////////////////////////
        // double SNRPeak, SNRFWHM;


        Double_t my_threshold    = thresholdValue[i_layer][i_chip][i_channel]; // 事前に求めた Threshold
        Double_t my_thresh_width = thresholdSigmaValue[i_layer][i_chip][i_channel];  // 事前に求めた S-curve 幅

        // チャンネル番号から MPPC の種類を判定してパラメータを切り替える例

        bool is_10um = i_layer>3 && i_layer<28;
        
        double init_mpv   = is_10um ? 100.0 : 500.0; // MPV初期値
        double fix_width  = is_10um ?  22.0 :  70.0; // 共通固定Width (プロットのMean/Peakより)
        double fit_max    = is_10um ? 800.0 : 2000.0; // フィット範囲上限

        // 切り替えたパラメータでフィットを呼び出す
        MIPFitResult res = FitMIPHistogram(h_MIP[i_layer][i_chip][i_channel], pix, my_threshold, my_thresh_width, fix_width, -1.0, init_mpv, 0.0, fit_max);
        gStyle->SetOptFit(1111);
        h_MIP[i_layer][i_chip][i_channel]->Draw();

        // --- 通常の解析（全パラメータ可変）---
        // MIPFitResult res = FitMIPHistogram(h_MIP[i_layer][i_chip][i_channel], my_threshold, my_thresh_width);

        // --- Threshold が MIP より大きく、ピークが見えないチャンネルの場合 ---
        // 他の正常チャンネルから求めた Width (例: 12.5) や Sigma (例: 18.0) を固定して呼び出す
        // MIPFitResult res = FitMIPHistogram(h_muon, my_threshold, my_thresh_width, 12.5, 18.0);

        // 抽出された値と誤差を出力
        if (res.status == 0) {
          _landauMPV = res.mpv;
          _landauWidth = res.width;
          _gausSigma = res.gsigma;
          _landauMPV_error = res.mpv_err;
          _landauWidth_error = res.width_err;
          _gausSigma_error = res.gsigma_err;
          _ChiSqr = res.chi2;
          _Ndf = res.ndf;

          if(res.mpv>0){
            MIPs[pix]->Fill(res.mpv);
            MIPs_width[pix]->Fill(res.width);
            MIPs_GauSigma[pix]->Fill(res.gsigma);
            MIPs_error[pix]->Fill(res.mpv_err);
            MIPs_width_error[pix]->Fill(res.width_err);
            MIPs_GauSigma_error[pix]->Fill(res.gsigma_err);
          }
          chi2ndf[pix]->Fill(res.chi2/res.ndf);
        } else {
          _landauMPV = 0;
          _landauWidth = 0;
          _gausSigma = 0;
          _landauMPV_error = 0;
          _landauWidth_error = 0;
          _gausSigma_error = 0;
          _ChiSqr = 0;
          _Ndf = 0;
        }

        // if(fitting){
        //   TF1 *fitsnr = langaufit(h_MIP[i_layer][i_chip][i_channel],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
        //   langaupro(fp,SNRPeak,SNRFWHM);
        //   // Global style settings
        //   gStyle->SetOptStat(1111);
        //   gStyle->SetOptFit(111);
        //   gStyle->SetLabelSize(0.03,"x");
        //   gStyle->SetLabelSize(0.03,"y");
        //   fitsnr->Draw("same");
        // } else {
        //   _landauWidth = fp[0];
        //   _landauMPV = fp[1];
        //   _gausSigma = fp[3];
        //   _ChiSqr = chisqr;
        //   _Ndf = ndf;
        //   _MIP_histo = h_MIP[i_layer][i_chip][i_channel];
        //   mip_tree->Fill();
        //   continue;

        // if(i_layer<4 || i_layer>28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+200, HG_ped[i_layer][i_chip][i_channel]+1000);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 50, 400);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 400, 800);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+50, HG_ped[i_layer][i_chip][i_channel]+400);

        // _landauWidth = langaus_peakPError;
        // _landauMPV = langaus_peakP;
        //_gausSigma = fun->GetParameter(3);
        // _ChiSqr = langaus_func->GetChisquare();
        // _Ndf = langaus_func->GetNDF();
        // _landauWidth = fp[0];
        // _landauMPV = fp[1];
        // _gausSigma = fp[3];
        // _ChiSqr = chisqr;
        // _Ndf = ndf;
        // cout<<"--------  "<<_landauMPV<<" , "<<_landauWidth<<" , "<<_gausSigma<<endl;
        // _MIP_histo = h_MIP[i_layer][i_chip][i_channel];
        // _temp_graph = g_temperature[i_layer][i_chip][i_channel];


        // fr[0]=0.3*h_MIP_tempCorrected[i_layer][i_chip][i_channel]->GetMean();
        // fr[1]=3.0*h_MIP_tempCorrected[i_layer][i_chip][i_channel]->GetMean();
        // if(fitting){
        // TF1 *fitsnr_temp = langaufit(h_MIP_tempCorrected[i_layer][i_chip][i_channel],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
        // langaupro(fp,SNRPeak,SNRFWHM);
        // // Global style settings
        // gStyle->SetOptStat(1111);
        // gStyle->SetOptFit(111);
        // gStyle->SetLabelSize(0.03,"x");
        // gStyle->SetLabelSize(0.03,"y");
        // }
        // if(i_layer<4 || i_layer>28) langaus(h_MIP_tempCorrected[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+200, HG_ped[i_layer][i_chip][i_channel]+1000);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 50, 400);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 400, 800);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP_tempCorrected[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+50, HG_ped[i_layer][i_chip][i_channel]+400);
        // _landauMPV_tempCorrected = langaus_peakP;
        // _landauWidth_tempCorrected = langaus_peakPError;
        // _ChiSqr_tempCorrected = langaus_func->GetChisquare();
        // _Ndf_tempCorrected = langaus_func->GetNDF();
        // _ChiSqr_tempCorrected = chisqr;
        // _Ndf_tempCorrected = ndf;
        //_gausSigma = fun->GetParameter(3);
        // _gausSigma_tempCorrected = 0;
        // _landauWidth_tempCorrected = fp[0];
        // _landauMPV_tempCorrected = fp[1];
        // _gausSigma_tempCorrected = fp[3];
        // _ChiSqr_tempCorrected = chisqr;
        // _Ndf_tempCorrected = ndf;
        // _MIP_histo_tempCorrected = h_MIP_tempCorrected[i_layer][i_chip][i_channel];
        
        mip_tree->Fill();
      }
      C_MIP[i_layer][i_chip]->Update();
      C_MIP[i_layer][i_chip]->Write(Form("Layer%d_Chip%d", i_layer, i_chip));
      C_MIP[i_layer][i_chip]->SaveAs(Form("../result/mip/figures/chips/layer%dchip%d.png",i_layer,i_chip));
    }
  }
  mip_tree->Write();

  TCanvas* c_MIPs = new TCanvas(Form("c_MIPs"), Form("c_MIPs"), 2560, 2560);
  TCanvas* c_MIPs_width = new TCanvas(Form("c_MIPs_width"), Form("c_MIPs_width"), 2560, 2560);
  TCanvas* c_MIPs_GauSigma = new TCanvas(Form("c_MIPs_GauSigma"), Form("c_MIPs_GauSigma"), 2560, 2560);
  TCanvas* c_chi2ndf = new TCanvas(Form("c_chi2ndf"), Form("c_chi2ndf"), 2560, 2560/2);
  c_MIPs->Divide(2,2);
  c_MIPs_width->Divide(2,2);
  c_MIPs_GauSigma->Divide(2,2);
  c_chi2ndf->Divide(2,1);
  for(int pix=0;pix<2;pix++){
    MIPs[pix]->Write();
    MIPs_width[pix]->Write();
    MIPs_GauSigma[pix]->Write();
    MIPs_error[pix]->Write();
    MIPs_width_error[pix]->Write();
    MIPs_GauSigma_error[pix]->Write();

    c_MIPs->cd();
    c_MIPs->cd(pix+1);
    MIPs[pix]->Draw();
    c_MIPs->cd(pix+3);
    MIPs_error[pix]->Draw();

    c_MIPs_width->cd();
    c_MIPs_width->cd(pix+1);
    MIPs_width[pix]->Draw();
    c_MIPs_width->cd(pix+3);
    MIPs_width_error[pix]->Draw();

    c_MIPs_GauSigma->cd();
    c_MIPs_GauSigma->cd(pix+1);
    MIPs_GauSigma[pix]->Draw();
    c_MIPs_GauSigma->cd(pix+3);
    MIPs_GauSigma_error[pix]->Draw();

    c_chi2ndf->cd();
    c_chi2ndf->cd(pix+1);
    // gPad->SetLogy(1);
    // chi2ndf[pix]->SetLogy();
    chi2ndf[pix]->Draw();
  }
  c_MIPs->SaveAs(Form("../result/mip/figures/mip.png"));
  c_MIPs_width->SaveAs(Form("../result/mip/figures/mip_width.png"));
  c_MIPs_GauSigma->SaveAs(Form("../result/mip/figures/mip_gausigma.png"));
  c_chi2ndf->SaveAs(Form("../result/mip/figures/chi2ndf.png"));


  // runごとのnumber of hitの確認用  20011, 12,13,14,18
  // TGraph *graph_nhit_run[5];
  // for(int iii=0; iii<5; iii++) {
  //   graph_nhit_run[iii] = new TGraph();
  //   for(int irawfile=0; irawfile<rawfilenum/2; irawfile++) {
  //     graph_nhit_run[iii]->SetPoint(irawfile, irawfile, Nhit_run[iii][irawfile]);
  //     // if(Nhit_run[iii][irawfile]>1000) cout << argv[irawfile*2+2] << endl;
  //   }
  //   graph_nhit_run[iii]->SetMarkerStyle(4);
  //   graph_nhit_run[iii]->Write(Form("graph_nhit_run_%d",iii));
  // }

}