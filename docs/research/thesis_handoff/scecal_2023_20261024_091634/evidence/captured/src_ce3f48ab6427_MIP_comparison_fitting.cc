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
#include "langaus_newConvoluted.C"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;

// ssa後にcutをかけてMPVを求めるマクロ
// pedestal を引いている

const bool cut = true;  // cutをするかどうかのbool

const bool fitting = false;


double comparison(TH1D *h1, TH1D *h2){
  int nbin = h1->GetXaxis()->GetNbins();
  double diff = 0;
  for(int ibin=0;ibin<nbin;ibin++){
    double dif = pow(h1->GetBinContent(ibin) - h2->GetBinContent(ibin), 2);
    diff += dif;
  }
  return diff;
}


int main(int argc, char* argv[])
{ 
    if(argc < 4){                                                     //エラー処理
        cout << "usage: ./MIP_comparison output.root input_data.root input_simulation1.root, input_simulation2.root ... " << endl;
        return 1;
    }

  // cout << "=====>  " << argv[1] << endl;
  
  // int rawfilenum = 2;
  int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  int total_entry_max=0;
  int threshold[rawfilenum];
	threshold[0] = -1;
  for(int i=1; i<rawfilenum; i++){
    filein[i] = new TFile(argv[i+2]);
    tree[i] = (TTree*) filein[i]->Get("T_Event");
    entry_max[i] = tree[i]->GetEntries();

    string txtname(argv[i+2]);
    txtname.erase(0,txtname.find_last_of("_")+1);
    cout << argv[i+2] << ",  " << stoi(txtname) << endl;
    threshold[i] = stoi(txtname);
  }

  // int cycleID, triggerID;
  // vector<int> *cellIDs = nullptr;
  // vector<int> *BCIDs = nullptr;
  // vector<int> *hitTags = nullptr;
  // vector<int> *gainTags = nullptr;
  // vector<double> *charges = nullptr;
  // vector<double> *times = nullptr;
  // vector<vector<double>> *temp = nullptr;
  // vector<double> *posX = nullptr;
  // vector<double> *posY = nullptr;
  // vector<double> *posZ = nullptr;

  int cellIDs, chnEntries, ndf;
  double landauMPV=0, landauWidth, gausSigma, ChiSqr;
  TH1F *histo = new TH1F();

  


  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;

  double mip_0[layerNu][chipNu][channelNu];
  double mip_100[layerNu][chipNu][channelNu];
  double entry_100[layerNu][chipNu][channelNu];






  int _threshold, _cellIDs, _chnEntries, _ndf;
  double _landauMPV=0, _landauWidth, _gausSigma, _ChiSqr;
  TH1D *h_data;
  TH1F *h_sim;
  TFile fileout(argv[1],"RECREATE");
  TTree* treeout = new TTree("T_Event", "after comparison");
  treeout->Branch("threshold", &_threshold);
  treeout->Branch("chnEntries", &_chnEntries);
  treeout->Branch("CellID", &_cellIDs);
  treeout->Branch("landauMPV", &_landauMPV);
  treeout->Branch("landauWidth", &_landauWidth);
  treeout->Branch("gausSigma", &_gausSigma);
  treeout->Branch("ChiSqr", &_ChiSqr);
  treeout->Branch("ndf", &_ndf);
  treeout->Branch("histo_simulation", &h_sim);
  

    TFile *jiaxuanMIPHistoFile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip_old/jiaxuan/all_auto_muon_v2_hist.root");
    // TFile *jiaxuanMIPFitFile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip_old/jiaxuan/all_auto_muon_v2_fit.root");
    // TDirectoryFile *tdict_layer[layerNu];
    // TDirectoryFile *tdict_layer_chip[layerNu][chipNu];
    TH2D *hMIPChn = (TH2D*)jiaxuanMIPHistoFile->Get("hMIPChn");
    double MIP_jiaxuan[layerNu][chipNu][channelNu];
    TH1D *MIP_histo_jiaxuan[layerNu][chipNu][channelNu];
    jiaxuanMIPHistoFile->cd();
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
              MIP_jiaxuan[i][j][k] = 0;
              mip_0[i][j][k] = 0;
              mip_100[i][j][k] = 0;
              entry_100[i][j][k] = 0;
                if(j==5&&k>29) continue;
                MIP_histo_jiaxuan[i][j][k] = (TH1D*)jiaxuanMIPHistoFile->Get(Form("layer%d/chip%d/Layer%d_Chip%d_ Chn%d",i,j,i,j,k));
                if(MIP_histo_jiaxuan[i][j][k] == nullptr) continue;

                int nbin = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetNbins();
                double histo_min = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetXmin();
                double histo_max = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetXmax();
                int binWidth = (histo_max - histo_min) / nbin;

                // if(nbin!=400 || histo_min!=0 || histo_max!=4000) cout << "  histo bin error " << i << ", " << j << ", " << k << endl;

                // MIP_histo_jiaxuan[i][j][k]->Scale(1.0/MIP_histo_jiaxuan[i][j][k]->GetEntries());
                // MIP_histo_jiaxuan[i][j][k] = MIP_histo_jiaxuan[i][j][k]->Scale(1.0/MIP_histo_jiaxuan[i][j][k]->GetEntries());
                MIP_jiaxuan[i][j][k] = hMIPChn->GetBinContent(i*chipNu + j + 1, k);
            }
        }
    }
    cout << "sps mip (jiaxuan) reading finished" << endl;




  vector<pair<double,int>> dif_channel[layerNu][chipNu][channelNu];




    char char_tmp[100];
  TH1D* h_MIP[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
        // if (3 < i_layer && i_layer < 28)
        //   // h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -50, 1000);
        //   h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 210, -50, 1000);
        // else
          // h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -100, 3000);
          h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 400, 0, 4000);
        h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
        h_MIP[i_layer][i_chip][i_channel]->SetLineWidth(2);
        h_MIP[i_layer][i_chip][i_channel]->SetLineColor(4);
      }
    }
  }
  // simulation
  for(int irawfile=1; irawfile<rawfilenum; irawfile++){
  // int irawfile = 1;
    // filein[irawfile] = new TFile("/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_100000_thre_0.root");
    // filein[irawfile] = new TFile(argv[i+2]);
    // tree[irawfile] = (TTree*) filein[irawfile]->Get("Digi_Raw_Hit");
    // entry_max[irawfile] = tree[irawfile]->GetEntries();
    


    tree[irawfile]->SetBranchAddress("chnEntries", &chnEntries);
    tree[irawfile]->SetBranchAddress("cellIDs", &cellIDs);
    tree[irawfile]->SetBranchAddress("landauMPV", &landauMPV);
    tree[irawfile]->SetBranchAddress("landauWidth", &landauWidth);
    tree[irawfile]->SetBranchAddress("gausSigma", &gausSigma);
    tree[irawfile]->SetBranchAddress("ChiSqr", &ChiSqr);
    tree[irawfile]->SetBranchAddress("ndf", &ndf);
    tree[irawfile]->SetBranchAddress("histo", &histo);

    _threshold = threshold[irawfile];

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);
      // if(cellIDs->size()==0) continue;

      // for(int i=0; i<(int)cellIDs->size(); i++){
        int _layer = cellIDs/1e6;
        int _chip = (cellIDs%1000000)/1e4;
        int _channel = cellIDs%100;
        // cout << cellIDs << "  ";
        // double _charges = charges->at(i);
        // double _times = times->at(i);
        // int _hitTags = hitTags->at(i);
        // int _gainTags = gainTags->at(i);
        // if(_hitTags==0) continue;

        // h_MIP[_layer][_chip][_channel]->Fill(_times);
        _chnEntries = chnEntries;
        _cellIDs = cellIDs;
        _landauMPV = landauMPV;
        _landauWidth = landauWidth;
        _gausSigma = gausSigma;
        _ChiSqr = ChiSqr;
        _ndf = ndf;
        h_sim = histo;
        treeout->Fill();

        // cout << _layer << " " << _chip << " " << _channel << endl;
        if(_threshold==0) mip_0[_layer][_chip][_channel] = landauMPV;
        if(_threshold==100){
          mip_100[_layer][_chip][_channel] = landauMPV;
          entry_100[_layer][_chip][_channel] = chnEntries;
        }
      // }
    }
    // triggerID = 0;
    // BCIDs = nullptr;
    // cellIDs = nullptr;
    // charges = nullptr;
    // times = nullptr;
    // hitTags = nullptr;
    // gainTags = nullptr;
    histo = nullptr;

    cout << "simulation mip reading finished  threshold " << threshold[irawfile] << endl;






    // _threshold = threshold[irawfile];
    // for(int i=0; i<layerNu;i++){
    //   for(int j=0; j<chipNu;j++){
    //     for(int k=0; k<channelNu;k++){
    //       if(j==5&&k>29) continue;
    //       if(h_MIP[i][j][k]->GetEntries()<500) continue;
    //       // cout << i << ", " << j << ", " << k << "   " << h_MIP[i][j][k]->GetEntries() << endl;
    //       h_MIP[i][j][k]->Scale(1.0/h_MIP[i][j][k]->GetEntries());
    //       if(MIP_histo_jiaxuan[i][j][k] == nullptr) continue;

    //       int range_maximum = (i>=4 && i<28) ? 2500 : 1200;
    //       h_MIP[i][j][k]->GetXaxis()->SetRangeUser(300,range_maximum);

    //       double difference = comparison(MIP_histo_jiaxuan[i][j][k], h_MIP[i][j][k]);
    //       // cout << "   " << i << ", " << j << ", " << k << "  difference :" << difference << endl;

          
    //       CellID_out = i*1e5 + j*1e4 + k;
    //       difference_out = difference;
    //       h_data = MIP_histo_jiaxuan[i][j][k];
    //       h_sim = h_MIP[i][j][k];

    //       treeout->Fill();
    //       h_MIP[i][j][k]->Reset();
    //       h_MIP[i][j][k]->ResetStats();
    //       dif_channel[i][j][k].push_back({difference,threshold[irawfile]});
    //     }
    //   }
    // }



    filein[irawfile]->Close();
  }
  fileout.cd();
  treeout->Write();


  TH1F *histo_ratio_10 = new TH1F("histo_ratio_10","histo_ratio_10",200,0,2);
  TH1F *histo_ratio_15 = new TH1F("histo_ratio_15","histo_ratio_15",200,0,2);
  for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
              // cout << mip_100[i][j][k] << " " << mip_0[i][j][k] << endl;
              if(entry_100[i][j][k]>1000 && mip_0[i][j][k]>0){
                if(i<4 || i>=28) histo_ratio_15->Fill(mip_100[i][j][k]/mip_0[i][j][k]);
                else histo_ratio_10->Fill(mip_100[i][j][k]/mip_0[i][j][k]);
              }
            }
        }
  }
  histo_ratio_10->Write();
  histo_ratio_15->Write();

}