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


void MIP_comparison()
{ 
    // if(argc != 4){                                                     //エラー処理
    //     cout << "usage: ./MIP_comparison output.root input_data.root input_simulation.root" << endl;
    //     return 1;
    // }

  // cout << "=====>  " << argv[1] << endl;
  
  int rawfilenum = 2;
  // int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
  TTree *tree[rawfilenum];
  int entry_max[rawfilenum];
  int total_entry_max=0;
	
  // for(int i=0; i<rawfilenum; i++){
    // filein[1] = new TFile(argv[3]);
    // tree[1] = (TTree*) filein[1]->Get("Digi_Raw_Hit");
    // entry_max[1] = tree[1]->GetEntries();
  // }


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

  


  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;







  TFile fileout("../results/mip/simulation/cosmicThreScan/100GeV/comparison.root","RECREATE");
  fileout.mkdir("chip_canvas");
  // fileout.mkdir("chip_canvas_tempCorrected");
  // fileout.mkdir("chip_canvas_run");
  // fileout.mkdir("chip_canvas_tempCorrected_run");
  // gain_histo->Write();
  

    TFile *jiaxuanMIPHistoFile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip_old/jiaxuan/all_auto_muon_v2_hist.root");
    // TFile *jiaxuanMIPFitFile = new TFile("/home/murata_t/ScECAL_BeamTest/ECAL/results/mip_old/jiaxuan/all_auto_muon_v2_fit.root");
    // TDirectoryFile *tdict_layer[layerNu];
    // TDirectoryFile *tdict_layer_chip[layerNu][chipNu];
    // TH2D *hMIPChn = (TH2D*)jiaxuanMIPHistoFile->Get("hMIPChn");
    TH1D *MIP_histo_jiaxuan[layerNu][chipNu][channelNu];
    jiaxuanMIPHistoFile->cd();
    for(int i=0; i<layerNu;i++){
        for(int j=0; j<chipNu;j++){
            for(int k=0; k<channelNu;k++){
                if(j==5&&k>29) continue;
                MIP_histo_jiaxuan[i][j][k] = (TH1D*)jiaxuanMIPHistoFile->Get(Form("layer%d/chip%d/Layer%d_Chip%d_ Chn%d",i,j,i,j,k));
                if(MIP_histo_jiaxuan[i][j][k] == nullptr) continue;

                int nbin = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetNbins();
                double histo_min = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetXmin();
                double histo_max = MIP_histo_jiaxuan[i][j][k]->GetXaxis()->GetXmax();
                int binWidth = (histo_max - histo_min) / nbin;

                if(nbin!=400 || histo_min!=0 || histo_max!=4000) cout << "  histo bin error " << i << ", " << j << ", " << k << endl;

                MIP_histo_jiaxuan[i][j][k]->Scale(1.0/MIP_histo_jiaxuan[i][j][k]->GetEntries());
            }
        }
    }
    cout << "sps mip (jiaxuan) reading finished" << endl;










  // simulation
  int irawfile = 1;
  char char_tmp[100];
  filein[1] = new TFile("/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_100000_thre_0.root");
    tree[1] = (TTree*) filein[1]->Get("Digi_Raw_Hit");
    entry_max[1] = tree[1]->GetEntries();
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


    tree[irawfile]->SetBranchAddress("TriggerID", &triggerID);
    tree[irawfile]->SetBranchAddress("BCID", &BCIDs);
    tree[irawfile]->SetBranchAddress("CellID", &cellIDs);
    tree[irawfile]->SetBranchAddress("LG_Charge", &charges);
    tree[irawfile]->SetBranchAddress("HG_Charge", &times);
    tree[irawfile]->SetBranchAddress("HitTag", &hitTags);
    tree[irawfile]->SetBranchAddress("GainTag", &gainTags);

    for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      tree[irawfile]->GetEntry(ientry);
      if(cellIDs->size()==0) continue;

      for(int i=0; i<(int)cellIDs->size(); i++){
        int _layer = cellIDs->at(i)/1e5;
        int _chip = (cellIDs->at(i)%100000)/1e4;
        int _channel = cellIDs->at(i)%100;
        double _charges = charges->at(i);
        double _times = times->at(i);
        int _hitTags = hitTags->at(i);
        int _gainTags = gainTags->at(i);

        if(_hitTags==0) continue;

        h_MIP[_layer][_chip][_channel]->Fill(_times);
      }
    }
    triggerID = 0;
    BCIDs = nullptr;
    cellIDs = nullptr;
    charges = nullptr;
    times = nullptr;
    hitTags = nullptr;
    gainTags = nullptr;
    // filein[irawfile]->Close();

    cout << "simulation mip reading finished" << endl;







  for(int i=0; i<layerNu;i++){
    for(int j=0; j<chipNu;j++){
      for(int k=0; k<channelNu;k++){
        if(j==5&&k>29) continue;
        if(h_MIP[i][j][k]->GetEntries()<500) continue;
        // cout << i << ", " << j << ", " << k << "   " << h_MIP[i][j][k]->GetEntries() << endl;

        h_MIP[i][j][k]->Scale(1.0/h_MIP[i][j][k]->GetEntries());
        if(MIP_histo_jiaxuan[i][j][k] == nullptr) continue;

        double difference = comparison(MIP_histo_jiaxuan[i][j][k], h_MIP[i][j][k]);

        cout << "   " << i << ", " << j << ", " << k << "  difference :" << difference << endl;

      }
    }
  }

  

}