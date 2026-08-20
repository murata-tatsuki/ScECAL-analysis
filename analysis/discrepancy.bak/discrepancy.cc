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
#include <set>
#include "TCanvas.h"
#include "langaus.C"

using namespace std;


// エネルギー分布を求めるマクロ
// calibを使う

int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files output.root N_data N_sim datafile_1.root...datafile_N_data.root datafile_1.root...datafile_N_sim.root" << endl;
        // cout << "usage: ./between_files output.root Nenegry E1...En E1Nfiles...EnNfiles input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;

  int N_data = atoi(argv[2]);
  int N_sim = atoi(argv[3]);
  
	// double cycleID, triggerID;
	// vector<int> *cellIDs = nullptr;
	// vector<int> *BCIDs = nullptr;
	// vector<int> *hitTags = nullptr;
	// vector<int> *gainTags = nullptr;
	// vector<double> *charges = nullptr;
	// vector<double> *times = nullptr;
	// vector<double> *temp = nullptr;
	// vector<double> *posX = nullptr;
	// vector<double> *posY = nullptr;
	// vector<double> *posZ = nullptr;

	// double cycleID_sim, triggerID_sim;
	// vector<int> *cellIDs_sim = nullptr;
	// vector<int> *BCIDs_sim = nullptr;
	// vector<int> *hitTags_sim = nullptr;
	// vector<int> *gainTags_sim = nullptr;
	// vector<double> *charges_sim = nullptr;
	// vector<double> *times_sim = nullptr;
	// vector<double> *temp_sim = nullptr;
	// vector<double> *posX_sim = nullptr;
	// vector<double> *posY_sim = nullptr;
	// vector<double> *posZ_sim = nullptr;
	// vector<int> *ssaTag_sim = nullptr;

  double EnergyDep;

  vector<int> *CellID_data = nullptr;
  vector<double> *Hit_Energy_data = nullptr;

  vector<int> *CellID_sim = nullptr;
  vector<double> *Hit_Energy_sim = nullptr;
	


	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;

  const int range_maximum = 10;
  // const int range_maximum_PS = 500;
  // const double binWidth_PS = 0.1;



	TFile fileout(argv[1],"RECREATE");
	// gain_histo->Write();


	// dataを入れるもの
  char char_tmp[100];
  TH1F* edep_data[layerNu][chipNu][channelNu];
  TH1F* edep_sim[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        sprintf(char_tmp, "data_Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
        edep_data[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, range_maximum*100,0,range_maximum);
        sprintf(char_tmp, "sim_Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
        edep_sim[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp,range_maximum*100,0,range_maximum);
      }
    }
  }

  // 除外するcellID（必要に応じて追加）
  set<int> excludedCellIDs;
  // excludedCellIDs.insert(10000);  // 例: Layer0 Chip1 Chn0 を除外

  TH1F* h_total_edep_data = new TH1F("total_edep_data", "Total Hit Energy per Event (data)", 12000, 0, 12000);
  TH1F* h_total_edep_sim = new TH1F("total_edep_sim", "Total Hit Energy per Event (sim)", 12000, 0, 12000);

	cout << "made save TTree" << endl;


  // vector<string> skipFiles{"ECAL_Run196_20230504_123821.root", "ECAL_Run197_20230504_135137.root", "ECAL_Run194_20230504_082851.root", "ECAL_Run195_20230504_104136.root", "AHCAL_Run133_20230501_081117.root"};


	// data をとってきてる
  TFile *filein_data[N_data];
  TTree *tree_data[N_data];
  for(int irawfile=0; irawfile<N_data; irawfile++){
    cout << "data : " << irawfile << "/" << N_data << endl;

    filein_data[irawfile] = new TFile(argv[4+irawfile]);
		tree_data[irawfile] = (TTree*) filein_data[irawfile]->Get("Calib_Hit");
	  // tree_data[irawfile]->SetBranchAddress("TotalEnergyDep", &EnergyDep);
	  tree_data[irawfile]->SetBranchAddress("CellID", &CellID_data);
	  tree_data[irawfile]->SetBranchAddress("Hit_Energy", &Hit_Energy_data);
	  for(int ientry=0; ientry<tree_data[irawfile]->GetEntries(); ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
	  	tree_data[irawfile]->GetEntry(ientry);

      double sum_edep = 0;
      for(int i=0; i<(int)CellID_data->size(); i++){
        int _cellID = CellID_data->at(i);
        int _layer = _cellID/1e5;
        int _chip = (_cellID%100000)/1e4;
        int _channel = _cellID%100;
        double _charges = Hit_Energy_data->at(i);

        edep_data[_layer][_chip][_channel]->Fill(_charges);
        if (!excludedCellIDs.count(_cellID)) sum_edep += _charges;
      }
      h_total_edep_data->Fill(sum_edep);
      
	  }
    delete filein_data[irawfile];
  }

  TFile *filein_sim[N_sim];
  TTree *tree_sim[N_sim];
  for(int irawfile=0; irawfile<N_sim; irawfile++){
    cout << "simulation : " << irawfile << "/" << N_sim << endl;

    filein_sim[irawfile] = new TFile(argv[4+irawfile]);
		tree_sim[irawfile] = (TTree*) filein_sim[irawfile]->Get("Calib_Hit");
	  // tree_sim[irawfile]->SetBranchAddress("TotalEnergyDep", &EnergyDep);
	  tree_sim[irawfile]->SetBranchAddress("CellID", &CellID_sim);
	  tree_sim[irawfile]->SetBranchAddress("Hit_Energy", &Hit_Energy_sim);
	  for(int ientry=0; ientry<tree_sim[irawfile]->GetEntries(); ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
	  	tree_sim[irawfile]->GetEntry(ientry);

      double sum_edep = 0;
      for(int i=0; i<(int)CellID_sim->size(); i++){
        int _cellID = CellID_sim->at(i);
        int _layer = _cellID/1e5;
        int _chip = (_cellID%100000)/1e4;
        int _channel = _cellID%100;
        double _charges = Hit_Energy_sim->at(i);

        edep_sim[_layer][_chip][_channel]->Fill(_charges);
        if (!excludedCellIDs.count(_cellID)) sum_edep += _charges;
      }
      h_total_edep_sim->Fill(sum_edep);
      
	  }
    delete filein_sim[irawfile];
  }

  cout << "fitting histograms ... " << endl;

  fileout.cd();
  // gStyle->SetOptFit();

  TCanvas* c_edep[layerNu][chipNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      c_edep[i_layer][i_chip] = new TCanvas( Form("edep_Layer%d_Chip%d", i_layer, i_chip), Form("Layer%d Chip%d", i_layer, i_chip), 1800, 1200);
      c_edep[i_layer][i_chip]->Divide(6, 6);

      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        c_edep[i_layer][i_chip]->cd(i_channel + 1);

        // TH1F* h_data = edep_data[i_layer][i_chip][i_channel];
        // TH1F* h_sim = edep_sim[i_layer][i_chip][i_channel];

        edep_data[i_layer][i_chip][i_channel]->SetLineColor(kBlue + 1);
        edep_data[i_layer][i_chip][i_channel]->SetLineWidth(2);
        sprintf(char_tmp, "Layer%d Chip%d Chn%d", i_layer, i_chip, i_channel);
        edep_data[i_layer][i_chip][i_channel]->SetTitleSize(0.1);
        edep_data[i_layer][i_chip][i_channel]->SetTitle(char_tmp);
        edep_data[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("[MeV]]");
        edep_data[i_layer][i_chip][i_channel]->GetYaxis()->SetTitle("[a.u.]");
        edep_data[i_layer][i_chip][i_channel]->GetXaxis()->SetTitleSize(0.08);
        edep_data[i_layer][i_chip][i_channel]->GetXaxis()->SetTitleOffset(0.6);
        edep_data[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.07);
        edep_data[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelOffset(0.005);
        edep_data[i_layer][i_chip][i_channel]->GetYaxis()->SetTitleSize(0.08);
        edep_data[i_layer][i_chip][i_channel]->GetYaxis()->SetTitleOffset(0.6);
        edep_data[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.07);
        edep_data[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelOffset(0.005);
        edep_data[i_layer][i_chip][i_channel]->Scale(1./edep_data[i_layer][i_chip][i_channel]->GetEntries());
        edep_sim[i_layer][i_chip][i_channel]->SetLineColor(kRed + 1);
        edep_sim[i_layer][i_chip][i_channel]->SetLineWidth(2);
        edep_sim[i_layer][i_chip][i_channel]->Scale(1./edep_data[i_layer][i_chip][i_channel]->GetEntries());

        const double max_y = max(edep_data[i_layer][i_chip][i_channel]->GetMaximum(), edep_sim[i_layer][i_chip][i_channel]->GetMaximum());
        if (max_y > 0) {
          edep_data[i_layer][i_chip][i_channel]->Scale(max_y*1.15/edep_data[i_layer][i_chip][i_channel]->GetEntries());
          edep_sim[i_layer][i_chip][i_channel]->Scale(max_y*1.15/edep_sim[i_layer][i_chip][i_channel]->GetEntries());
        }

        edep_data[i_layer][i_chip][i_channel]->Draw("hist");
        edep_sim[i_layer][i_chip][i_channel]->Draw("hist same");
      }
      c_edep[i_layer][i_chip]->Write();
    }
  }

  h_total_edep_data->Write();
  h_total_edep_sim->Write();



}