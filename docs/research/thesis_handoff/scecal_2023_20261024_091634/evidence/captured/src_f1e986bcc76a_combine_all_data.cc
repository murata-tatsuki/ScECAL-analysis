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




int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./combine_all_data output.root input1.root input1.root  ..." << endl;
        return 1;
		}

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;
  
	int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	int total_entry_max=0;
	
  for(int i=0; i<rawfilenum; i++){
		// filein[i] = new TFile(argv[filename[i].second]);
		filein[i] = new TFile(argv[i+2]);
		tree[i] = (TTree*) filein[i]->Get("Calib_Hit");
		entry_max[i] = tree[i]->GetEntries();
	}

	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;


	TFile fileout(argv[1],"RECREATE");



	double TotalEnergyDep;
	vector<int> *CellID=nullptr;


  // TH1F * energyHisto = new TH1F("Ehisto","Ehisto",1000,0,1000);
  // TH2F *E_nhit = new TH2F("E_nhit",";MeV;numer of hit",200,0,200, 200,0,200);
	TH1F * energyHisto = new TH1F("Ehisto","Ehisto",2000,0,2000);
  TH2F *E_nhit = new TH2F("E_nhit",";MeV;numer of hit",2000,0,2000, 3000,0,3000);




	// data をとってきてる
	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;
    vector<vector<int>> hitStrip;
	
		// ssa
		tree[irawfile]->SetBranchAddress("TotalEnergyDep", &TotalEnergyDep);
		tree[irawfile]->SetBranchAddress("CellID", &CellID);

		for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);

      energyHisto->Fill(TotalEnergyDep);
			E_nhit->Fill(TotalEnergyDep,CellID->size());
		}
		filein[irawfile]->Close();
  }

  energyHisto->Write();
	E_nhit->Write();


}