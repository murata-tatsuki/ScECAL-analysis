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
        cout << "usage: ./combine_all_data output.root input_ssa.root input_calib.root  ..." << endl;
        return 1;
		}

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;
  
	int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	int total_entry_max=0;
	
  // for(int i=0; i<rawfilenum; i++){
	// 	// filein[i] = new TFile(argv[filename[i].second]);
	// 	filein[i] = new TFile(argv[i+2]);
	// 	tree[i] = (TTree*) filein[i]->Get("Calib_Hit");
	// 	entry_max[i] = tree[i]->GetEntries();
	// }

	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;


	TFile fileout(argv[1],"RECREATE");



	double TotalEnergyDep;
	
	int triggerID;
	vector<int> *cellIDs = nullptr;
	vector<double> *charges = nullptr;

	int triggerID_ssa;
	vector<int> *cellIDs_ssa = nullptr;
	vector<double> *charges_ssa = nullptr;
	vector<int> *ssaTag_ssa = nullptr;


  TH1F * energyHisto = new TH1F("Ehisto","Ehisto;MeV;",1000,0,1000);
  TH2F *E_nhit = new TH2F("E_nhit",";MeV;numer of hit",200,0,200, 200,0,200);




	vector<vector<int>> hitStrip;

	// data をとってきてる
	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;
    // vector<vector<int>> hitStrip;
		

		// ssa
		filein[irawfile] = new TFile(argv[irawfile+2]);
		tree[irawfile] = (TTree*) filein[irawfile]->Get("SSA_Hit");

		tree[irawfile]->SetBranchAddress("Event_Num", &triggerID_ssa);
		tree[irawfile]->SetBranchAddress("CellID", &cellIDs_ssa);
		tree[irawfile]->SetBranchAddress("ssaTag", &ssaTag_ssa);

		for(int ientry=0; ientry<tree[irawfile]->GetEntries(); ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);
      vector<int> v;
      int layerHit[32] = {0};
      int nlayerHit = 0;

      if(cellIDs_ssa->size()==0) continue;
      v.push_back((int)triggerID_ssa);
			for(int i=0; i<(int)cellIDs_ssa->size(); i++){
				// cout << "   " << i << "/" << cellIDs_ssa->size() << endl;
				int _layer = cellIDs_ssa->at(i)/1e5;
				int _chip = (cellIDs_ssa->at(i)/10000)%100;
				int _channel = cellIDs_ssa->at(i)%100;
				int _ssaTag = ssaTag_ssa->at(i);
				if(_ssaTag==0) continue;
				auto result = find(v.begin()+1, v.end(), (int)cellIDs_ssa->at(i));
				if (result == v.end()) {
          layerHit[_layer]++;
					v.push_back((int)cellIDs_ssa->at(i));
				}
			}
      for(int i=0;i<32;i++) if(layerHit[i]>0) nlayerHit++;
      if(v.size()>=22 && v.size()<64 && nlayerHit>=22) hitStrip.push_back(v);
		}
		triggerID_ssa = 0;
		cellIDs_ssa = nullptr;
		ssaTag_ssa = nullptr;
		filein[irawfile]->Close();



		// calib
		irawfile++;
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << "   " << argv[irawfile+2] << endl;

		filein[irawfile] = new TFile(argv[irawfile+2]);
		tree[irawfile] = (TTree*) filein[irawfile]->Get("Calib_Hit");

		tree[irawfile]->SetBranchAddress("Event_Num", &triggerID);
		tree[irawfile]->SetBranchAddress("CellID", &cellIDs);
		tree[irawfile]->SetBranchAddress("Hit_Energy", &charges);
		tree[irawfile]->SetBranchAddress("TotalEnergyDep", &TotalEnergyDep);
int coutcheck=0;
int runNo=-1;
if(irawfile==runNo) cout << coutcheck++ << endl;
		for(int ientry=0; ientry<tree[irawfile]->GetEntries(); ientry++){
      // if(irawfile==17) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);
if(irawfile==runNo) cout << coutcheck++ << endl;
      if(cellIDs->size()==0) continue;
if(irawfile==runNo) cout << coutcheck++ << endl;


      int itr = 0;
      for(itr=0;itr<(int)hitStrip.size();itr++){
        if(hitStrip[itr][0]==triggerID) break;
      }
if(irawfile==runNo) cout << coutcheck++ << endl;
if(irawfile==runNo) cout << itr << endl;

			double _totalEnergy = 0;
			int _totalHit = 0;
			for(int i=0; i<(int)cellIDs->size(); i++){
				if(irawfile==runNo){
					cout << cellIDs->at(i) << endl;
					cout << hitStrip[itr].size() << endl;
					for(int iiitr=0;iiitr<(int)hitStrip[itr].size();iiitr++){
						cout << hitStrip[itr][iiitr] << endl;
					}
				}
				auto result = find(hitStrip[itr].begin()+1, hitStrip[itr].end(), (int)cellIDs->at(i));
				if(result == hitStrip[itr].end()) {
					if(irawfile==runNo) cout << "skipped";
					continue;
				}
				
				// int _layer = cellIDs->at(i)/1e5;
				// int _chip = (cellIDs->at(i)/10000)%100;
				// int _channel = cellIDs->at(i)%100;
				_totalEnergy += charges->at(i);
				_totalHit++;
			}
if(irawfile==runNo) cout << coutcheck++ << endl;
			// energyHisto->Fill(TotalEnergyDep);
			// E_nhit->Fill(TotalEnergyDep,cellIDs->size());
			if(_totalHit>22){
				energyHisto->Fill(_totalEnergy);
				E_nhit->Fill(_totalEnergy,_totalHit);
			}
		}
    triggerID = 0;
		TotalEnergyDep = 0;
		cellIDs = nullptr;
		charges = nullptr;
		filein[irawfile]->Close();

		hitStrip.clear();

  }


	fileout.cd();
  energyHisto->Write();
	E_nhit->Write();


}