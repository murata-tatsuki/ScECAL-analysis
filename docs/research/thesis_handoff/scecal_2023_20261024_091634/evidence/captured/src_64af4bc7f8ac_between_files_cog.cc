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
#include <unordered_set>
#include <cmath>
#include "TLegend.h"
#include "TCanvas.h"
#include "langaus.C"				// langaus(chargeH, &fitFunc, &peakP, &peakPError);
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// エネルギー分布を求めるマクロ
// ssaを使う


const unordered_set<int> bad_channels = {
  // bad channel の CellID をここに列挙
};


int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files output.root Nenegry E1...En E1Nfiles...EnNfiles input1.root input2.root  ..." << endl;
        return 1;
    }

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;

  int Nenegry = atoi(argv[2]);
  double energy[Nenegry];
  int energy_files[Nenegry];
  for(int i=0;i<Nenegry;i++){
    energy[i] = atoi(argv[3+i]);
    if(energy[i]<1) energy[i] = 0.5;
    energy_files[i] = atoi(argv[3+Nenegry+i]);
  }
  
  
	int rawfilenum = argc - 2 - 2*Nenegry;
  int irawfilenum = 3 + 2*Nenegry;


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

	// double cycleID_ssa, triggerID_ssa;
	// vector<int> *cellIDs_ssa = nullptr;
	// vector<int> *BCIDs_ssa = nullptr;
	// vector<int> *hitTags_ssa = nullptr;
	// vector<int> *gainTags_ssa = nullptr;
	// vector<double> *charges_ssa = nullptr;
	// vector<double> *times_ssa = nullptr;
	// vector<double> *temp_ssa = nullptr;
	// vector<double> *posX_ssa = nullptr;
	// vector<double> *posY_ssa = nullptr;
	// vector<double> *posZ_ssa = nullptr;
	// vector<int> *ssaTag_ssa = nullptr;

  vector<int> *cellID = nullptr;
  vector<int> *ssaTag = nullptr;
  vector<double> *Hit_Energy = nullptr;
  vector<double> *Hit_X = nullptr;
	vector<double> *Hit_Y = nullptr;
	vector<double> *Hit_Z = nullptr;
	


	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;

  const int range_maximum = 12000;
  const int range_maximum_PS = 500;
  const double binWidth_PS = 0.1;



	TFile fileout(argv[1],"RECREATE");
  fileout.mkdir("all__channel");
  fileout.mkdir("best_channel");
	// gain_histo->Write();



	// dataを入れるもの
  TH1F* energy_deposition_bestChannels[Nenegry];
  TH1F* energy_deposition[Nenegry];
  TH2F* hitmap[Nenegry][layerNu];
  TH2F* cog[Nenegry][layerNu];
  for (int i = 0; i < Nenegry; ++i) {
    if(energy[i]>5) energy_deposition[i] = new TH1F(Form("energy_deposition_%d",i),Form("%d GeV energy deposition",(int)energy[i]),range_maximum,0,range_maximum);
    else if(energy[i]<1) energy_deposition[i] = new TH1F(Form("energy_deposition_%d",i),"0.5 GeV energy deposition",range_maximum_PS/binWidth_PS,0,range_maximum_PS);
    else energy_deposition[i] = new TH1F(Form("energy_deposition_%d",i),Form("%d GeV energy deposition",(int)energy[i]),range_maximum_PS/binWidth_PS,0,range_maximum_PS);
    energy_deposition[i]->SetXTitle("Reconstructed Energy [MeV]");

    if(energy[i]>5) energy_deposition_bestChannels[i] = new TH1F(Form("energy_deposition_bestChannels_%d",i),Form("%d GeV energy deposition",(int)energy[i]),range_maximum,0,range_maximum);
    else if(energy[i]<1) energy_deposition_bestChannels[i] = new TH1F(Form("energy_deposition_bestChannels_%d",i),"0.5 GeV energy deposition",range_maximum_PS/binWidth_PS,0,range_maximum_PS);
    else energy_deposition_bestChannels[i] = new TH1F(Form("energy_deposition_bestChannels_%d",i),Form("%d GeV energy deposition",(int)energy[i]),range_maximum_PS/binWidth_PS,0,range_maximum_PS);
    energy_deposition_bestChannels[i]->SetXTitle("Detected Energy [MeV]");

    for (int ilayer = 0; ilayer < layerNu; ++ilayer) {
      // hitmap[i] = new TH2F(Form("hitmap_%dGeV_layer",i,ilayer),Form("hitmap %d GeV layer %d",(int)energy[i],ilayer),range_maximum,0,range_maximum);
      cog[i][ilayer] = new TH2F(Form("cog_%dGeV_layer_%d",i,ilayer),Form("center of gravity %d GeV layer %d",(int)energy[i],ilayer),226,-113,113,226,-113,113);
    }
	}
  
	

	// runごとにつくるもの
	


  



	// 保存形式
  



	cout << "made save TTree" << endl;



  // vector<string> skipFiles{"ECAL_Run196_20230504_123821.root", "ECAL_Run197_20230504_135137.root", "ECAL_Run194_20230504_082851.root", "ECAL_Run195_20230504_104136.root", "AHCAL_Run133_20230501_081117.root"};


	// data をとってきてる
	for(int ienegry=0; ienegry<Nenegry; ienegry++){
		if(Nenegry>1) cout << energy[ienegry] << " GeV" << endl;

    TFile *filein[energy_files[ienegry]];
    TTree *tree[energy_files[ienegry]];
    int entry_max[energy_files[ienegry]];
    for(int i=0; i<energy_files[ienegry]; i++){
	  	filein[i] = new TFile(argv[irawfilenum+i]);
	  	tree[i] = (TTree*) filein[i]->Get("SSA_Hit");
	  	entry_max[i] = tree[i]->GetEntries();
    }

    for(int irawfile=0; irawfile<energy_files[ienegry]; irawfile++){
      cout << irawfile << "/" << energy_files[ienegry] << endl;

      // string fileName = argv[irawfilenum+irawfile];
      // auto result = find(skipFiles.begin(), skipFiles.end(), fileName);
      // if (result != skipFiles.end()) continue;
      
		  tree[irawfile]->SetBranchAddress("CellID", &cellID);
		  tree[irawfile]->SetBranchAddress("Hit_Energy", &Hit_Energy);
		  tree[irawfile]->SetBranchAddress("ssaTag", &ssaTag);
		  tree[irawfile]->SetBranchAddress("Hit_X", &Hit_X);
		  tree[irawfile]->SetBranchAddress("Hit_Y", &Hit_Y);
		  tree[irawfile]->SetBranchAddress("Hit_Z", &Hit_Z);

	  	for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
        // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
	  		tree[irawfile]->GetEntry(ientry);
        
        double cog_x[32] = {0};
        double cog_y[32] = {0};
        double nhits[32] = {0};
        for(int ihit=0;ihit<Hit_Energy->size();ihit++){
          int _layerID = cellID->at(ihit)/1e5;
          if(_layerID>=30 || ssaTag->at(ihit)==0) continue;
          const int currentCellId = cellID->at(ihit);
          if(bad_channels.find(currentCellId) != bad_channels.end()) continue;
          
          cog_x[_layerID] += Hit_X->at(ihit) * Hit_Energy->at(ihit);
          cog_y[_layerID] += Hit_Y->at(ihit) * Hit_Energy->at(ihit);
          nhits[_layerID]++;
        }

        for (int ilayer = 0; ilayer < layerNu; ++ilayer) {
          cog[ienegry][ilayer]->Fill(cog_x[ilayer]/nhits[ilayer], cog_y[ilayer]/nhits[ilayer]);
        }
	  	}

      cellID = nullptr;
      Hit_Energy = nullptr;
      ssaTag = nullptr;
      Hit_X = nullptr;
      Hit_Y = nullptr;
      Hit_Z = nullptr;
      delete filein[irawfile];
	  }

    irawfilenum += energy_files[ienegry];
  }

  cout << "fitting histograms ... " << endl;

  fileout.cd();
  gStyle->SetOptFit(1111111);













  fileout.cd("all__channel");
  TCanvas *C_cog[Nenegry];
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    C_cog[ienegry] = new TCanvas(Form("C_cog_%d",ienegry), Form("center of gravity %d GeV",(int)energy[ienegry]), 1) ;
    C_cog[ienegry]->Divide(6,5);
    C_cog[ienegry]->cd();
    for (int ilayer = 0; ilayer < 30; ++ilayer) {
      C_cog[ienegry]->cd(ilayer+1);
      cog[ienegry][ilayer]->Draw("colz");
    }
    C_cog[ienegry]->Write(Form("C_cog_%dGeV",(int)energy[ienegry]));
  }


}