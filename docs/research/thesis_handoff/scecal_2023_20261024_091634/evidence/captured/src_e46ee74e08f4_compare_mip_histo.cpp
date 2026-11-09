// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/src/EBUdecode.cpp"
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
#include "langaus.C"				// langaus(chargeH, &fitFunc, &peakP, &peakPError);

// 適当なchannel(今は20011,12,13,14,18)のMIPのhistoを温度補正有り無しで比べるマクロ

void compare_mip_histo(){
  const bool cut = true;	// cutをするかどうかのbool
  const bool chip_cut_nai = false;	// 同一chipに他にhitが"ない"channelの集計	1chn
  const bool chip_cut_aru = false;	// 同一chipに他にhitが"ある"channelの集計	2chn
  const bool auto_gain = false;	// auto gain mode の解析用
  const bool combined_gain = true;	// 両方のgain modeの解析用

  

  ifstream ifs("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/ECAL/results/MIP_gain_allfilename.txt");
	string line, stream, field;
	vector<string> fileNames;
  int rawfilenum = 0;

	while(getline(ifs, line)){
    istringstream stream(line);
    while(getline(stream, field, ' ')){
      fileNames.push_back(field);
      rawfilenum++;
    }
	}
	ifs.close();
  
  
  TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	// vector<pair<string,int>> filename;
	int total_entry_max=0;


  for(int i=0; i<rawfilenum; i++){
		// filein[i] = new TFile(argv[filename[i].second]);
		filein[i] = new TFile(fileNames[i].c_str());
		tree[i] = (TTree*) filein[i]->Get("T_Event");
		entry_max[i] = tree[i]->GetEntries();

    // string txtname(argv[filename[i].second]);
    // // txtname.erase(0,txtname.find_last_of("/")+1);
    // // txtname.erase(0,txtname.find_first_of("_")+1);
    // // txtname.erase(0,txtname.find_first_of("_")+1);
    // // txtname.erase(txtname.find_last_of("."),5);
    // // txtname.erase(8,1);
    // // txtname.erase(0,4);
    // // cout << stod(txtname.c_str()) / 1000000.0 << endl;
    // // rawFileDate[i] = stod(txtname.c_str()) / 1000000.0;
    // txtname.erase(0,txtname.find_last_of("/")+1);
    // txtname.erase(0,8);
    // txtname.erase(txtname.find_first_of("_"),txtname.size());
    // rawFileDate[i] = atoi(txtname.c_str());

	}


	double cycleID, triggerID;
	vector<int> *cellIDs = nullptr;
	vector<int> *BCIDs = nullptr;
	vector<int> *hitTags = nullptr;
	vector<int> *gainTags = nullptr;
	vector<double> *charges = nullptr;
	vector<double> *times = nullptr;
	// vector<double> *temp = nullptr;
	// vector<double> *posX = nullptr;
	// vector<double> *posY = nullptr;
	// vector<double> *posZ = nullptr;

	double cycleID_ssa, triggerID_ssa;
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





	TFile *fileinPed;
  fileinPed = new TFile("/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/ped/fit/pedestal/ECAL_Run23_20230518_153022.root");
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




	TFile fileout("mip_compare.root","RECREATE");
	// gain_histo->Write();



	// dataを入れるもの
  
	

	// runごとにつくるもの


	



	// // MIPがきちんと評価できてないチャンネル
	// // 統計不足
	// 0,0,0
	// 0,3,odd
	// 0,3,32
	// 0,3,34
	// 0,4,0
	// 0,4,1
	// 0,4,2
	// 0,4,3
	// 0,4,5

	// // noiseと分離ができてない
	// 0,0,3

	// // fitがうまくできてない
	// 0,1,35


  



	// 保存形式
  double _landauMPV, _landauWidth, _gausSigma, _ChiSqr;
	int _chn_Entries;
	int _cellID_fit;
	int _Ndf;  
	TH1F *_MIP_histo = new TH1F();

	TTree* mip_tree = new TTree("T_Event", "after calibration");
	mip_tree->Branch("chnEntries", &_chn_Entries);
	mip_tree->Branch("cellIDs", &_cellID_fit);
	mip_tree->Branch("landauMPV", &_landauMPV);
	mip_tree->Branch("landauWidth", &_landauWidth);
	mip_tree->Branch("gausSigma", &_gausSigma);
	mip_tree->Branch("ChiSqr", &_ChiSqr);
	mip_tree->Branch("ndf", &_Ndf);
	mip_tree->Branch("histo", &_MIP_histo);


  Double_t langaus_peakP, langaus_peakPError;
	TF1 *langaus_func;



	cout << "made save TTree" << endl;



  char char_tmp[100];
  TH1F* h_MIP[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
				if (3 < i_layer && i_layer < 28)
					h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -50, 1000);
				else
					h_MIP[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 150, -100, 3000);
				h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
				h_MIP[i_layer][i_chip][i_channel]->SetLineWidth(2);
				h_MIP[i_layer][i_chip][i_channel]->SetLineColor(4);
			}
		}
	}

	// runごとのnumber of hit の確認用
	int Nhit_run[5][rawfilenum/2];		// 20011, 12,13,14,18 の5つのchannelをまずは確認する(中心のchannel)
	for(int irawfile=0; irawfile<rawfilenum/2; irawfile++) {
		for(int iii=0; iii<5; iii++) {
			Nhit_run[iii][irawfile] = 0;
		}
	}






	// data をとってきてる
	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;
    vector<vector<int>> hitStrip;
	
		// ssa
		tree[irawfile]->SetBranchAddress("triggerID", &triggerID_ssa);
		tree[irawfile]->SetBranchAddress("BCIDs", &BCIDs_ssa);
		tree[irawfile]->SetBranchAddress("cellIDs", &cellIDs_ssa);
		tree[irawfile]->SetBranchAddress("charges", &charges_ssa);
		tree[irawfile]->SetBranchAddress("times", &times_ssa);
		tree[irawfile]->SetBranchAddress("hitTags", &hitTags_ssa);
		tree[irawfile]->SetBranchAddress("gainTags", &gainTags_ssa);
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
      v.push_back((int)BCIDs_ssa->at(0));
      v.push_back((int)triggerID_ssa);
			for(int i=0; i<(int)cellIDs_ssa->size(); i++){
				// cout << "   " << i << "/" << cellIDs_ssa->size() << endl;
				int _layer = cellIDs_ssa->at(i)/1000000;
				int _chip = (cellIDs_ssa->at(i)/10000)%100;
				int _channel = cellIDs_ssa->at(i)%100;
				double _charges = charges_ssa->at(i);
				double _times = times_ssa->at(i);
				int _hitTags = hitTags_ssa->at(i);
				// double _x = posX->at(i);
				// double _y = posY->at(i);
				// double _z = posZ->at(i);
				int _ssaTag = ssaTag_ssa->at(i);
				if(_hitTags==0 || _ssaTag==0) continue;
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


  // decode
	// for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		irawfile++;
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

		tree[irawfile]->SetBranchAddress("triggerID", &triggerID);
		tree[irawfile]->SetBranchAddress("BCIDs", &BCIDs);
		tree[irawfile]->SetBranchAddress("cellIDs", &cellIDs);
		tree[irawfile]->SetBranchAddress("charges", &charges);
		tree[irawfile]->SetBranchAddress("times", &times);
		tree[irawfile]->SetBranchAddress("hitTags", &hitTags);
		tree[irawfile]->SetBranchAddress("gainTags", &gainTags);
		// t[irawfile]->SetBranchAddress("temp", &temp);
		// t[irawfile]->SetBranchAddress("posX", &posX);
		// t[irawfile]->SetBranchAddress("posY", &posY);
		// t[irawfile]->SetBranchAddress("posZ", &posZ);

		for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      // if(irawfile==17) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);
      if(cellIDs->size()==0) continue;

      int itr = 0;
      for(itr=0;itr<(int)hitStrip.size();itr++){
        if(hitStrip[itr][0]==BCIDs->at(0) && hitStrip[itr][1]==triggerID) break;
      }
      if(cut) if(itr==hitStrip.size()) continue;

			// int nchip = cellIDs->size()/channelNu;
			int nhit_chip[32][6] = {0};
			if(chip_cut_nai || chip_cut_aru){
				// for(int i=0;i<nchip;i++) nhit_chip[i] = 0;
				for(int i=0; i<(int)cellIDs->size(); i++){
					int _layer = cellIDs->at(i)/1000000;
					int _chip = (cellIDs->at(i)/10000)%100;
					int _hitTags = hitTags->at(i);
					if(_hitTags!=0) nhit_chip[_layer][_chip]++;
				} 
			}

			for(int i=0; i<(int)cellIDs->size(); i++){
        auto result = find(hitStrip[itr].begin()+2, hitStrip[itr].end(), (int)cellIDs->at(i));
        if(cut) if(result == hitStrip[itr].end()) continue;
				// if(chip_cut && nhit_chip[i/36]>1) continue; 
				// cout << "   " << i << "/" << cellIDs->size() << endl;
				int _layer = cellIDs->at(i)/1000000;
				int _chip = (cellIDs->at(i)/10000)%100;
				int _channel = cellIDs->at(i)%100;
				double _charges = charges->at(i);
				double _times = times->at(i);
				int _hitTags = hitTags->at(i);
				int _gainTags = gainTags->at(i);
				// double _x = posX->at(i);
				// double _y = posY->at(i);
				// double _z = posZ->at(i);

				if(_hitTags==0) continue;
				if(chip_cut_nai && nhit_chip[_layer][_chip]>1) continue; 
				if(chip_cut_aru && nhit_chip[_layer][_chip]<=1) continue; 

				// runごとのnumber of hitの確認用
				if(cellIDs->at(i)==20011) Nhit_run[0][(int)irawfile/2]++;
				if(cellIDs->at(i)==20012) Nhit_run[1][(int)irawfile/2]++;
				if(cellIDs->at(i)==20013) Nhit_run[2][(int)irawfile/2]++;
				if(cellIDs->at(i)==20014) Nhit_run[3][(int)irawfile/2]++;
				if(cellIDs->at(i)==20018) Nhit_run[4][(int)irawfile/2]++;

        // if(_gainTags==1) h_MIP[_layer][_chip][_channel]->Fill(_charges);
				if(auto_gain && _gainTags == 1) h_MIP[_layer][_chip][_channel]->Fill(_charges);
				if(combined_gain) {
					if(irawfile>=156 && _gainTags == 1) h_MIP[_layer][_chip][_channel]->Fill(_charges);
					if(irawfile<156) h_MIP[_layer][_chip][_channel]->Fill(_times);
				}
				if(!auto_gain && !combined_gain)  h_MIP[_layer][_chip][_channel]->Fill(_times);
			}
		}
    triggerID = 0;
		BCIDs = nullptr;
		cellIDs = nullptr;
		charges = nullptr;
		times = nullptr;
		hitTags = nullptr;
		gainTags = nullptr;
	}

  cout << "saving histograms ... " << endl;


	// data の保存をしている
  TCanvas* C_MIP[layerNu][chipNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    cout << i_layer << "/" << layerNu << endl;
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
      C_MIP[i_layer][i_chip]->Divide(6, 6);

      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if (i_chip == 5 && i_channel > 29) continue;
        _chn_Entries = h_MIP[i_layer][i_chip][i_channel]->GetEntries();
        _cellID_fit = i_layer * 1e6 + i_chip * 1e4 + i_channel;

        if(i_layer<4 || i_layer>28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+200, HG_ped[i_layer][i_chip][i_channel]+1000);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 50, 400);
        // if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, 400, 800);
        if(i_layer>=4 && i_layer<=28) langaus(h_MIP[i_layer][i_chip][i_channel],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+50, HG_ped[i_layer][i_chip][i_channel]+400);

				_landauWidth = langaus_peakPError;
				_landauMPV = langaus_peakP;
				//_gausSigma = fun->GetParameter(3);
				_gausSigma = 0;
				_ChiSqr = langaus_func->GetChisquare();
				_Ndf = langaus_func->GetNDF();
				// cout<<"--------  "<<_landauMPV<<" , "<<_landauWidth<<" , "<<_gausSigma<<endl;
				_MIP_histo = h_MIP[i_layer][i_chip][i_channel];
        mip_tree->Fill();

        C_MIP[i_layer][i_chip]->cd(i_channel + 1);
      	h_MIP[i_layer][i_chip][i_channel]->Draw();
        langaus_func->Draw("same");
      }
      C_MIP[i_layer][i_chip]->Write();
    }
  }
  mip_tree->Write();




	// runごとのnumber of hitの確認用  20011, 12,13,14,18
	TGraph *graph_nhit_run[5];
	for(int iii=0; iii<5; iii++) {
		graph_nhit_run[iii] = new TGraph();
		for(int irawfile=0; irawfile<rawfilenum/2; irawfile++) {
			graph_nhit_run[iii]->SetPoint(irawfile, irawfile, Nhit_run[iii][irawfile]);
			// if(Nhit_run[iii][irawfile]>1000) cout << argv[irawfile*2+2] << endl;
		}
		graph_nhit_run[iii]->SetMarkerStyle(4);
		graph_nhit_run[iii]->Write(Form("graph_nhit_run_%d",iii));
	}




  
}