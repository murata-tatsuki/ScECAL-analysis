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
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/src/EBUdecode.cxx"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;


// ssa後にcutをかけてMPVを求めるマクロ
// pedestal を引いている

const bool cut = true;	// cutをするかどうかのbool
const bool chip_cut_nai = false;	// 同一chipに他にhitが"ない"channelの集計	1chn
const bool chip_cut_aru = false;	// 同一chipに他にhitが"ある"channelの集計	2chn
const bool auto_gain = false;	// auto gain mode の解析用
const bool combined_gain = true;	// 両方のgain modeの解析用
const bool temperature_correction = true;	// 温度補正「も」するかどうか


int main(int argc, char* argv[])
{ 
    if(argc < 3 && argc%2==1){                                                     //エラー処理
        cout << "usage: ./between_files output.root input1.root input2.root  ..." << endl;
        return 1;
    }

		if(chip_cut_nai && chip_cut_aru){
			cout << "chip_cut_nai && chip_cut_aru are both true is forbidden" << endl;
			return 1;
		}

		if(auto_gain && combined_gain){
			cout << "auto_gain && combined_gain are both true is forbidden" << endl;
			return 1;
		}

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;
  
	int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	// vector<pair<string,int>> filename;
	int total_entry_max=0;

  // double rawFileDate[rawfilenum];

	// for(int i=0; i<rawfilenum; i++){
	// 	string txtname(argv[i+2]);
	// 	txtname.erase(0,txtname.find_last_of("/")+1);
	// 	// txtname.erase(txtname.find_last_of("."),txtname.end()-1);
	// 	filename.push_back({txtname, i+2});
	// }
  // sort(filename.begin(),filename.end());
  for(int i=0; i<rawfilenum; i++){
		// filein[i] = new TFile(argv[filename[i].second]);
		filein[i] = new TFile(argv[i+2]);
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
	vector<vector<double>> *temp = nullptr;
	// vector<double> *posX = nullptr;
	// vector<double> *posY = nullptr;
	// vector<double> *posZ = nullptr;


	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;




	// pedestal
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





	TFile fileout(argv[1],"RECREATE");
	// fileout.mkdir("chip_canvas");
	// fileout.mkdir("chip_canvas_tempCorrected");
	// fileout.mkdir("chip_canvas_run");
	// fileout.mkdir("chip_canvas_tempCorrected_run");
	// gain_histo->Write();



	// dataを入れるもの
  
	

	// runごとにつくるもの


	// 保存形式
  double _slope, _intercept, _slopeError, _interceptError, _slope_tempCorrected, _intercept_tempCorrected, _slopeError_tempCorrected, _interceptError_tempCorrected, _ChiSqr;
	int _chn_Entries;
	int _cellID_fit;
	int _Ndf;
	TGraph *_hl_graph = new TGraph();
	TGraph *_hl_graph_tempCorrected = new TGraph();

	TTree* hl_tree = new TTree("T_Event", "after calibration");
  hl_tree->Branch("cellIDs",&_cellID_fit);
	hl_tree->Branch("chnEntries", &_chn_Entries);
  hl_tree->Branch("slope",&_slope);
  hl_tree->Branch("intercept",&_intercept);
  hl_tree->Branch("slopeError",&_slopeError);
  hl_tree->Branch("interceptError",&_interceptError);
  hl_tree->Branch("slope",&_slope_tempCorrected);
  hl_tree->Branch("intercept",&_intercept_tempCorrected);
  hl_tree->Branch("slopeError",&_slopeError_tempCorrected);
  hl_tree->Branch("interceptError",&_interceptError_tempCorrected);
  hl_tree->Branch("chiSquare",&_ChiSqr);
  hl_tree->Branch("NDF",&_Ndf);
  hl_tree->Branch("hlGraph", &_hl_graph);
  hl_tree->Branch("hlGraph_tempCorrected", &_hl_graph_tempCorrected);


  Double_t langaus_peakP, langaus_peakPError;
	TF1 *langaus_func;



	cout << "made save TTree" << endl;


  char char_tmp[100];
  // TGraph* g_hl[layerNu][chipNu][channelNu];
  // TGraph* g_hl_tempCorrected[layerNu][chipNu][channelNu];
  // int g_hl_count[layerNu][chipNu][channelNu];
  // // TH1F* h_MIP_run[layerNu][chipNu][channelNu];
  // // TH1F* h_MIP_tempCorrected_run[layerNu][chipNu][channelNu];
  // for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
	// 	for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
	// 		for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
	// 			sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);

	// 			g_hl[i_layer][i_chip][i_channel] = new TGraph();
	// 			g_hl_tempCorrected[i_layer][i_chip][i_channel] = new TGraph();

  //       g_hl_count[i_layer][i_chip][i_channel] = 0;
	// 		}
	// 	}
	// }





	// data をとってきてる

  cout << "saving histograms ... " << endl;


	// data の保存をしている
  TCanvas* C_hl[layerNu][chipNu];
  TCanvas* C_hl_tempCorrected[layerNu][chipNu];
  TF1 *fit_pol1 = new TF1("fit_pol1","[0]*x+[1]",1000,2500);
  fit_pol1->SetParameter(0,0.03);
  fit_pol1->SetParameter(1,300);
  fit_pol1->SetParLimits(0,0,1);
  fit_pol1->SetParLimits(1,100,500);
  fit_pol1->SetLineColor(2);
  fit_pol1->SetLineWidth(2);
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    cout << i_layer << "/" << layerNu << endl;
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      cout << i_chip << "/" << chipNu << endl;
      sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
      C_hl[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
      C_hl[i_layer][i_chip]->Divide(6, 6);
			sprintf(char_tmp, "layer%d_chip%d_tempCorrected", i_layer, i_chip);
      C_hl_tempCorrected[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
      C_hl_tempCorrected[i_layer][i_chip]->Divide(6, 6);

      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        if (i_chip == 5 && i_channel > 29) continue;




        TGraph* g_hl = new TGraph();
        TGraph* g_hl_tempCorrected = new TGraph();
        int g_hl_count = 0;
        
      	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
      		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

      		tree[irawfile]->SetBranchAddress("triggerID", &triggerID);
      		tree[irawfile]->SetBranchAddress("BCIDs", &BCIDs);
      		tree[irawfile]->SetBranchAddress("cellIDs", &cellIDs);
      		tree[irawfile]->SetBranchAddress("charges", &charges);
      		tree[irawfile]->SetBranchAddress("times", &times);
      		tree[irawfile]->SetBranchAddress("hitTags", &hitTags);
      		tree[irawfile]->SetBranchAddress("gainTags", &gainTags);
      		tree[irawfile]->SetBranchAddress("temp", &temp);
      		// t[irawfile]->SetBranchAddress("posX", &posX);
      		// t[irawfile]->SetBranchAddress("posY", &posY);
      		// t[irawfile]->SetBranchAddress("posZ", &posZ);

      		for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
            // if(irawfile==17) cout << ientry << "/" << entry_max[irawfile] << endl;
      			tree[irawfile]->GetEntry(ientry);
            if(cellIDs->size()==0) continue;

      			for(int i=0; i<(int)cellIDs->size(); i++){
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
      				if(!(_layer==i_layer && _chip==i_chip && _channel==i_channel)) continue;

      				double *_position = EBUdecode(_layer,_chip,_channel);
              double SiPMtemp;
              if(temp->at(_layer).size()!=0) SiPMtemp = tempReconstruction(_layer, _position, temp->at(_layer));
              else SiPMtemp=20;
      				double slope_now = gain_slope[_layer][_chip][_channel];
      				double intercept_now = gain_intercept[_layer][_chip][_channel];
      				double _charge_corrected = (_charges - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];
      				double _time_corrected = (_times - HG_ped[_layer][_chip][_channel]) * (slope_now*20.0+intercept_now)/(slope_now*SiPMtemp+intercept_now)  + HG_ped[_layer][_chip][_channel];

      				// h_temperature[_layer][_chip][_channel]->Fill(SiPMtemp);

      				// if(_layer==0 && _chip==2 && _channel==11) cout << slope_now << ", " << SiPMtemp << ", " << intercept_now << ", " << _times << ", " << HG_ped[_layer][_chip][_channel] << ", " << _time_corrected << ", " << slope_now*SiPMtemp+intercept_now << "," << slope_now*20.0+intercept_now << endl;

      				if(_charges<=600 && _times>HG_ped[_layer][_chip][_channel]+50) g_hl->SetPoint(g_hl_count,_times,_charges);
      				if(_charge_corrected<=600 && _time_corrected>HG_ped[_layer][_chip][_channel]+50) g_hl_tempCorrected->SetPoint(g_hl_count,_time_corrected,_charge_corrected);
              g_hl_count++;
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
      	}






























        _chn_Entries = g_hl->GetN();
        _cellID_fit = i_layer * 1e6 + i_chip * 1e4 + i_channel;

        if(i_layer>3&&i_layer<28){
          g_hl->GetXaxis()->SetRangeUser(350,1500);
          g_hl->GetYaxis()->SetRangeUser(300,600);
          g_hl->Fit(fit_pol1,"qnROB=0.75","",1000,2500);
          // fit_pol1->SetRange(1000,2500);
        }
        else{
          g_hl->GetXaxis()->SetRangeUser(100,4000);
          g_hl->GetYaxis()->SetRangeUser(350,600);
          g_hl->Fit(fit_pol1,"qnROB=0.75","",1000,2500);
          // fit_pol1->SetRange(1000,2500);
        }
        C_hl[i_layer][i_chip]->cd(i_channel + 1);
        g_hl->Draw("AP");
        fit_pol1->Draw("same");
        _slope = fit_pol1->GetParameter(0);
        _intercept = fit_pol1->GetParameter(1);
        _slopeError = fit_pol1->GetParError(0);
        _interceptError = fit_pol1->GetParError(1);
        _hl_graph = g_hl;
        _ChiSqr = fit_pol1->GetChisquare();
        _Ndf = fit_pol1->GetNDF();
        

        if(i_layer>3&&i_layer<28){
          g_hl_tempCorrected->GetXaxis()->SetRangeUser(350,1500);
          g_hl_tempCorrected->GetYaxis()->SetRangeUser(300,600);
          g_hl_tempCorrected->Fit(fit_pol1,"qnROB=0.75","",500,1400);
          fit_pol1->SetRange(100,1000);
        }
        else{
          g_hl_tempCorrected->GetXaxis()->SetRangeUser(100,4000);
          g_hl_tempCorrected->GetYaxis()->SetRangeUser(350,600);
          g_hl_tempCorrected->Fit(fit_pol1,"qnROB=0.75","",600,2600);
          fit_pol1->SetRange(200,2200);
        }
        C_hl_tempCorrected[i_layer][i_chip]->cd(i_channel + 1);
        g_hl_tempCorrected->Draw("AP");
        fit_pol1->Draw("same");
        _slope_tempCorrected = fit_pol1->GetParameter(0);
        _intercept_tempCorrected = fit_pol1->GetParameter(1);
        _slopeError_tempCorrected = fit_pol1->GetParError(0);
        _interceptError_tempCorrected = fit_pol1->GetParError(1);
				_hl_graph_tempCorrected = g_hl_tempCorrected;

        hl_tree->Fill();

        g_hl->Delete();
        g_hl_tempCorrected->Delete();
      }
      C_hl[i_layer][i_chip]->Write();
      C_hl_tempCorrected[i_layer][i_chip]->Write();
    }
  }
  hl_tree->Write();



}
