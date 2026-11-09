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

double* EBUdecode_itr(int LayerIDs, int ChipIDs, int ChannelIDs){
	// static const int layerNu = 32;
	static const int chipNu = 6;
	static const int chnNu = 36;
	int decodeID[chipNu][chnNu] = { 0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13, 55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,
																	24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,
																	149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
																	150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
																	191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181, 139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
																	137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, 210, 210, 210, 210, 210, 210 };

	int ScintillatorIDs = decodeID[ChipIDs][ChannelIDs];
	double layerZ;

	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	int _yID = ScintillatorIDs / rowNu;
	int _xID = ScintillatorIDs % rowNu;
	static double _position[3] = { 0 };
	double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
	double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;

	// for prototype test
	if (LayerIDs % 2 == 0) {
		_position[0] = -y0;
		_position[1] = -x0;
	}
	if (LayerIDs % 2 == 1) {
		_position[0] = -x0;
		_position[1] = -y0;
	}
	if (LayerIDs % 2 == 0)
		layerZ = 1 + LayerIDs / 2 * 19.9;
	//else layerZ = 12.95+(LayerIDs-1)/2*19.9;
	else
		layerZ = 12.2 + (LayerIDs - 1) / 2 * 19.9;
	_position[2] = layerZ;
	// for grouped CR test
	// if(LayerIDs%4==0) {
	//     _position[0] = -y0;
	//     _position[1] = -x0;
	// }
	// if(LayerIDs%4==1) {
	//     _position[0] = x0;
	//     _position[1] = y0;
	// }
	// if(LayerIDs%4==2) {
	//     _position[0] = y0;
	//     _position[1] = x0;
	// }
	// if(LayerIDs%4==3) {
	//     _position[0] = -x0;
	//     _position[1] = -y0;
	// }
	// _position[2] = layerZ[LayerIDs];

	return _position;
}

void check_threshold(){

  const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int doublelayerNu = 2;
	const int rowNu = 42;
	const int columnNu = 5;
  const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction

  double constantFraction = 0.15;
  double thresholds[layerNu][chipNu][channelNu];

  // pedestal 
  TFile *fileinPed;
  fileinPed = new TFile("/megraid01/users/data_beamtest/ECAL_data/analysed_oldFormat/2023/ps/ped/fit/pedestal/ECAL_Run23_20230518_153022.root");
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
        thresholds[i_layer][i_chip][i_channel] = 0;
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


  double landauMPV, landauWidth, gausSigma, ChiSqr, landauWidth_tempCorrected, gausSigma_tempCorrected, landauMPV_tempCorrected, ChiSqr_tempCorrected;
	int chnEntries, cellIDs, ndf, ndf_tempCorrected;
	TH1F *histo = new TH1F();
	TH1F *histo_tempCorrected = new TH1F();
	TGraphErrors *tempGraph = new TGraphErrors();

	T_Event->SetBranchAddress("chnEntries", &chnEntries);
	T_Event->SetBranchAddress("cellIDs", &cellIDs);
	T_Event->SetBranchAddress("landauMPV", &landauMPV);
	T_Event->SetBranchAddress("landauWidth", &landauWidth);
	T_Event->SetBranchAddress("gausSigma", &gausSigma);
	T_Event->SetBranchAddress("ChiSqr", &ChiSqr);
	T_Event->SetBranchAddress("ndf", &ndf);
	T_Event->SetBranchAddress("histo", &histo);
	T_Event->SetBranchAddress("landauWidth_tempCorrected", &landauWidth_tempCorrected);
	T_Event->SetBranchAddress("gausSigma_tempCorrected", &gausSigma_tempCorrected);
	T_Event->SetBranchAddress("landauMPV_tempCorrected", &landauMPV_tempCorrected);
	T_Event->SetBranchAddress("ChiSqr_tempCorrected", &ChiSqr_tempCorrected);
	T_Event->SetBranchAddress("ndf_tempCorrected", &ndf_tempCorrected);
	T_Event->SetBranchAddress("histo_tempCorrected", &histo_tempCorrected);
	T_Event->SetBranchAddress("tempGraph", &tempGraph);



  TH1F *thre_hist = new TH1F("thre","thre",200,0,1000);

  for(int ientry=0;ientry<T_Event->GetEntries();ientry++){
		T_Event->GetEntry(ientry);

			int _layer = cellIDs/1000000;
			int _chip = (cellIDs/10000)%100;
			int _channel = cellIDs%100;

      double *_position = EBUdecode_itr(_layer,_chip,_channel);
      double _x = *_position;
      double _y = *(_position+1);
      double _z = *(_position+2);
      int _xID;
	    int _yID;
      if(_layer%2==0) {
	      _xID = (_x + _yInterval * (columnNu - 1) / 2.) / _yInterval + 0.5;
        _yID = (_y + _xInterval * (rowNu - 1) / 2.) / _xInterval + 0.5;
      }
      if(_layer%2==1) {
        _xID = (_x + _xInterval * (rowNu - 1) / 2.) / _xInterval + 0.5;
	      _yID = (_y + _yInterval * (columnNu - 1) / 2.) / _yInterval + 0.5;
      }


      histo_tempCorrected->Rebin(5);
      // int binWidth = histo_tempCorrected->GetXaxis()->GetBinWidth();
      int nbin = histo_tempCorrected->GetXaxis()->GetNbins();
      double histo_min = histo_tempCorrected->GetXaxis()->GetXmin();
      double histo_max = histo_tempCorrected->GetXaxis()->GetXmax();
      int binWidth = (histo_max - histo_min) / nbin;
      // if( binWidth != (histo_max-histo_min)/nbin ) cout << "something is wrong at entry:" << ientry << " cellID:" << cellIDs << "  " << nbin << ", " << histo_min << ", " << histo_max << endl;

      double maximum = histo_tempCorrected->GetBinContent(histo_tempCorrected->GetMaximumBin());
      int threBin = histo_tempCorrected->FindFirstBinAbove(maximum * constantFraction);
      double thre = binWidth * threBin + histo_min;
      thresholds[_layer][_chip][_channel] = thre;

      thre_hist->Fill(thre);


  }

  thre_hist->Draw();




  // TH2F *MIP = new TH2F("MIP","MIP", 192,0,192, 36,0,36);
  // TH2F *MIP_woPed = new TH2F("MIP_woPed","MIP_woPed", 192,0,192, 36,0,36);
  // TH2F *MIP_diff = new TH2F("MIP_diff","MIP_diff", 192,0,192, 36,0,36);
  // TH1F *MIP_diff_1d = new TH1F("MIP_diff_1d","MIP_diff_1d", 300,-1500,1500);
  // TH1F *MIP_15 = new TH1F("MIP_15","MIP_15", 90,0,900);
  // TH1F *MIP_10 = new TH1F("MIP_10","MIP_10", 90,0,900);
  // TGraph *MIP_entry_15 = new TGraph();
  // TH2F *NuEntry[32];
  // for(int i=0;i<32;i++){
  //   if(i%2==0) NuEntry[i] = new TH2F(Form("NuEntry_%d",i),Form("NuEntry_%d",i), 5,-113.5,113.5, 42,-113.5,113.5);
  //   if(i%2==1) NuEntry[i] = new TH2F(Form("NuEntry_%d",i),Form("NuEntry_%d",i), 42,-113.5,113.5, 5,-113.5,113.5);
  // }
  // int count_out_of_range = 0;
  // int count_in_range = 0;
  // int count_dead = 0;
  // int count_MIP_entry_15 = 0;
	// for(int ientry=0;ientry<T_Event->GetEntries();ientry++){
	// 	T_Event->GetEntry(ientry);

	// 		int _layer = cellIDs/1000000;
	// 		int _chip = (cellIDs/10000)%100;
	// 		int _channel = cellIDs%100;

  //     double *_position = EBUdecode_itr(_layer,_chip,_channel);
  //     double _x = *_position;
  //     double _y = *(_position+1);
  //     double _z = *(_position+2);
  //     int _xID;
	//     int _yID;
  //     if(_layer%2==0) {
	//       _xID = (_x + _yInterval * (columnNu - 1) / 2.) / _yInterval + 0.5;
  //       _yID = (_y + _xInterval * (rowNu - 1) / 2.) / _xInterval + 0.5;
  //     }
  //     if(_layer%2==1) {
  //       _xID = (_x + _xInterval * (rowNu - 1) / 2.) / _xInterval + 0.5;
	//       _yID = (_y + _yInterval * (columnNu - 1) / 2.) / _yInterval + 0.5;
  //     }

  //     // NuEntry[_layer]->SetBinContent(_xID+1,_yID+1, chnEntries);

  //     // if(_layer==0 && _chip==3) MIP_entry_15->SetPoint(count_MIP_entry_15++, chnEntries, landauMPV - HG_pedMean);

  //     // MIP->SetBinContent(_layer*6 + _chip+1, _channel+1, landauMPV);
  //     // if(landauMPV <300) count_dead++;
  //     // if(landauMPV - HG_pedMean>0) {
  //     //   MIP_woPed->SetBinContent(_layer*6 + _chip+1, _channel+1, landauMPV - HG_pedMean);
  //     //   if(_layer<4 || (_layer>=28&&_layer<30)) {
  //     //     MIP_15->Fill(landauMPV - HG_pedMean);
  //     //     if(landauMPV - HG_pedMean < 300) cout << _layer << "," << _chip << "," << _channel << endl;
  //     //   }
  //     //   if(_layer>=4 && _layer<28) MIP_10->Fill(landauMPV - HG_pedMean);
  //     // }
  //     // // else cout << _layer << "," << _chip << "," << _channel << endl;

  //     // MIP_diff_1d->Fill(landauMPV_normal - landauMPV);
  //     // if(landauMPV_normal - landauMPV > 250 || landauMPV_normal - landauMPV < -250){
  //     //   count_out_of_range++;
  //     //   // continue;
  //     // }
  //     // if(landauMPV_normal - landauMPV > 50 && landauMPV_normal - landauMPV < -50){
  //     //   count_in_range++;
  //     // }
  //     // MIP_diff->SetBinContent(_layer*6 + _chip+1, _channel+1, landauMPV_normal - landauMPV);

  //     // // MIP_diff_1d->Fill(landauMPV_normal - landauMPV_auto);
  //     // // if(landauMPV_normal - landauMPV_auto > 250 || landauMPV_normal - landauMPV_auto < -250){
  //     // //   count_out_of_range++;
  //     // //   continue;
  //     // // }
  //     // // if(landauMPV_normal - landauMPV_auto > 50 && landauMPV_normal - landauMPV_auto < -50){
  //     // //   count_in_range++;
  //     // // }
  //     // // MIP_diff->SetBinContent(_layer*6 + _chip+1, _channel+1, landauMPV_normal - landauMPV_auto);

	// }
  // MIP_diff->SetStats(kFALSE);
  // MIP_diff->SetXTitle("layer*6+chip");
  // MIP_diff->SetYTitle("channel");
	// MIP_diff->Draw("colz");
  // MIP->SetStats(kFALSE);
  // MIP->SetXTitle("layer*6+chip");
  // MIP->SetYTitle("channel");
  // MIP_woPed->SetStats(kFALSE);
  // MIP_woPed->SetXTitle("layer*6+chip");
  // MIP_woPed->SetYTitle("channel");
  // MIP_woPed->Draw("colz");
  // MIP_10->SetXTitle("MIP [ADC]");
  // MIP_15->SetXTitle("MIP [ADC]");
  // MIP_15->Draw();
  // MIP_entry_15->SetMarkerStyle(4);
  // cout << count_out_of_range << "," << count_in_range << "," << count_dead << endl;

}