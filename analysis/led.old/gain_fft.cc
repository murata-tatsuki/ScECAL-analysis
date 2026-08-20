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
#include "TMultiGraph.h"
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
#include "TLegend.h"

#include <chrono>
#include <time.h>
#include "TStyle.h"

// #include "EBUdecode.h"
// #include "EBUdecode.cxx"
// #include "langaus.C"
#include "fft_fullrange_tempcorr_plusNphoto_v5.C"

using namespace std;


int* position_itr(int _layer, double _x, double _y, int canvas_division){		// returns a itr of the channel with x,y coordinate(0~4, 0~41) and the canvas dividion number and the canvas position
	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	const double _zInitialPosition = 1.;
	const double _zIntervalEven = 12.7; // mm
	const double _zIntervalOdd = 19.9; // mm
	double _xMin = -rowNu*_xInterval/2.;
	double _xMax = rowNu*_xInterval/2.;
	double _yMin = -columnNu*_yInterval/2.;
	double _yMax = columnNu*_yInterval/2.;
	double _zMin = -10; double _zMax = 320; int _zBin = 165;
	double coormin = -(rowNu-1)*_xInterval/2.;
	double accurateBinY[10] = {-113.3,-68.3, -67.9,-22.9, -22.5,22.5, 22.9,67.9, 68.3,113.3};
	double accurateBinX[84] = {0};
	
	double xNu_ = _layer%2==0 ? columnNu : rowNu;
	double yNu_ = _layer%2==0 ? rowNu : columnNu;
	double xMax_ = _layer%2==0 ? _yMax : _xMax;
	double yMax_ = _layer%2==0 ? _xMax : _yMax;
	double xInterval_ = _layer%2==0 ? _yInterval : _xInterval;
	double yInterval_ = _layer%2==0 ? _xInterval : _yInterval;

	int x_itr = (_x + xMax_) / xInterval_ + 0.2;
	int y_itr = (_y + yMax_) / yInterval_ + 0.2;


	static int __position[4] = { 0 };
	__position[0] = x_itr;
	__position[1] = y_itr;
	__position[2] = _layer%2==0 ? ((41-y_itr) * columnNu + x_itr) / (210/canvas_division)  :  x_itr/(rowNu/canvas_division);
	__position[3] = _layer%2==0 ?  (41-y_itr) * columnNu + x_itr + 1 - __position[2]*(210/canvas_division)  :  (4-y_itr)*rowNu/canvas_division + x_itr + 1 - __position[2]*(rowNu/canvas_division);



	return __position;
}


int main(int argc, char* argv[])
{ 
		if(argc < 3){																										 //エラー処理
				cout << "usage: ./compiled output.root decode_1.root	...	decode_n.root" << endl;
				return 1;
		}


		//cout << argv[1] << endl;
	cout << "=====>	" << argv[1] << endl;
	
	int rawfilenum = (argc - 2);

	TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	int total_entry_max=0;

	// for(int i=0; i<rawfilenum; i++){
	// 	filein[i] = new TFile(argv[i+2]);
	// 	tree[i] = (TTree*) filein[i]->Get("Raw_Hit");
	// 	entry_max[i] = tree[i]->GetEntries();
	// 	total_entry_max += entry_max[i];
	// }



	
	int Run_Num, Event_Time, Event_Num, DetectorID;
	vector<int> *CellID = nullptr;
	vector<int> *HitTag = nullptr;
	vector<vector<double>> *Temperature = nullptr;
	vector<double> *HG_Charge = nullptr;
	vector<double> *LG_Charge = nullptr;
	


	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	const double _zInitialPosition = 1.;
	const double _zIntervalEven = 12.7; // mm
	const double _zIntervalOdd = 19.9; // mm
	double _xMin = -rowNu*_xInterval/2.;
	double _xMax = rowNu*_xInterval/2.;
	double _yMin = -columnNu*_yInterval/2.;
	double _yMax = columnNu*_yInterval/2.;
	double _zMin = -10; double _zMax = 320; int _zBin = 165;
	double coormin = -(rowNu-1)*_xInterval/2.;
	double accurateBinY[10] = {-113.3,-68.3, -67.9,-22.9, -22.5,22.5, 22.9,67.9, 68.3,113.3};
	double accurateBinX[84] = {0};
	






	TFile fileout(argv[1],"RECREATE");




	int histo_energy_max = 500;
	int histo_energy_min = 0;
	int histo_nhit_max = 300;
	int histo_nhit_min = 0;
	



	// 保存形式

	// TH2F *gainMap[30];			// channelごとのgain
	// TH2F *mipMap[30];			// channelごとのMIP
	TH2F *hit_2d[30];					// channelごとのhit数
	for(int i=0;i<30;i++){
		double xNu_ = i%2==0 ? columnNu : rowNu;
		double yNu_ = i%2==0 ? rowNu : columnNu;
		double xMax_ = i%2==0 ? _yMax : _xMax;
		double yMax_ = i%2==0 ? _xMax : _yMax;
		
		// gainMap[i] = new TH2F(Form("gainMap_%d",i),Form("gain at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
		// mipMap[i] = new TH2F(Form("mipMap_%d",i),Form("mip at layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
		hit_2d[i] = new TH2F(Form("hit_2d_%d",i),Form("hit map layer %d",i), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
	}


	char char_tmp[100];
	TH1F* adc[layerNu][chipNu][channelNu];
	TH1F* adc_1hit[layerNu][chipNu][channelNu];
	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
				adc[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 1600, 400, 2000);
				adc[i_layer][i_chip][i_channel]->SetTitleSize(0.16);
				adc[i_layer][i_chip][i_channel]->SetStats(0);
				adc[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
				adc[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
				sprintf(char_tmp, "Layer%d_Chip%d_Chn%d_1hit", i_layer, i_chip, i_channel);
				adc_1hit[i_layer][i_chip][i_channel] = new TH1F(char_tmp, char_tmp, 1600, 400, 2000);
				adc_1hit[i_layer][i_chip][i_channel]->SetTitleSize(0.16);
				adc_1hit[i_layer][i_chip][i_channel]->SetStats(0);
				adc_1hit[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.08);
				adc_1hit[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.08);
			}
		}
	}


	cout << "made save TTree" << endl;


	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;

		// for(int i=0; i<rawfilenum; i++){
			filein[irawfile] = new TFile(argv[irawfile+2]);
			tree[irawfile] = (TTree*) filein[irawfile]->Get("Raw_Hit");
			entry_max[irawfile] = tree[irawfile]->GetEntries();
			total_entry_max += entry_max[irawfile];
		// }
	
		// tree[irawfile]->SetBranchAddress("Run_Num", &Run_Num);
		// tree[irawfile]->SetBranchAddress("Event_Time", &Event_Time);
		tree[irawfile]->SetBranchAddress("CellID", &CellID);
		tree[irawfile]->SetBranchAddress("HitTag", &HitTag);
		tree[irawfile]->SetBranchAddress("HG_Charge", &HG_Charge);
		// tree[irawfile]->SetBranchAddress("LG_Charge", &LG_Charge);
		// tree[irawfile]->SetBranchAddress("Temperature", &Temperature);

		for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
			if(ientry%100000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);

			int hitLayers[32] = {0};

			for(int i=0; i<(int)CellID->size(); i++){
				int _layer = CellID->at(i)/1e5;
				if(HitTag->at(i)==1) hitLayers[_layer]++;
			}

			for(int i=0; i<(int)CellID->size(); i++){
				int _layer = CellID->at(i)/1e5;
				int _chip = (CellID->at(i)/10000)%10;
				int _channel = CellID->at(i)%100;
				int _hitTags = HitTag->at(i);

				if(_hitTags==0) continue;

				if(hitLayers[_layer]==1) adc_1hit[_layer][_chip][_channel]->Fill(HG_Charge->at(i));
				adc[_layer][_chip][_channel]->Fill(HG_Charge->at(i));
			}

		}
		Run_Num=0, Event_Time=0, Event_Num=0, DetectorID=0;
		CellID = nullptr;
		HG_Charge = nullptr;
		LG_Charge = nullptr;
		HitTag = nullptr;
		Temperature = nullptr;
	
		filein[irawfile]->Close();


	}


	Double_t langaus_peakP, langaus_peakPError;
	TF1 *langaus_func;

	double _landauMPV, _landauWidth, _gausSigma, _landauMPV_tempCorrected, _landauWidth_tempCorrected, _gausSigma_tempCorrected, _ChiSqr, _ChiSqr_tempCorrected;
  int _chn_Entries;
  int _cellID_fit;
  int _Ndf, _Ndf_tempCorrected;  
  // TH1F *_MIP_histo = new TH1F();
  // TH1F *_MIP_histo_tempCorrected = new TH1F();
  // TGraphErrors *_temp_graph = new TGraphErrors();

	TTree* mip_tree = new TTree("T_Event", "after calibration");
  mip_tree->Branch("chnEntries", &_chn_Entries);
  mip_tree->Branch("cellIDs", &_cellID_fit);
  mip_tree->Branch("landauMPV", &_landauMPV);
  mip_tree->Branch("landauWidth", &_landauWidth);
  mip_tree->Branch("gausSigma", &_gausSigma);
  mip_tree->Branch("ChiSqr", &_ChiSqr);
  mip_tree->Branch("ndf", &_Ndf);
  // mip_tree->Branch("histo", &_MIP_histo);
  // mip_tree->Branch("landauWidth_tempCorrected", &_landauWidth_tempCorrected);
  // mip_tree->Branch("gausSigma_tempCorrected", &_gausSigma_tempCorrected);
  // mip_tree->Branch("landauMPV_tempCorrected", &_landauMPV_tempCorrected);
  // mip_tree->Branch("ChiSqr_tempCorrected", &_ChiSqr_tempCorrected);
  // mip_tree->Branch("ndf_tempCorrected", &_Ndf_tempCorrected);
  // mip_tree->Branch("histo_tempCorrected", &_MIP_histo_tempCorrected);
  // mip_tree->Branch("tempGraph", &_temp_graph);


	

	cout << "saving raw histograms" << endl;
	fileout.mkdir("raw_histogram");
	fileout.cd("raw_histogram");
	for(int i=0; i<layerNu;i++){
		cout << i << "/" << layerNu << endl;
		for(int j=0; j<chipNu;j++){
			for(int k=0; k<channelNu;k++){
			// 	if (j == 5 && k > 29) continue;
      //   _chn_Entries = adc_1hit[i][j][k]->GetEntries();
      //   _cellID_fit = i * 1e5 + j * 1e4 + k;

			// 	/////////////////////////// langaus setting 
      //   double fr[2];
      //   double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      //   fr[0]=0.3*adc_1hit[i][j][k]->GetMean();
      //   fr[1]=3.0*adc_1hit[i][j][k]->GetMean();
      //   pllo[0]=0.5; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
      //   plhi[0]=5.0; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=5.0;
      //   sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;
      //   double chisqr;
      //   int    ndf;
      //   sv[1] = i<4 || i>28 ? 45 : 20;
      //   sv[3] = i<4 || i>28 ? 100 : 30;
      //   pllo[0] = 0;
      //   pllo[1] = i<4 || i>28 ? 400+200 : 400+100;
      //   pllo[3] = 0;
      //   plhi[0] = 150;
      //   plhi[1] = i<4 || i>28 ? 400+1000 : 400+600;
      //   plhi[3] = 250;
      //   ///////////////////////////
      //   double SNRPeak, SNRFWHM;
      //   TF1 *fitsnr = langaufit(adc_1hit[i][j][k],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
      //   langaupro(fp,SNRPeak,SNRFWHM);
      //   // Global style settings
      //   gStyle->SetOptStat(1111);
      //   gStyle->SetOptFit(111);
      //   gStyle->SetLabelSize(0.03,"x");
      //   gStyle->SetLabelSize(0.03,"y");
      //   fitsnr->Draw("same");

      //   // if(i_layer<4 || i_layer>28) langaus(adc_1hit[i][j][k],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+200, HG_ped[i_layer][i_chip][i_channel]+1000);
      //   // if(i_layer>=4 && i_layer<=28) langaus(adc_1hit[i][j][k],&langaus_func,&langaus_peakP,&langaus_peakPError, 50, 400);
      //   // if(i_layer>=4 && i_layer<=28) langaus(adc_1hit[i][j][k],&langaus_func,&langaus_peakP,&langaus_peakPError, 400, 800);
      //   // if(i_layer>=4 && i_layer<=28) langaus(adc_1hit[i][j][k],&langaus_func,&langaus_peakP,&langaus_peakPError, HG_ped[i_layer][i_chip][i_channel]+50, HG_ped[i_layer][i_chip][i_channel]+400);

      //   // _landauWidth = langaus_peakPError;
      //   // _landauMPV = langaus_peakP;
      //   //_gausSigma = fun->GetParameter(3);
      //   // _ChiSqr = langaus_func->GetChisquare();
      //   // _Ndf = langaus_func->GetNDF();
      //   _landauWidth = fp[0];
      //   _landauMPV = fp[1];
      //   _gausSigma = fp[3];
      //   _ChiSqr = chisqr;
      //   _Ndf = ndf;
      //   // cout<<"--------  "<<_landauMPV<<" , "<<_landauWidth<<" , "<<_gausSigma<<endl;
      //   // _MIP_histo = adc_1hit[i][j][k];
      //   // _temp_graph = g_temperature[i_layer][i_chip][i_channel];
			// 	mip_tree->Fill();

				adc[i][j][k]->Write(Form("Layer%d_Chip%d_Chn%d",i,j,k));
				adc_1hit[i][j][k]->Write(Form("Layer%d_Chip%d_Chn%d_1hit",i,j,k));
			}
		}
	}
	fileout.cd();
	// mip_tree->Write();
	cout << "saved raw histograms" << endl;



	int canvas_division = 7;
	int c_num = 42/canvas_division;
	TCanvas *C_adc[layerNu][canvas_division];
	TCanvas *C_adc_1hit[layerNu][canvas_division];
	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		

		int xNu_ = i_layer%2==0 ? columnNu : rowNu;
		int yNu_ = i_layer%2==0 ? rowNu : columnNu;

		int xNum_per_1canvas = i_layer%2==0 ? columnNu : rowNu/canvas_division;
		int yNum_per_1canvas = i_layer%2==0 ? rowNu/canvas_division : columnNu;

		for (int i_div = 0; i_div < canvas_division; ++i_div) {
			C_adc[i_layer][i_div] = new TCanvas(Form("layer%d_%d_adc", i_layer, i_div), Form("layer%d_%d_adc", i_layer, i_div), 1);
			C_adc[i_layer][i_div]->Divide(xNum_per_1canvas, yNum_per_1canvas);

			C_adc_1hit[i_layer][i_div] = new TCanvas(Form("layer%d_%d_adc_1hit", i_layer, i_div), Form("layer%d_%d_adc_1hit", i_layer, i_div), 1);
			C_adc_1hit[i_layer][i_div]->Divide(xNum_per_1canvas, yNum_per_1canvas);
		}
	}
	

	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				if(i_chip==5 && i_channel>29) continue;

				double *_position = EBUdecode(i_layer,i_chip,i_channel);
				int *_position_itr = position_itr(i_layer, _position[0], _position[1], canvas_division);

				// if(i_layer<2) cout << _position[0] << ", " << _position[1] << ", " << _position[2] << "    " << _position_itr[0] << ", " << _position_itr[1] << ", " << _position_itr[2] << ", " << _position_itr[3] << endl;
				

				C_adc[i_layer][_position_itr[2]]->cd();
				C_adc[i_layer][_position_itr[2]]->cd(_position_itr[3]);
				// gStyle->SetTitleAlign(33);
				// gStyle->SetTitle(.99);
				gStyle->SetTitleFontSize(0.1);
				adc[i_layer][i_chip][i_channel]->Draw();

				C_adc_1hit[i_layer][_position_itr[2]]->cd();
				C_adc_1hit[i_layer][_position_itr[2]]->cd(_position_itr[3]);
				gStyle->SetTitleFontSize(0.1);
				adc_1hit[i_layer][i_chip][i_channel]->Draw();
				
			}
		}
	}

	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		for (int i_div = 0; i_div < canvas_division; ++i_div) {
			C_adc[i_layer][i_div]->Write();
			C_adc[i_layer][i_div]->SaveAs(Form("pic/layer%d_%d_adc.png", i_layer, i_div));
			C_adc_1hit[i_layer][i_div]->Write();
			C_adc_1hit[i_layer][i_div]->SaveAs(Form("pic/layer%d_%d_adc_1hit.png", i_layer, i_div));
		}
	}





}
