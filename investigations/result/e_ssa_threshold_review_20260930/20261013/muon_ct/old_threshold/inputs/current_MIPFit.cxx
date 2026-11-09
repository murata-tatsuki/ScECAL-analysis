// This code package is for 
#include <fstream>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TMath.h>
#include "TH1.h"
#include "TF1.h"
#include "TVirtualPad.h"
#include "TPad.h"
//#include "SimpleCalibration.h"
#include "MIPFit.h"
Double_t langaufun(Double_t* x, Double_t* par)
{
	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation),
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
	Double_t mpshift = -0.22278298;			 // Landau maximum location

	// Control constants
	Double_t np = 100.0; // number of convolution steps
	Double_t sc = 5.0;	 // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow, xupp;
	Double_t step;
	Double_t i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0];

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp - xlow) / np;

	// Convolution integral of Landau and Gaussian by sum
	for (i = 1.0; i <= np / 2; i++) {
		xx = xlow + (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);

		xx = xupp - (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
	}

	return (par[2] * step * sum * invsq2pi / par[3]);
}

using namespace std;

void MIPCalibration::ReadTree(TTree* tree)
{
	tree->SetBranchAddress("Run_Num",&_run_no);
	tree->SetBranchAddress("Event_Time",&_event_time);
	tree->SetBranchAddress("CycleID", &_cycleID);
	tree->SetBranchAddress("TriggerID", &_triggerID);
	tree->SetBranchAddress("CellID", &_cellID);
	tree->SetBranchAddress("BCID", &_bcid);
	tree->SetBranchAddress("HitTag", &_hitTag);
	tree->SetBranchAddress("GainTag", &_gainTag);
	tree->SetBranchAddress("HG_Charge", &_chargehi);//for high gain charge
	tree->SetBranchAddress("LG_Charge", &_chargelo);//for low gain charge 
	tree->SetBranchAddress("Hit_Time", &_hit_time);//for AutoGain mode time information
	tree->SetBranchAddress("Temperature", &_temp);//ECAL temperature data is with scientific data
}

void MIPCalibration::SetTree(TTree* tree)
{
	tree->Branch("ChnEntries", &_chn_entries);
	tree->Branch("CellID", &_cellID_mip);
	tree->Branch("LandauMPV", &_landauMPV);
	tree->Branch("LandauWidth", &_landauWidth);
	tree->Branch("GauSigma", &_gausSigma);
	tree->Branch("ChiSquare", &_ChiSqr);
	tree->Branch("NDF", &_Ndf);
}

void MIPCalibration::PedestalExtract(string pedestal_name)
{
    cout<<"=========== PedestalExtract ============"<<endl;
    // initial the branch of pedestal fitting
    _pedCellID=0;
    _pedMeanCharges=0;
    _pedSigmaCharges=0;

    // *********************************************************************************
    //          extract pedestal fit parameters
    // *********************************************************************************
    fitfile = new TFile(pedestal_name.c_str());
    fitfile->GetObject("ChnLevel",fNtuple1); // channel level fit parameter
    //fitfile->GetObject("ChipLevel",fNtuple1); // chip level fit parameter
    //fitfile->GetObject("CellLevel",fNtuple1); // cell level fit parameter
    if(!fitfile){ cout<<" Get pedestal fit file failed !!!"<<endl;return;}
    fNtuple1 ->SetBranchAddress("cellIDs",&_pedCellID);
    fNtuple1 ->SetBranchAddress("pedMeanCharges",&_pedMeanCharges);
    fNtuple1 ->SetBranchAddress("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple1 ->SetBranchAddress("pedMeanTimes",&_pedMeanTimes);
    fNtuple1 ->SetBranchAddress("pedSigmaTimes",&_pedSigmaTimes);

    for(int entry=0; entry!=fNtuple1->GetEntries(); ++entry)
    {
        fNtuple1->GetEntry(entry);
        //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(unsigned ihit=0; ihit!=_pedCellID->size(); ++ihit)
        {
            int layerID = _pedCellID->at(ihit)/1e5;
            int chipID = _pedCellID->at(ihit)%100000/1e4;
          //  int memoryCells = _pedCellID->at(ihit)%10000/1e2;
            int channelID = _pedCellID->at(ihit)%100;
            pedMeanChargesChnLevel[layerID][chipID][channelID] = _pedMeanCharges->at(ihit);
            pedSigmaChargesChnLevel[layerID][chipID][channelID] = _pedSigmaCharges->at(ihit);
            pedMeanTimesChnLevel[layerID][chipID][channelID] = _pedMeanTimes->at(ihit);
            pedSigmaTimesChnLevel[layerID][chipID][channelID] = _pedSigmaTimes->at(ihit);
        }
    }
}

int MIPCalibration::MIPFit(string input_name, string hist_name, string mip_name)
{
	std::cout << "********** MIPFit now **********" << std::endl;
	string txtname(input_name);
	txtname.erase(0, txtname.find_last_of("/") + 1);
	string histname(hist_name);
	histname.insert(histname.length(), txtname);
	string mipname(mip_name);
	mipname.insert(mipname.length(), txtname);
	// *********************************************************************************
	cout << "MIP fit " << input_name << " ..." << endl;
	TFile* infile = TFile::Open(input_name.c_str(), "READ");
	if(!infile)	cout << "read file failed" << endl;
	TTree* mipTreeRead = (TTree*)infile->Get("Raw_Hit");
	if(!mipTreeRead)	cout << "Get tree failed" << endl;
	ReadTree(mipTreeRead);
	
	TFile* mipfile = new TFile(mipname.c_str(), "RECREATE");
	TTree* miptree = new TTree("MIP_Fit", "MIP_Calibration");
	SetTree(miptree);
	
	char char_tmp[100];
	char h_mip_chip_name[100];
	TH1D* h_MIP[layerNu][chipNu][channelNu];
	TH1D* h_MIP_Chip[layerNu][chipNu];
	
	// Set Format
	for (int i_layer = 0; i_layer < layerNu; ++i_layer)
	{
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) 
		{
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) 
			{
				sprintf(char_tmp, "Layer%d_Chip%d_ Chn%d", i_layer, i_chip, i_channel);
			//	if (3 < i_layer && i_layer < 28)
			//		h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 150, -50, 1000);
			//	else
			//		h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 150, -100, 3000);
				h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 150, -100, 3000);
				h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
				h_MIP[i_layer][i_chip][i_channel]->SetLineWidth(2);
				h_MIP[i_layer][i_chip][i_channel]->SetLineColor(4);
			}
			sprintf(h_mip_chip_name, "Layer%d_Chip%d", i_layer, i_chip);
			if(i_layer>=4 || i_layer <=27)
			{
				h_MIP_Chip[i_layer][i_chip] = new TH1D(h_mip_chip_name, h_mip_chip_name, 800, -100, 3100);
			}
			else
			{
				h_MIP_Chip[i_layer][i_chip] = new TH1D(h_mip_chip_name, h_mip_chip_name, 400, -100, 3100);
			}
			h_MIP_Chip[i_layer][i_chip]->GetXaxis()->SetTitle("ADC");
			h_MIP_Chip[i_layer][i_chip]->SetLineWidth(2);
			h_MIP_Chip[i_layer][i_chip]->SetLineColor(4);
		}
	}
	
	//Read Data
	int totalEntries = mipTreeRead->GetEntries();
	for (int entry = 0; entry < totalEntries; ++entry) 
	{
		mipTreeRead->GetEntry(entry);
		if(entry%2000==0)	cout << " [ Process : MIP Fit : ]  event : " << entry << endl;
		for (unsigned i_hit = 0; i_hit != _cellID->size(); ++i_hit)
		{
			if(_hitTag->at(i_hit)==0)	continue;
			int layerID = _cellID->at(i_hit) / 1e5;
			int chipID = _cellID->at(i_hit) % 100000 / 1e4;
			int channelID = _cellID->at(i_hit) % 100;
			//cout << "layer : chip : channel" << layerID << " " << chipID << " " << channelID << endl; 
			//double charge = _gainTag->at(i_hit)==1 ?  _chargehi->at(i_hit) : _chargelo->at(i_hit); 
			//double time = _hit_time->at(i_hit);	
			double time =  _chargehi->at(i_hit);
			//double charge = _chargelo->at(i_hit); 

			// incident angle correction
			// time = 2. / _realLength * time;

			// temperature dependence correction
			//time = (1+(_temp_mip->at(i_hit)-20)*0.03)*time;
			// if (layerID < 3 || layerID > 28)
			// 	time = (1 + (_temp_mip->at(i_hit) - 20) * 0.0122) * time;
			// else
			// 	time = (1 + (_temp_mip->at(i_hit) - 20) * 0.0073) * time;

			// MIP response decline correction
			//time = time/(1-0.17*(double)entry/(double)totalEntries);
			//int period1 = 20*13000;
			//int period2 = 75*13000;
			//if(entry<=period1) time = time/(1-0.1*(double)entry/(double)period1); // MIP response decline correction
			//else if(entry<=period2) time = time/(1-0.077*(double)(entry-period1)/(double)(period2-period1))*14./12.8; // MIP response decline correction
			//else  time = time/(1-0.017*(double)(entry-period2)/(double)(totalEntries-period2))*14/11.9; // MIP response decline correction
			h_MIP[layerID][chipID][channelID]->Fill(time);
			h_MIP_Chip[layerID][chipID]->Fill(time);
		}
	}
	infile->Close();
	
	int _blindChn = 0;
	TH1D* hBlind = new TH1D("hBlind", "", 32, 0, 32);
	TH1D* hChn = new TH1D("hChn", "", 3200000, 0, 3200000);
	
	TFile* histfile = new TFile(histname.c_str(), "RECREATE");
	TCanvas* C_MIP[layerNu][chipNu];
	TCanvas* C_MIP_chip[layerNu];
	
	for(int i_layer = 0; i_layer < layerNu; ++i_layer)
	{
		cout << " [ LandauGaus Fitting : ] Layer : " << i_layer << endl;
		for(int i_chip = 0; i_chip < chipNu; ++i_chip) 
		{
			histfile->cd();
			sprintf(char_tmp, "layer%d/chip%d", i_layer, i_chip);
			histfile->mkdir(char_tmp);
			histfile->cd(char_tmp);
			
			for (int i_channel = 0; i_channel < channelNu; ++i_channel)
			{
			//	_cellID_mip = i_layer * 1e5 + i_chip * 1e4 + i_channel;
				if(i_chip == 5 && i_channel > 29)	continue;
				if((h_MIP[i_layer][i_chip][i_channel]->GetEntries()) == 0)
				{
					_blindChn++;
					hBlind->Fill(i_layer);
					hChn->Fill(_cellID_mip);
					continue;
				}
			//	_chn_entries = h_MIP[i_layer][i_chip][i_channel]->GetEntries();
			/*
				if (3 < i_layer && i_layer < 28)
				{
					TF1* fun = new TF1("langau", langaufun, 50, 1000, 4);
					fun->SetParLimits(0, 10, 50);
					fun->SetParLimits(1, 30, 300);
					fun->SetParLimits(3, 5, 100);
					fun->SetParameter(0, 20);
					fun->SetParameter(1, 200);
					fun->SetParameter(2, 1000);
					fun->SetParameter(3, 100);

					double initPeak = 120;
					double initLX = initPeak - 40;
					double initRX = initPeak + 200;
					//double initLX = 100;
					//double initRX = 500;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					
					double peakX = fun->GetMaximumX(initLX, initRX);
					double peakY = fun->GetMaximum();
					double rangeLX = peakX-40;
					double rangeRX = peakX+200;
					//rangeLX = fun->GetX(peakY/2.,rangeLX,peakX);
					//rangeRX = fun->GetX(peakY/4.5,peakX,rangeRX);
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "Q", "", rangeLX, rangeRX);
					
					peakX = fun->GetMaximumX(initLX, initRX);
					peakY = fun->GetMaximum();
					//rangeLX = peakX-40;
					//rangeRX = peakX+200;
					_landauWidth = fun->GetParameter(0);
					_gausSigma = fun->GetParameter(3);
					rangeLX = peakX - sqrt(pow(_landauWidth, 2) + pow(_gausSigma, 2));
					rangeRX = peakX + 5 * sqrt(pow(_landauWidth, 2) + pow(_gausSigma, 2));
					//rangeLX = fun->GetX(peakY/2.,rangeLX,peakX);
					//rangeRX = fun->GetX(peakY/4.5,peakX,rangeRX);
					//cout<<i_channel<<" : "<<initPeak<<" , "<<peakX<<", "<<rangeLX<<", "<<rangeRX<<endl;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "Q", "", rangeLX, rangeRX);

					_landauWidth = fun->GetParameter(0);
					_landauMPV = fun->GetParameter(1);
					_gausSigma = fun->GetParameter(3);
					_ChiSqr = fun->GetChisquare();
					_Ndf = fun->GetNDF();
					//cout<<"--------  "<<_landauMPV<<" , "<<_landauWidth<<" , "<<_gausSigma<<endl;
				}
				else
				{
					TF1* fun = new TF1("langau", langaufun, 50, 2500, 4);
					fun->SetParLimits(0, 10, 100);
					fun->SetParLimits(1, 100, 1500);
					fun->SetParLimits(3, 10, 500);

					double initPeak = h_MIP[i_layer][i_chip][i_channel]->GetBinCenter(h_MIP[i_layer][i_chip][i_channel]->GetMaximumBin());
					double initLX = (initPeak - 300) > 200 ? (initPeak - 300) : 200;
					double initRX = initPeak + 1100;
					fun->SetParameter(0, 50);
					fun->SetParameter(1, 800);
					fun->SetParameter(2, 10000);
					fun->SetParameter(3, 160);
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					double peakX = fun->GetMaximumX(initLX, initRX);
					double peakY = fun->GetMaximum();
					double rangeLX = (peakX - 300) > 150 ? (peakX - 300) : 150;
					double rangeRX = peakX + 1000;
					//rangeLX = fun->GetX(peakY/2.,rangeLX,peakX);
					//rangeRX = fun->GetX(peakY/4.5,peakX,rangeRX);
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", rangeLX, rangeRX);
					_landauWidth = fun->GetParameter(0);
					_landauMPV = fun->GetParameter(1);
					_gausSigma = fun->GetParameter(3);
					_ChiSqr = fun->GetChisquare();
					_Ndf = fun->GetNDF();
				}

			*/

				//if(_landauMPV>4000 || _landauMPV <0)	{cout<<"WARNING : "<<i_layer<<" : "<<i_chip<<" : "<<i_channel<<" LandauGaus fit failed"<<endl;continue;}
				miptree->Fill();
				h_MIP[i_layer][i_chip][i_channel]->Write();
			}
			// Chip Level Fit
			_cellID_mip = i_layer * 10 + i_chip ;
			if((i_layer>=4 && i_layer <=27) || i_layer == 30 || i_layer==31)
			{
				TF1* fun = new TF1("langau", langaufun, 50, 2500, 4);
				//fun->SetParLimits(0, 10, 100);
				fun->SetParLimits(0, 10, 50);
				fun->SetParLimits(1, 300, 800);
				fun->SetParLimits(3, 0, 10);
				//fun->SetParameter(0, 50);
				//fun->SetParameter(1, 800);
				fun->SetParameter(0, 20);
				fun->SetParameter(1, 500);
				fun->SetParameter(2, 10000);
				fun->SetParameter(3, 2);

				double initPeak = h_MIP_Chip[i_layer][i_chip]->GetBinCenter(h_MIP_Chip[i_layer][i_chip]->GetMaximumBin());
				//double initLX = (initPeak - 300) > 200 ? (initPeak - 300) : 200;
				//double initRX = initPeak + 1100;
				double initLX = (initPeak - 300) > 200 ? (initPeak - 300) : 200;
				double initRX = initPeak + 1100;
				h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", initLX, initRX);
				h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", initLX, initRX);
				double peakX = fun->GetMaximumX(initLX, initRX);
				double peakY = fun->GetMaximum();
				double rangeLX = (peakX - 300) > 150 ? (peakX - 300) : 150;
				double rangeRX = peakX + 1000;
				//rangeLX = fun->GetX(peakY/2.,rangeLX,peakX);
				//rangeRX = fun->GetX(peakY/4.5,peakX,rangeRX);
				//h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", rangeLX, rangeRX);
				h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", rangeLX , rangeRX);
				_landauWidth = fun->GetParameter(0);
				_landauMPV = fun->GetParameter(1);
				_gausSigma = fun->GetParameter(3);
				_ChiSqr = fun->GetChisquare();
				_Ndf = fun->GetNDF();
			}
			else
			{
				TF1* fun = new TF1("langau", langaufun, 50, 2500, 4);
				//fun->SetParLimits(0, 10, 100);
				//fun->SetParameter(0, 50);
				//fun->SetParameter(1, 800);
				fun->SetParameter(0, 55);
				fun->SetParameter(1, 799);
				fun->SetParameter(2, 10000);
				fun->SetParameter(3, 13);
				
				fun->SetParLimits(0, 35, 75);
				fun->SetParLimits(1, 500, 1300);
				fun->SetParLimits(2, 5000, 30000);
				fun->SetParLimits(3, 0, 30);

				double initPeak = h_MIP_Chip[i_layer][i_chip]->GetBinCenter(h_MIP_Chip[i_layer][i_chip]->GetMaximumBin());
				//double initLX = (initPeak - 300) > 200 ? (initPeak - 300) : 200;
				//double initRX = initPeak + 1100;
				//double initLX = initPeak > 650 ? (initPeak-150) : 500;
				double initRX = initPeak + 1400;
				//h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", 600, initRX);
				h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", 600, initRX);
				//double peakX = fun->GetMaximumX(initLX, initRX);
				//double peakY = fun->GetMaximum();
				//double rangeLX = (peakX - 300) > 200 ? (peakX - 300) : 200;
				//double rangeRX = peakX + 1000;
				//rangeLX = fun->GetX(peakY/2.,rangeLX,peakX);
				//rangeRX = fun->GetX(peakY/4.5,peakX,rangeRX);
				h_MIP_Chip[i_layer][i_chip]->Fit(fun, "QR", "", 700,2200);
				_landauWidth = fun->GetParameter(0);
				_landauMPV = fun->GetParameter(1);
				_gausSigma = fun->GetParameter(3);
				_ChiSqr = fun->GetChisquare();
			}
			sprintf(char_tmp, "layer%d", i_layer);
			histfile->cd(char_tmp);
			sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
			C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
			C_MIP[i_layer][i_chip]->Divide(6, 6);
			for (int i_channel = 0; i_channel < channelNu; ++i_channel)
			{
				C_MIP[i_layer][i_chip]->cd(i_channel + 1);
				h_MIP[i_layer][i_chip][i_channel]->Draw();
			}
			C_MIP[i_layer][i_chip]->Write();

		}

		char canv_name[100];
		sprintf(canv_name, "layer%d", i_layer);
		C_MIP_chip[i_layer] = new TCanvas(canv_name, canv_name, 1);
		C_MIP_chip[i_layer]->Divide(3,2);
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
			cout << "layer: " << i_layer << endl;
			if(i_chip<=6) gStyle->SetOptFit(1111);
			C_MIP_chip[i_layer]->cd();
			C_MIP_chip[i_layer]->cd(i_chip + 1);
			//h_MIP_Chip[i_layer][i_chip]->DrawClone();
		/*	if((i_layer>=4 || i_layer <=27) || i_layer == 30 || i_layer==31)
			{
				h_MIP_Chip[i_layer][i_chip]->GetXaxis()->SetRangeUser(300,1600);
			}
			else
			{
				h_MIP_Chip[i_layer][i_chip]->GetXaxis()->SetRangeUser(300,2200);
			}
		*/
			h_MIP_Chip[i_layer][i_chip]->GetXaxis()->SetRangeUser(300,2200);
			h_MIP_Chip[i_layer][i_chip]->Draw();
		}
		C_MIP_chip[i_layer]->Write();
	}
	
	cout << "Total Blind Chn : " << _blindChn << " **************" << endl;
	histfile->cd();
	hBlind->Write();
	hChn->Write();
	mipfile->cd();
	miptree->Write();

	histfile->Close();
	mipfile->Close();
}
