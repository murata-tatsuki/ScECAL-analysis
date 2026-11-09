#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TKey.h>
#include "MCDigi.h"
#include "Extract.h"

using namespace std; 

Double_t langaufun(Double_t* x, Double_t* par) {
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
	Double_t invsq2pi = 0.3989422804014; // (2 /pi)^(-1/2)
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

void MCDigi::InDigiBranchClear() {
	if(!_digiCellID.empty()) 		_digiCellID.clear();
	if(!_digiBCID.empty())	 		_digiBCID.clear();
	if(!_digiHitTag.empty()) 		_digiHitTag.clear();
	if(!_digiGainTag.empty()) 		_digiGainTag.clear();
	if(!_digiHitHighADC.empty()) 	_digiHitHighADC.clear();
	if(!_digiHitLowADC.empty()) 	_digiHitLowADC.clear();
	if(!_digiHitTime.empty()) 		_digiHitTime.clear();
	if(!_digiTemperature.empty()) _digiTemperature.clear();
	//if(!_digiHitEnergy.empty()) 	_digiHitEnergy.clear();
	//if(!_digiHitX.empty()) 			_digiHitX.clear();
	//if(!_digiHitY.empty()) 			_digiHitY.clear();
	//if(!_digiHitZ.empty()) 			_digiHitZ.clear();
}

Point MCDigi::CellID2Position(int cellID) {
	int layerNo = cellID / (1e5);
	int chipNo 	= (cellID%100000) / (1e4);
	int chanNo  = (cellID%100);
	int layerOrderNo = decodeID[chipNo][chanNo];// layerOrderNo: 0-209 
	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	int _xID = layerOrderNo % rowNu; // 0-41
	int _yID = layerOrderNo / rowNu; // 0-4
	static double _position[3] = { 0 };
	double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
	double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;
	if (layerNo % 2 == 0) {
		_position[0] = -y0; 
		_position[1] = -x0;
		_position[2] = 1 + layerNo / 2 * 19.9;
	} else {
		_position[0] = -x0; 
		_position[1] = -y0;
		_position[2] = 12.2 + (layerNo - 1) / 2 * 19.9;
		//_position[2] = 12.95+(layerNo-1)/2*19.9;
	}
	Point point(_position[0], _position[1], _position[2]);
	return point;
}

int MCDigi::CopyNo2CellID(int copyNo) {// copyNo: 0-6719 
	int layerNo, decodeIDNow, chipNo, chanNo;
	layerNo  = copyNo / 210;
	decodeIDNow = copyNo % 210;
	bool isFound = false;
	for(int i_chip = 0 ; i_chip < CHIPNO ; i_chip++) {
		for(int i_chan = 0 ; i_chan < CHANNO; i_chan++) {
			if(i_chip==5 && i_chan>29) continue;
			if(decodeID[i_chip][i_chan]==decodeIDNow) {
				chipNo = i_chip;
				chanNo = i_chan;
				isFound = true;
				break;
			} else	
				continue;
		}
		if(isFound) break;
	}
	return (layerNo*1e5 + chipNo*1e4 + chanNo);
}

void MCDigi::PixelAva(int npx, int npy) {
	//cout << "===== Pixel Ava Start ========" << endl;
	if (isa[npx][npy]) return;
	isa[npx][npy]=1;
	n_fired++;
	if (gRandom->Uniform(0,1) >= pCrosstalk) return;
	
	// pick one valid pixel
	int deltax, deltay;
	bool isIn = false;//suppose not in 
	
	int randomx = 1;
	int randomy = 1;
	while (randomx == 1 && randomy == 1) {// not self
		randomx = gRandom->Integer(3);
		randomy = gRandom->Integer(3);
		//cout << "random : " << randomx << " " << randomy << endl;
	}
	deltax = randomx - 1;
	deltay = randomy - 1;
	//cout << "delta : " << deltax << " " << deltay << endl;
	//cout << "npx : npy : length : " << npx << "	" << npy << " " << length << endl;		
	if (npx>(length-1-deltax) || npx<(-deltax) || npy>(length-1-deltay) || npy<(-deltay)) {
		isIn = false;
		return;
	} else {
		isIn = true;
	}
	
	if (isIn)	PixelAva(npx+deltax,npy+deltay);
	//cout << "===== Pixel Ava End ========" << endl;
	return;
}

int MCDigi::PixelCorrection(int nPhoton, int nPixel, double pde, double pCT) {
    length = sqrt(nPixel);
	int area = length * length;
	pCrosstalk = pCT;	
	
	n_fired = 0;
	isa.resize(length, vector<int>(length, 0));
	for (int i = 0; i < nPhoton; i++) {
		if (pde < 1 && gRandom->Uniform(0,1) > pde) continue;
		//cout << " ========== Photon : " << i << endl;
		bool isPixel = false;
		int npx, npy;
		npx = gRandom->Integer(length);
		npy = gRandom->Integer(length);
		PixelAva(npx,npy);
		//cout << "End Photon : " << i << endl;
		//cout << "nFired : " << n_fired << endl;
	}
	//cout << "=============" << endl;
	isa.clear();
	return n_fired;
}


Result MCDigi::StepDigi(double energyDep, int cellID, Extract& ex)
{
	int layerNo = cellID/100000;
	int chipNo = (cellID%100000) / 10000;
	int chanNo = cellID % 100;
	int SiPMType = ((layerNo < 4) || (layerNo >27)) ? 1 : 0;
	
	// The Poisson distribution of scintillation process is not considered
    //sPhoton = gRandom->Poisson(energyDep * LightYieldperMeV); 10000p.e./MeV
	double MIPs = energyDep / MIPEnergy;
	if (ex.SiPMGain[layerNo][chipNo][chanNo] <= 1) {
		ex.SiPMGain[layerNo][chipNo][chanNo] = (SiPMType == 0) ? 8 : 22.5;
	} 
	// Keep the MIP value accepted/supplemented by the calibration policy.
	if (!std::isfinite(ex.MIPPeakValue[layerNo][chipNo][chanNo]) ||
	    ex.MIPPeakValue[layerNo][chipNo][chanNo] <= 0.) {
		throw std::runtime_error("Invalid adopted MIP factor for CellID " + std::to_string(cellID));
	}
	//cout << "ex.MIPPeakValue : ex.SiPMGain" << "	" << ex.MIPPeakValue[layerNo][chipNo][chanNo] << "	" << ex.SiPMGain[layerNo][chipNo][chanNo] << endl;
	
	// consider about croostalk
	// if (SiPMType == 0 ) {
	// 	ex.MIPPeakValue[layerNo][chipNo][chanNo] /= 1.1;
	// } else {
	// 	ex.MIPPeakValue[layerNo][chipNo][chanNo] /= 1.15;
	// }

	double pCrosstalk10 = 0.1 ; double pCrosstalk15 = 0.15 ;
	double pct_SiPMType = SiPMType == 0 ? pCrosstalk10 : pCrosstalk15;
	
	// double lightYieldperMIP = (ex.MIPPeakValue[layerNo][chipNo][chanNo]) / ex.SiPMGain[layerNo][chipNo][chanNo]; // p.e/MIP
	double lightYieldperMIP = (1-pct_SiPMType) * (ex.MIPPeakValue[layerNo][chipNo][chanNo]) / ex.SiPMGain[layerNo][chipNo][chanNo]; // p.e/MIP
	
	double photonNo = -1.;
	double pde10  = 0.1;   double pde15 = 0.25;
	double lightOutputNonUni = 0.042;
	photonNo = (SiPMType == 0) ? (MIPs*lightYieldperMIP/pde10) : (MIPs*lightYieldperMIP/pde15);
	photonNo = floor(photonNo);
	photonNo = gRandom->Poisson(photonNo);
	photonNo = gRandom->Uniform((1-lightOutputNonUni)*photonNo, (1+lightOutputNonUni)*photonNo);
	photonNo = round(photonNo);
	
	double peNo   = -1.;
	double nSeed  = 0. ;
	int nPixel10 = 10000;  int nPixel15 = 4489;
	// double pCrosstalk10 = 0.1 ; double pCrosstalk15 = 0.15 ;
	//double pCrosstalk10 = 0 ; double pCrosstalk15 = 0;
	//Note : nFired = nPixel * (1 - exp(-nSeed / (nPixel+0.) ));peNo = round(nFired);
	if (SiPMType == 0 ) {
		peNo = PixelCorrection(photonNo, nPixel10, pde10, pCrosstalk10);
		//peNo = round(photonNo*pde10);
		//peNo = round(nPixel10 * (1 - exp(- photonNo*pde10 / (nPixel10+0.) )));
		//photon2pe[layerNo][chipNo][chanNo]->Fill(photonNo,peNo);
	} else {
		peNo = PixelCorrection(photonNo, nPixel15, pde15, pCrosstalk15);
		//peNo = round(photonNo*pde15);
		//peNo = round(nPixel15 * (1 - exp(- photonNo*pde15 / (nPixel15+0.) )));
		//photon2pe[layerNo][chipNo][chanNo]->Fill(photonNo,peNo);
	}
	//double peNo = MIPs * lightYieldperMIP;
	double ADCMean  = peNo * ex.SiPMGain[layerNo][chipNo][chanNo];
	//double SiPMSigma2 = ex.SiPMGainError[layerNo][chipNo][chanNo] / ex.SiPMGain[layerNo][chipNo][chanNo];
	//double ADCSigma   = ADCMean * sqrt(1/peNo + pow(SiPMSigma2,2));
	double ADCSigma	= 3 * sqrt(peNo);
    ADCMean = gRandom->Gaus(ADCMean, ADCSigma);
	//cout << "peNo : ADCMean : ADCSigma " << peNo << " " << ADCMean << " " << ADCSigma << endl;
	
	double ADCEff = ADCMean;
    double highADC, lowADC;
	double highSwitch = ex.HighADCSwitch[layerNo][chipNo][chanNo];
	int lowSwitch  = 3000;
	//HLRatio~30
    highADC = round(gRandom->Gaus(  (ADCEff + ex.pedMeanHigh[layerNo][chipNo][chanNo]), ex.pedSigmaHigh[layerNo][chipNo][chanNo]));
	lowADC = round(gRandom->Gaus(  (ADCEff / ex.HLRatio[layerNo][chipNo][chanNo] + ex.HLIntercept[layerNo][chipNo][chanNo] + ex.pedMeanLow[layerNo][chipNo][chanNo]) , ex.pedSigmaLow[layerNo][chipNo][chanNo]));
    //double highADCtest =  (ADCEff + ex.pedMeanHigh[layerNo][chipNo][chanNo]);
	//double lowADCtest = (ADCEff / ex.HLRatio[layerNo][chipNo][chanNo] + ex.pedMeanLow[layerNo][chipNo][chanNo]);
	//if(energyDep > 0.8 *0.305) {
	//	test11->Fill(highADCtest-highADC);
	//	test12->Fill(lowADCtest-lowADC);
	//	test13->Fill(highADC/highADCtest-1);
	//	test14->Fill(lowADC/lowADCtest-1);
	//	test21->Fill(highADC,highADCtest);
	//	test22->Fill(lowADC,lowADCtest);
	//}
	
	
	if (highADC > highSwitch) {
		highADC = 4000;
	}
	if (lowADC  > lowSwitch ) {
		lowADC  = 3000;
	}
	
	//Result thisResult(highADCtest/1. , lowADCtest/1.);
	Result thisResult(highADC/1. , lowADC/1.);
    return thisResult;
}

bool MCDigi::ApplyThreshold(double adc_raw, double thre, double threSigma, double pedestal) {
    // 毎イベント、閾値自体を Ath を中心とするガウス分布から生成
    double Ath_event = gRandom->Gaus(thre+pedestal, threSigma);

    // 単純な閾値判定
    return (adc_raw >= Ath_event);
}

//void MCDigi::SimuDigi(string inFileName, string outDecodeFileName, string outCalibFileName, string pedFileName, string HLRatioFileName, string MIPFileName, string SiPMGainFileName)
void MCDigi::SimuDigi(string inFileName, string outFileName, Extract& ex)
{
	gRandom->SetSeed(0);
	cout << "========= MC Digi Start ==========" << endl; 
	
	inFile = TFile::Open(TString(inFileName), "READ");
	inTree = (TTree*)inFile->Get("MC_Truth");
	if(!inTree) cout << "No Input Tree!" << endl;
	inTree->SetBranchAddress("CellID",&_readCellID);
    inTree->SetBranchAddress("Hit_Energy",&_readHitEnergy);
    inTree->SetBranchAddress("Hit_X",&_readHitX);
    inTree->SetBranchAddress("Hit_Y",&_readHitY);
    inTree->SetBranchAddress("Hit_Z",&_readHitZ);

    outFile = new TFile(TString(outFileName), "RECREATE");
	if(!outFile) cout << "FILE IS NOT CREATED!" << endl;
    outTree = new TTree("Raw_Hit", "after digitization");
	if(!outTree) cout << "TREE IS NOT CREATED!" << endl;
    outTree->Branch("Run_Num",		&_digiRunNo			);
    outTree->Branch("Event_Time",	&_digiEventTime		);
    outTree->Branch("CycleID",		&_digiCycleID		);
    outTree->Branch("TriggerID",	&_digiTriggerID		);
    outTree->Branch("CellID",		&_digiCellID		);
    outTree->Branch("BCID",			&_digiBCID			);
    outTree->Branch("HitTag", 		&_digiHitTag		);
    outTree->Branch("GainTag", 		&_digiGainTag		);
    outTree->Branch("HG_Charge", 	&_digiHitHighADC	);
    outTree->Branch("LG_Charge", 	&_digiHitLowADC		);
    outTree->Branch("Hit_Time", 	&_digiHitTime		);
    outTree->Branch("Temperature",&_digiTemperature		);
   	
	/*for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
			photon2pe[il][ip][ic] = new TH2D(Form("%d_%d_%d",il,ip,ic),Form("%d_%d_%d",il,ip,ic),20000,0,100000,1000,0,5000);
			}
		}
	}*/
	//TH1D *Edep = new TH1D("Edep","Edep[MeV]",1000,0,5000);
	//TH2D *HitNo2EnergyDep = new TH2D("HitNo2EnergyDep","HitNo2EnergyDep",750,0,3000,1000,0,1000);
	//test11 = new TH1D("higha","higha",1000,-500,500);
	//test12 = new TH1D("lowa","lowa",1000,-500,500);
	//test13 = new TH1D("high","high",100,-1,1);
	//test14 = new TH1D("low","low",100,-1,1);
	//test21 = new TH2D("testhigh","testhigh",4000,0,4000,4000,0,4000);
	//test22 = new TH2D("testlow","testlow",3000,0,3000,3000,0,3000);
	for (size_t i_entry = 0; i_entry < inTree->GetEntries(); i_entry++) {
	//for (size_t i_entry = 0; i_entry < 3; i_entry++) {
		if (i_entry%5000 == 0) cout << "Process : event : " << i_entry << endl;
		//cout << "======= Process : event : " << i_entry << endl;
        inTree->GetEntry(i_entry);
		double EnergyDep = 0;
		int hitChannelNuEventLevel = 0;
		InDigiBranchClear();
		_digiRunNo 		= 0;// default settings
		_digiEventTime 	= i_entry;// default settings
		_digiCycleID 	= i_entry;// default settings
		_digiTriggerID	= i_entry;// default settings
        for (size_t i_vec = 0; i_vec < _readHitEnergy->size(); i_vec++)
        {
			// if (Hit_Time->at(i_vec) - Tof < TimeCut)
			// {
					//cout << "=========== i_vec : " << i_vec << endl; 
					int layerNo = _readCellID->at(i_vec) / 100000;
					if (layerNo>=LAYERNO)	continue;// 
					int chipNo = (_readCellID->at(i_vec)%100000) / 10000;
					int chanNo = _readCellID->at(i_vec) % 100;
					//EnergyDep += _readCellID->at(i_vec) ;
					//hitChannelNuEventLevel++;
					
					
					Result result = StepDigi(_readHitEnergy->at(i_vec), _readCellID->at(i_vec), ex);
					double hADC = result.hADC;
					double lADC = result.lADC;
					// if(layerNo>3 && layerNo<28) {
					// 	if(hADC<ex.thresholdValue_with_pedestal[layerNo][chipNo][chanNo]) continue;	// checking threshold 
					// } else {
					// 	if(hADC<ex.thresholdValue[layerNo][chipNo][chanNo]+ex.pedMeanHigh[layerNo][chipNo][chanNo]) continue;	// checking threshold 
					// }
					// if(hADC<ex.thresholdValue_with_pedestal[layerNo][chipNo][chanNo]) continue;	// checking threshold 
					// hADC -= ex.thresholdValue[layerNo][chipNo][chanNo];	// checking threshold 
					// lADC -= ex.thresholdValue[layerNo][chipNo][chanNo];	// checking threshold 
					if(!ApplyThreshold(hADC, ex.thresholdValue[layerNo][chipNo][chanNo], ex.thresholdSigmaValue[layerNo][chipNo][chanNo], ex.pedMeanHigh[layerNo][chipNo][chanNo])) continue;
					if(hADC<ex.pedMeanHigh[layerNo][chipNo][chanNo]) continue;									// checking pedestal 
					
					
					_digiCellID.push_back(_readCellID->at(i_vec));
					_digiBCID.push_back(-1);// default settings
					_digiHitTag.push_back(1);
					_digiGainTag.push_back(1);
					_digiHitHighADC.push_back(hADC);
					_digiHitLowADC.push_back(lADC);
					_digiHitTime.push_back(-1.);
			// }
		}
		_digiTemperature.assign(30, vector<double>(16, 20));
		//Edep->Fill(EnergyDep);//exclude muon-like
		//HitNo2EnergyDep->Fill(EnergyDep,hitChannelNuEventLevel);
        outTree->Fill();
		//cout << "End : event : " << i_entry << endl;
    }
	cout << "========= MC Digi Finish ==========" << endl; 
	
    outFile->cd();
	//test11->Write();
	//test12->Write();
	//test13->Write();
	//test14->Write();
	//test21->Write();
	//test22->Write();
    outTree->Write();
	
	/*for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
			outFile->mkdir(Form("%d/%d",il,ip));
			outFile->cd(Form("%d/%d",il,ip));
            for (int ic=0; ic<CHANNO; ic++) {
				photon2pe[il][ip][ic]->Write();
				delete photon2pe[il][ip][ic];
			}
		}
	}*/
	//Edep->SetLineWidth(2);
	//Edep->Write();
	//HitNo2EnergyDep->SetOption("colz");
	//HitNo2EnergyDep->GetXaxis()->SetTitle("Energy[MeV]");
	//HitNo2EnergyDep->GetYaxis()->SetTitle("HitNo");
	//HitNo2EnergyDep->Write();
    // inFile->Close();
	// delete test11;
	// delete test12;
	// delete test13;
	// delete test14;
	// delete test21;
	// delete test22;
	delete outTree;
  delete outFile;
}
/*
int MCDigi::CellID2HistNo(int cellID) {
	int layerNo = cellID / 1e5;
	int chipNo  = (cellID%100000) / 1e4;
	int chanNo  = cellID% 100;
	if (chipNo==5 && chanNo > 29) {
		cout << "!!!!! MORE THAN 30 CHANNELS in CHIP5 !!!!! " << endl;
		return -1;
	}
	return (layerNo*210 + chipNo*36 + chanNo);
}

vector<int> deadchannel = {
	910021, 1300008, 2200018,

	30028, 50003, 50004, 140015, 200007, 230032, 240018, 250000, 300008, 410012, 
	430035, 440030, 440034, 500000, 500022, 510035, 520003, 530033, 630019, 630021, 
	630025, 730027, 730031, 730034, 830019, 830023, 830033, 840035, 850019, 900015, 
	930033, 1020002, 1020024, 1100000, 1120034, 1130029, 1220024, 1220026, 1240028, 1250001, 
	1250002, 1330017, 1330019, 1430025, 1440017, 1500000, 1500002, 1530000, 1530030, 1540011, 
	1540035, 1610005, 1620016, 1630019, 1630029, 1630031, 1650010, 1700003, 1710009, 1730027,
	1750021, 1810011, 1810015, 1850010, 1910030, 1910032, 2000017, 2010035, 2020006, 2030033,
	2100000, 2100002, 2110027, 2110032, 2130002, 2230029, 2240009, 2330027, 2330029, 2330031,
	2500014, 2510035, 2530027, 2530034, 2530035, 2610004, 2630031, 2630033, 2700009, 2730032,
	2730034, 2800021, 2800032, 2810025, 2810029, 2820006, 2820026, 2830031, 2840015, 2840023,
	2920023, 2930019, 2950012, 3010023, 3030027, 3050015, 3130018,

	10030, 10031, 10032, 10033, 10034, 10035, 20033, 20034, 20035, 30013, 
	30015, 30017, 30019, 30021, 30023, 30025, 30027, 30028, 30029, 30030, 
	30031, 30032, 30033, 30034, 30035, 40003, 40015, 40017, 40024, 40025,
	40026, 40027, 40028, 40029, 40030, 40032, 40034, 40035, 50003, 50004, 
	130035, 140015, 200007, 210029, 210030, 210031, 210032, 210033, 210034, 210035, 
	220034, 220035, 230030, 230031, 230032, 230033, 230034, 230035, 240018, 240026,
	240027, 240028, 240029, 240030, 240032, 240034, 240035, 250000, 300008, 330033,
	410012, 410029, 410030, 410031, 410032, 410033, 410034, 410035, 420033, 420034, 

	420035, 430015, 430019, 430023, 430025, 430027, 430029, 430030, 430031, 430032,
	430033, 430034, 430035, 440015, 440024, 440025, 440026, 440027, 440028, 440029, 
	440030, 440032, 440033, 440034, 440035, 500022, 510034, 510035, 520003, 520035,
	530031, 530033, 530035, 540032, 540034, 610029, 610030, 610031, 610032, 610033,
	610034, 610035, 620033, 620034, 620035, 630013, 630015, 630017, 630019, 630021,
	630023, 630025, 630027, 630028, 630029, 630030, 630031, 630032, 630033, 630034,
	630035, 640015, 640017, 640021, 640023, 640024, 640025, 640026, 640027, 640028,
	640029, 640030, 640032, 640033, 640034, 640035, 700034, 710032, 710034, 730027,
	730031, 730033, 730035, 810008, 810029, 810030, 810031, 810032, 810033, 810034,
	810035, 820033, 820034, 820035, 830015, 830019, 830023, 830025, 830027, 830029,
	830030, 830031, 830032, 830033, 830034, 830035, 840024, 840025, 840026, 840027,
	840028, 840029, 840030, 840032, 840034, 840035, 850019, 900015, 930027, 930033, 930035, 
	1010027, 1010029, 1010030, 1010031, 1010032, 1010033, 1010034, 1010035, 1020002, 1020024, 
	1020034, 1020035, 1030015, 1030019, 1030025, 1030027, 1030028, 1030029, 1030030, 1030031, 
	1030032, 1030033, 1030034, 1030035, 1040015, 1040026, 1040027, 1040028, 1040029, 1040030,
	1040032, 1040033, 1040034, 1040035, 1110034, 1130029, 1130031, 1130033, 1130035, 1210029,
	1210030, 1210031, 1210032, 1210033, 1210034, 1210035, 1220024, 1220026, 1220034, 1220035,
	1230019, 1230021, 1230023, 1230025, 1230027, 1230029, 1230030, 1230031, 1230032, 1230033,
	1230034, 1230035, 1240013, 1240015, 1240017, 1240024, 1240025, 1240026, 1240027, 1240028,
	1240029, 1240030, 1240032, 1240034, 1240035, 1250001, 1250002, 1330017, 1330019, 1330033,

	1330035, 1340034, 1410025, 1410027, 1410029, 1410030, 1410031, 1410032, 1410033, 1410034,
	1410035, 1420033, 1420034, 1420035, 1430013, 1430019, 1430021, 1430023, 1430025, 1430027,
	1430029, 1430030, 1430031, 1430032, 1430033, 1430034, 1430035, 1440002, 1440009, 1440011,
	1440015, 1440017, 1440021, 1440024, 1440025, 1440026, 1440027, 1440028, 1440029, 1440030,
	1440032, 1440033, 1440034, 1440035, 1510034, 1530000, 1530031, 1530033, 1530035, 1540011,
	1610005, 1610029, 1610030, 1610031, 1610032, 1610033, 1610034, 1610035, 1620016, 1620033,
	1620034, 1620035, 1630019, 1630027, 1630029, 1630030, 1630031, 1630032, 1630033, 1630034,
	1630035, 1640015, 1640026, 1640027, 1640028, 1640029, 1640030, 1640032, 1640034, 1640035,
	1650010, 1700003, 1710009, 1710034, 1730031, 1730033, 1730035, 1750021, 1810011, 1810015,
	1810029, 1810030, 1810031, 1810032, 1810033, 1810034, 1810035, 1820033, 1820034, 1820035,
	1830019, 1830023, 1830025, 1830027, 1830028, 1830029, 1830030, 1830031, 1830032, 1830033,
	1830034, 1830035, 1840015, 1840023, 1840025, 1840026, 1840027, 1840028, 1840029, 1840030,
	1840032, 1840033, 1840034, 1840035, 1850010, 1910030, 1910032, 1910034, 1910035, 1930025,
	1930031, 1930033, 1930035, 1940034, 2010020, 2010027, 2010029, 2010030, 2010031, 2010032,
	2010033, 2010034, 2010035, 2020006, 2020033, 2020034, 2020035, 2030017, 2030019, 2030023,
	2030025, 2030027, 2030029, 2030030, 2030031, 2030032, 2030033, 2030034, 2030035, 2040015,
	2040023, 2040024, 2040025, 2040026, 2040027, 2040028, 2040029, 2040030, 2040032, 2040033,
	2040034, 2040035, 2050009, 2100002, 2110027, 2110032, 2130002, 2130025, 2130033, 2130035,
	2140034, 2210029, 2210030, 2210031, 2210032, 2210033, 2210034, 2210035, 2220034, 2220035,
	2230019, 2230027, 2230029, 2230030, 2230031, 2230032, 2230033, 2230034, 2230035, 2240009,
	2240025, 2240026, 2240027, 2240028, 2240029, 2240030, 2240032, 2240035, 2330027, 2330029,
	2330031, 2330033, 2330035, 2410029, 2410030, 2410031, 2410032, 2410033, 2410034, 2410035,
	2420015, 2420033, 2420034, 2420035, 2430021, 2430023, 2430025, 2430027, 2430029, 2430030,
	2430031, 2430032, 2430033, 2430034, 2430035, 2440015, 2440023, 2440025, 2440026, 2440027,
	2440028, 2440029, 2440030, 2440032, 2440034, 2440035, 2500014, 2510032, 2510034, 2510035,
	2530021, 2530031, 2530033, 2530034, 2530035, 2610004, 2610027, 2610029, 2610030, 2610031,
	2610032, 2610033, 2610034, 2610035, 2620033, 2620034, 2620035, 2630013, 2630017, 2630019,
	2630021, 2630023, 2630025, 2630027, 2630028, 2630029, 2630030, 2630031, 2630032, 2630033,
	2630034, 2630035, 2640009, 2640015, 2640023, 2640025, 2640026, 2640027, 2640028, 2640029,
	2640030, 2640032, 2640034, 2640035, 2700009, 2710034, 2730031, 2730033, 2730035, 2740034,
	2800021, 2800032, 2810025, 2810029, 2810030, 2810031, 2810032, 2810033, 2810034, 2810035,
	2820006, 2820026, 2820033, 2820034, 2820035, 2830027, 2830029, 2830030, 2830031, 2830032,
	2830033, 2830034, 2830035, 2840015, 2840023, 2840025, 2840026, 2840027, 2840028, 2840029,
	2840030, 2840032, 2840034, 2840035, 2920023, 2930019, 2930033, 2930035, 2940011, 2940021,
	2940034, 2950012, 3010023, 3010029, 3010030, 3010031, 3010032, 3010033, 3010034, 3010035,
	3020033, 3020034, 3020035, 3030019, 3030027, 3030028, 3030029, 3030030, 3030031, 3030032,
	3030033, 3030034, 3030035, 3040015, 3040025, 3040026, 3040027, 3040028, 3040029, 3040030,
	3040032, 3040034, 3040035, 3050015, 3130018, 3130035

};
*/

/*
void MCDigi::Process(string inFileName, string outFileName, string outHistFileName)
{  
    inFile = TFile::Open(TString(inFileName),"READ");
    if(!inFile) cout << "Open file failed!"<< endl;
	TIter next(inFile->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)next())) {
		if (strcmp(key->GetClassName(), "TTree") == 0) { 
			inTree = (TTree*)inFile->Get(key->GetName()); 
			cout << "Tree name: " << inTree->GetName() << endl; 
		}
	}
    if(!inTree) cout << "Open tree failed!"<< endl;
	inTree->SetBranchAddress("CellID",&_readCellID);
	if(strcmp(inTree->GetName(),"Simu_Digi_Hit")==0)
	{
    	inTree->SetBranchAddress("Hit_High_ADC",&_readHitHighADC);
    	inTree->SetBranchAddress("Hit_Low_ADC",&_readHitLowADC);
    	inTree->SetBranchAddress("Hit_High_Switch",&_readHitHighSwitch);
	}
    inTree->SetBranchAddress("Hit_Energy",&_readHitEnergy);
    inTree->SetBranchAddress("Hit_X",&_readHitX);
    inTree->SetBranchAddress("Hit_Y",&_readHitY);
    inTree->SetBranchAddress("Hit_Z",&_readHitZ);
	
	outFile 	= new TFile(TString(outFileName),"RECREATE");
	outHistFile = new TFile(TString(outHistFileName),"RECREATE");

	hTotal = new TH1D("AllEnergyDep"      , "AllEnergyDep", 1000, 0, 5000);	
	hEven  = new TH1D("EvenLayerEnergyDep", ""			  , 10000, 0, 10000);
	hOdd   = new TH1D("OddLayerEnergyDep" , ""			  , 10000, 0 ,10000);
	gStyle->SetOptStat(1);
	for (int i_cell = 0; i_cell < CELLNO; i_cell++)
	{
		int layerNo = i_cell / 210;
		int chipNo  = (i_cell%210) / 36;
		int chanNo  = (i_cell%210) % 36;
		ADCCell[i_cell] = new TH1D(TString(Form("ADC_Layer%d_Chip%d_Chan%d",layerNo,chipNo,chanNo)),	TString(Form("ADC_Layer%d_Chip%d_Chan%d",layerNo,chipNo,chanNo)),			400,	0, 4000);
		//MIPCell[i_cell] = new TH1D("MIPCell",	"MIPCell",	  500, 0,	200);
		energyCell[i_cell] = new TH1D(TString(Form("Energy_Layer%d_Chip%d_Chan%d",layerNo,chipNo,chanNo)),	TString(Form("Energy_Layer%d_Chip%d_Chan%d",layerNo,chipNo,chanNo)),	100,	0, 	2);
	}

	hitEnergy2D  = new TH2D("HitEnergy", "HitEnergy", LAYERNO, 0, LAYERNO, 150, 0, 300);
	hitNo2D		 = new TH2D("HitNo"    , "HitNo"    , LAYERNO, 0, LAYERNO, 120, 0, 120);
	hitEnergyPro = new TProfile("HitEnergyProfile", "HitEnergyProfile", LAYERNO, 0, LAYERNO);
	hitNoPro     = new TProfile("HitNoProfile"    , "HitNoProfile"    , LAYERNO, 0, LAYERNO);
    
	for (int i_entry = 0; i_entry < inTree->GetEntries(); i_entry++)
    //for (int i_entry = 0; i_entry < 2; i_entry++)
    {
		double totalSum = 0;
		double evenSum  = 0;
		double oddSum   = 0;
		for_each(layerHitEnergy	, layerHitEnergy + LAYERNO	, [](double& val) { val = 0; });
		for_each(layerHitNo		, layerHitNo + LAYERNO		, [](double& val) { val = 0; });

		inTree->GetEntry(i_entry);
		for (size_t i_vec = 0; i_vec < _readCellID->size(); i_vec++)
		{
			int cellIDNow = _readCellID->at(i_vec);
			int layerNo = cellIDNow / 100000;
			
			if (IsDeadChan(cellIDNow)) continue; // jiaxuan : Deadchannel cut
			if (1)
			//if (_readHitEnergy->at(i_vec) > 0.5*MIPEnergy) // jiaxuan : MIPEnergy cut
         	{
				layerHitNo[layerNo]++;
				layerHitEnergy[layerNo] += _readHitEnergy->at(i_vec);
				int histNo = CellID2HistNo(cellIDNow);
				int SiPMType = ((layerNo < 4) || (layerNo >27)) ? 1 : 0;
				// jiaxuan : unfinished
				if(strcmp(inTree->GetName(),"Simu_Digi_Hit")==0)
					ADCCell[histNo]->Fill(_readHitHighADC->at(i_vec));
				energyCell[histNo]->Fill(_readHitEnergy->at(i_vec));
				totalSum += _readHitEnergy->at(i_vec);
				//if ((_readCellID->at(i_vec) / 100000) % 2 == 0) 
				//	evenSum += _readHitEnergy->at(i_vec);
				//else
				//	oddSum  += _readHitEnergy->at(i_vec);
			}
		}
		for (int i_layer = 0; i_layer < LAYERNO; i_layer++)
		{
			hitEnergy2D  ->Fill(i_layer, layerHitEnergy[i_layer]);
			hitNo2D	     ->Fill(i_layer, layerHitNo[i_layer]);
			hitEnergyPro ->Fill(i_layer, layerHitEnergy[i_layer]);
			hitNoPro	 ->Fill(i_layer, layerHitNo[i_layer]);
		}
		hTotal->Fill(totalSum);
		//hEven->Fill(evenSum);
		//hOdd->Fill(oddSum);
    }

	double constant = hTotal->GetBinContent(hTotal->GetMaximumBin());
	double mean     = hTotal->GetMean();
	double sigma    = hTotal->GetRMS();
	//cout << constant << " " << mean << " " << sigma << endl;
	TF1 *fgaus = new TF1("fgaus","gaus");
	fgaus->SetParameters(constant, mean, sigma);
	fgaus->SetLineWidth(3);
	hTotal->Fit(fgaus);
	// TH1F *h2 = new TH1F();
	// hTotal->Copy(*h2);
	// TF1 *f2 = new TF1("f2","crystalball");
	// f2->SetLineColor(4);
	// f2->SetParameters(constant,mean,sigma,2,1);
	// f2->SetLineWidth(3);
	// h2->Fit(f2); 
	
	outHistFile->cd();
	
	for (int i_layer = 0; i_layer < LAYERNO; i_layer++)
	{
		for (int i_chip = 0; i_chip < CHIPNO; i_chip++)
		{
			if (strcmp(inTree->GetName(),"Simu_Digi_Hit")==0 ) {
				//cout << "Tree Name is Simu_Digi_Hit" << endl;
				outHistFile->mkdir(TString(Form("ADCCell/layer%d/chip%d",i_layer,i_chip)));
				outHistFile->cd(TString(Form("ADCCell/layer%d/chip%d",i_layer,i_chip)));
				for (int i_chan = 0; i_chan < CHANNO; i_chan++)
				{
					if(i_chip==5 && i_chan>29) continue;
					int histNo = i_layer*210+i_chip*36+i_chan;
					if(ADCCell[histNo]->GetEntries()>100)
					{
						TF1 *fun = new TF1("landauGaus",langaufun,300,2000,4);
						fun->SetParName(0, "LandauWidth");
						fun->SetParName(1, "LandauMPV");
						fun->SetParName(2, "Area");
						fun->SetParName(3, "GausSigma");
						if (i_layer>3 && i_layer<28) {
							fun->SetParameter(0, 20);
							fun->SetParameter(1, 500);
							fun->SetParameter(2, 10000);
							fun->SetParameter(3, 2);
							fun->SetParLimits(0, 10, 100);
							fun->SetParLimits(1, 350, 800);
							fun->SetParLimits(3, 0, 50);
							ADCCell[histNo]->GetXaxis()->SetRangeUser(400,1000);
							fun->Draw("same");
							gStyle->SetOptFit(1111);
						} else {
							fun->SetParameter(0, 65);
							fun->SetParameter(1, 900);
							fun->SetParameter(2, 25000);
							fun->SetParameter(3, 50);
							fun->SetParLimits(0, 20, 200);
							fun->SetParLimits(1, 500, 1300);
							fun->SetParLimits(3, 0, 200);
							ADCCell[histNo]->GetXaxis()->SetRangeUser(500,1800);
							fun->Draw("same");
							gStyle->SetOptFit(1111);
						}
						ADCCell[histNo]->Fit(fun,"QR","",300,2000);
						ADCCell[histNo]->SetLineWidth(2);
						ADCCell[histNo]->Draw();
						ADCCell[histNo]->Write();
					}
				}
			}
			
			outHistFile->mkdir(TString(Form("energyCell/layer%d/chip%d",i_layer,i_chip)));
			outHistFile->cd(TString(Form("energyCell/layer%d/chip%d",i_layer,i_chip)));
			for (int i_chan = 0; i_chan < CHANNO; i_chan++)
			{
				if(i_chip==5 && i_chan>29) continue;
				int histNo = i_layer*210 + i_chip*36 + i_chan;
				// cout << i_layer << " " << i_chip << " " << i_chan << endl;
				// cout << "histNo :" << histNo << endl;
				if(energyCell[histNo]->GetEntries()>100) 
					energyCell[histNo]->Write();
			}
		}
	}

	outFile->cd();
	gStyle->SetOptStat(0);
	
	// hOdd->SetLineColor(kBlue);
	// hOdd->SetLineWidth(2);
	// hOdd->GetXaxis()->SetTitle("MeV");
    // hOdd->GetYaxis()->SetTitle("Entries");
	// hEven->SetLineColor(kBlue);
	// hEven->SetLineWidth(2);
	// hEven->GetXaxis()->SetTitle("MeV");
	// hEven->GetYaxis()->SetTitle("Entries");
	
	hTotal->SetLineColor(kBlue);
	hTotal->SetLineWidth(2);
	hTotal->GetXaxis()->SetTitle("MeV");
	hTotal->GetYaxis()->SetTitle("Entries");
	hTotal->Write();
	//hEven->Write();
	//hOdd->Write();
	
	TCanvas * c1 = new TCanvas("HitEnergy", "HitEnergy", 1000, 800);
	TCanvas * c2 = new TCanvas("HitNo"    , "HitNo"    , 1000, 800);
	TCanvas * c3  = new TCanvas("EventEnergy","EventEnergy",800,600);
	//TLegend * lg3 = new TLegend(0.16,0.75,0.42,0.85);
	
	hitEnergy2D->SetOption("colz");
	hitEnergy2D->GetXaxis()->SetTitle("LayerNo");
	hitEnergy2D->GetYaxis()->SetTitle("Energy/MeV");
	hitEnergyPro->SetLineWidth(2);
	hitEnergyPro->SetLineColor(kRed);
	c1->cd();
	gPad->SetLogz();
	hitEnergy2D->Draw();
	hitEnergyPro->DrawCopy("sameA");
	hitEnergyPro->Write();
	
	hitNo2D->SetOption("colz");
	hitNo2D->GetXaxis()->SetTitle("LayerNo");
	hitNo2D->GetYaxis()->SetTitle("HitNo");
	hitNoPro->SetLineWidth(2);
	hitNoPro->SetLineColor(kRed);
	c2->cd();
	gPad->SetLogz();
	hitNo2D->Draw();
	hitNoPro->DrawCopy("sameA");
	hitNoPro->Write();
	
	//lg3->AddEntry(fgaus,"Gaussian Fit");
	//lg3->AddEntry(f2,"Crystal_Ball_Func Fit");
	c3->cd();
	c3->SetGridy();
	hTotal->SetTitle("EnergyDep");
	hTotal->GetXaxis()->SetRangeUser(mean-5*fabs(sigma), mean+5*fabs(sigma));
	hTotal->Draw();
	fgaus->Draw("sameL");
	//h2->Draw("same");
	//f2->Draw("sameL");
	//lg3->Draw("same");
	
	TPaveText *pt1 = new TPaveText(.16,.65,.42,.8,"blNDC");
	TText * tt1 = pt1->AddText("Gaussian Fit : ");
	TText * tt2 = pt1->AddText(Form("Mean : %.2f #pm %.2f", fgaus->GetParameter(1), fgaus->GetParError(1)));
	TText * tt3 = pt1->AddText(Form("Sigma : %.2f #pm %.2f", fgaus->GetParameter(2), fgaus->GetParError(2)));
	// TText * tt30 = pt1->AddText("");
	// TLine * line = pt1->AddLine(0,0,0,0);
	// line->SetLineColor(1);
	// line->SetLineWidth(4);
	// TText * tt4 = pt1->AddText("Crystal_Ball_Func Fit : ");
	// TText * tt5 = pt1->AddText(Form("Mean : %.2f #pm %.2f", f2->GetParameter(1), f2->GetParError(1)));
	// TText * tt6 = pt1->AddText(Form("Sigma : %.2f #pm %.2f", f2->GetParameter(2), f2->GetParError(2))); 
	pt1->Draw();
	c1->Write();
	c2->Write();
	c3->Write();
	outHistFile->Close();
	outFile->Close();
	delete c1;
	delete c2;
	delete c3;

	cout << "=============== !!! Process FINISHED !!! ==============" << endl;

}
*/

