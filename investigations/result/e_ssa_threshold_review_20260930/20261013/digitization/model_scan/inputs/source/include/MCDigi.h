#ifndef _MCDigi_H
#define _MCDigi_H

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "Extract.h"

using namespace std; 

static constexpr const double MIPEnergy = 0.305; // MeV
static const int LAYERNO = 30;
static const int CHIPNO  = 6 ;
static const int CHANNO  = 36;
static const int CELLNO = 6300;
// static constexpr const Double_t MIPResponse[2] = {7., 11.};   // p.e. per  MIP
// static constexpr const Double_t PEChargeMean[2] = {14., 22.7};// ADC  per p.e.
// static constexpr const Double_t PEChargeSigma[2] = {5.6, 10.};
// static constexpr const Double_t PEChargeNoise = 3;

struct Point{
	double posX;
	double posY;
	double posZ;
	Point(double a=0., double b=0., double c=0.) {
		posX = a;
		posY = b;
		posZ = c;
	}
};

struct Result{
	double hADC;
	double lADC;
	Result(double a=0., double b=0.) {
		hADC = a;
		lADC = b;
	}
};

static constexpr const int decodeID[CHIPNO][CHANNO] =
{  // decodeID[chipNo][chanNo] = copyNo (in a layer)
	0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13,
	55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,

	24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 
	33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,

	149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 
	108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,

	150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 
	159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
	
	191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181,
	139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
	
	137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 
	93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, -1, -1, -1, -1, -1, -1
}; // chip0_chan1 -> 42 

class MCDigi {
public:
	TFile* inFile;
	TTree* inTree;
	TFile* outFile;
	TTree* outTree;

	//TH2D* hitEnergy2D;
	//TH2D* hitNo2D;
	//TProfile* hitEnergyPro;
	//TProfile* hitNoPro;
	//TH1D* ADCCell[CELLNO];
	//TH1D* energyCell[CELLNO];
	//TH2D * photon2pe[LAYERNO][CHIPNO][CHANNO]; 

	double layerHitEnergy[LAYERNO];
	double layerHitNo[LAYERNO];
	
	TH1D * test11;
	TH1D * test12;
	TH1D * test13;
	TH1D * test14;
	TH2D * test21;
	TH2D * test22;

	int length;
	int area;
	//vector<vector<int> > isp;
	int n_fired;
	double pCrosstalk;
	vector<vector<int> > isa;

public:	
	void InDigiBranchClear();
	Point CellID2Position(int cellID);
	int CopyNo2CellID(int copyNo); // copyNo: 0-6719
	//bool IsDeadChan(int cellID);
	void PixelAva(int npx, int npy); 
	int PixelCorrection(int nPhoton, int nPixel, double pde, double pCT);
	Result StepDigi(double energyDep, int cellID, Extract& ex);
	bool ApplyThreshold(double adc_raw, double thre, double threSigma, double pedestal);
	void SimuDigi(string inFileName, string outFileName, Extract& ex); 
	// void SimuDigi(string inFileName, string outFileName, string pedFileName, string HLRatioFileName, string MIPFileName, string SiPMGainFileName);
	// void Process(string inFileName, string outFileName, string outHistFileName);

private:
	//	0	1	2
	//	7	T	3
	//	6	5	4
	

private:
	Int_t _digiRunNo, _digiEventTime, _digiCycleID, _digiTriggerID;
	vector<int>	   _digiCellID;
	vector<int>	   _digiBCID;
	vector<int>	   _digiHitTag;
	vector<int>	   _digiGainTag;
	vector<double> _digiHitHighADC;
	vector<double> _digiHitLowADC;
	vector<double> _digiHitTime;
	vector<vector<double>> _digiTemperature;
	
	vector<int>*	_readCellID = nullptr;
	vector<double>* _readHitHighADC = nullptr;
	vector<double>* _readHitLowADC = nullptr;
	vector<double>* _readHitHighSwitch= nullptr;
	vector<double>* _readHitEnergy = nullptr;
	vector<double>* _readHitX = nullptr;
	vector<double>* _readHitY = nullptr;
	vector<double>* _readHitZ = nullptr;

};

#endif
