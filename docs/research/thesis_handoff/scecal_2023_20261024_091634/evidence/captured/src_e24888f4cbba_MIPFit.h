#ifndef TreeCalibration_h_
#define TreeCalibration_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>

using namespace std;

class MIPCalibration
{
    // for draw some histrogram
    private:
		// from auto gain mode decode 
    	Int_t _run_no, _event_time, _cycleID, _triggerID;
    	vector< int >* _cellID;
    	vector< int >* _bcid;
    	vector< int >* _hitTag;
    	vector< int >* _gainTag;
    	vector< double >* _chargehi;
    	vector< double >* _chargelo;
    	vector< double >* _hit_time;
    	std::vector<vector<double> >* _temp;

		double _landauMPV, _landauWidth, _gausSigma, _ChiSqr;
		int _Ndf, _cellID_mip, _chn_entries;

    private:
		static const int layerNu =32;
		static const int chipNu = 6;
		static const int channelNu = 36;
		static const int memoryNu = 16;
		
		// for fitting parameters
		TFile* fitfile;
		TTree* fNtuple1;
		std::vector<int>* _pedCellID;
		std::vector<float>* _pedMeanCharges;
		std::vector<float>* _pedSigmaCharges;
		std::vector<float>* _pedMeanTimes;
		std::vector<float>* _pedSigmaTimes;
		double pedMeanChargesChnLevel[layerNu][chipNu][channelNu];
		double pedSigmaChargesChnLevel[layerNu][chipNu][channelNu];
		double pedMeanTimesChnLevel[layerNu][chipNu][channelNu];
		double pedSigmaTimesChnLevel[layerNu][chipNu][channelNu];

	public:
    	MIPCalibration(){};
    	int MIPFit(string input_name,string hist_name,string mip_name);
		void ReadTree(TTree* tree);
		void SetTree(TTree* tree);
		//void Init();
		//void Reset();
		void PedestalExtract(string pedestal_name);

};

#endif
