#ifndef _AUTOGAINDECODE_H
#define _AUTOGAINDECODE_H

#include<TFile.h>
#include<TTree.h>
#include<fstream>
#include<iostream>
#include <TCanvas.h>
#include <TH2.h>

using namespace std;

class raw2RootAuto
{
	private:
    	static const int channel_FEE = 73;//(36charges+36times + BCIDs )*16column+ ChipID
    	static const int layer_No = 32;
    	static const int chip_No = 6;
    	static const int channel_No = 36;
    	static const int temp_No = 16;

    public:
		bool chip_buffer_flag;
		bool IsAva;
		char char_tmp[1000];
		int trigger_period;
		int run_no, event_time, Event_No;
    	int IsDataAvailable(int layerID, int chipID);
    	void BranchClear();
    	void SetTreeBranch(TTree *tree);
    	int CatchSPIROCBag(ifstream &f_in,vector<int> &buffer_v, int &cycleID,int &triggerID);
    	//int FillChipBuffer(vector<int> &buffer_v,int cycleID,int triggerID,int layer_id);
    	int DecodeAEvent(vector<int> &chip_v,vector<int> &tmpr_v,int layer_ID);
    	int AutoGainDecode(string binary_name,string raw_name);
   
    private:
    	vector< int > _buffer_v;
		vector< int > _chip_v[layer_No][chip_No];
    	vector< int > _tmpr_v[layer_No];
    	int _temp_checkfirst[layer_No];
   
    	Int_t _run_no, _event_time, _cycleID, _triggerID;
    	vector< int > _cellID;
    	vector< int > _bcid;
    	vector< int > _hitTag;
    	vector< int > _gainTag;
    	vector< double > _chargehi;
    	vector< double > _chargelo;
    	vector< double > _hit_time;
    	std::vector<vector<double> > _temp;
};

#endif
