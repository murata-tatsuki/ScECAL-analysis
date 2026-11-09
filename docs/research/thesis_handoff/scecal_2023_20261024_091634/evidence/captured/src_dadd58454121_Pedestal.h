#ifndef TreePedestal_h_
#define TreePedestal_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include"TH2.h"


using namespace std;

class treePedestal
{
	private:
        void pedInit();
        void Reset();
    	TH2D* low_mean;
    	TH2D* high_mean;
    	TH2D* low_rms;
    	TH2D* high_rms;
    	TH2D* high_mean_chip;
    	TH2D* low_mean_chip;
    	TH2D* high_rms_chip;
    	TH2D* low_rms_chip;
    private:
        static const int layerNu =32;
        static const int chipNu = 6;
        static const int channelNu = 36;
        static const int memoryNu = 16;
        std::vector<int> _pedCellID;
        std::vector<double> _pedMeanCharges;
        std::vector<double> _pedSigmaCharges;
        std::vector<double> _pedMeanTimes;
        std::vector<double> _pedSigmaTimes;
        int _leftC,_righC;
    	//note: che
        //std::vector<double> _pedMean;
        //std::vector<double> _pedSigma;
		//int _leftT,_righT;
    public:
        TFile* rawfile;
        TTree* tree;
        int _run_num,_event_time,_cycleID,_triggerID;
        std::vector<int>* _cellID = nullptr;
        std::vector<int>* _hitTag = nullptr;
        std::vector<int>* _gainTag = nullptr;
        std::vector<double>* _lg_charge = nullptr;
        std::vector<double>* _hg_charge = nullptr;
        
        ~treePedestal();
        int Pedestal(string root_name,string pedestal_name,string fitdat_name);
        TF1* gausfun;
        TFile *fitfile;
        TFile *pedestalfile;
        // Memory Cells Level
        //TTree *fNtuple1;
        //TH1D *hCharge1[layerNu][chipNu][channelNu][memoryNu];
        //TH1D *hTime1[layerNu][chipNu][channelNu][memoryNu];
        //TObjArray *Hlist1;
        // Channels Level
        TTree *fNtuple2;
        TH1D *hCharge2[layerNu][chipNu][channelNu];
        TH1D *hTime2[layerNu][chipNu][channelNu];
        TObjArray *Hlist2h[layerNu][chipNu];
        TObjArray *Hlist2l[layerNu][chipNu];
        // Chips Level
        //TTree *fNtuple3;
        //TH1D *hCharge3[layerNu][chipNu];
        //TH1D *hTime3[layerNu][chipNu];
        //TObjArray *Hlist3;

};

#endif
