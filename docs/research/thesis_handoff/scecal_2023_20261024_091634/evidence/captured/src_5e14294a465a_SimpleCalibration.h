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

class treeCalibration
{
    // for draw some histrogram
    private:
    void Init();
    void Reset();
    void PedestalExtract(string pedestal_name);
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
    
    // for EventAlign
    TFile* rawfile;
    TTree* fNtuple2;
    int _run_num, _event_time, _triggerID;
    std::vector<int>*    _cellID=nullptr;
    std::vector<int>*    _hitTag=nullptr;
    std::vector<int>*    _gainTag=nullptr;
    std::vector<double>* _hg_charge=nullptr;
    std::vector<double>* _lg_charge=nullptr;
    std::vector<double>* _hit_time=nullptr;
    std::vector<vector<double> >* _tempLayer=nullptr;

    // for Calibration
    TFile* califile;
    TTree* tree;
	int _run_num_new, _event_time_new, _event_num, _detectorID;
    std::vector<int> _cellID_new;
    std::vector<double> _hg_charge_new;
    std::vector<double> _lg_charge_new;
    std::vector<double> _temp;
    std::vector<double> _posX_cali;
    std::vector<double> _posY_cali;
    std::vector<double> _posZ_cali;
    TH2D *h_GainRatio[layerNu][chipNu][channelNu];
    //MIP Fitting
    
public:
    int Calibration(string root_name,string pedestal_name,string cali_name);
    treeCalibration(){};
    int MIPFit(string input_name,string hist_name,string mip_name, int temperatureCalib, int angularCalib);

};

#endif
