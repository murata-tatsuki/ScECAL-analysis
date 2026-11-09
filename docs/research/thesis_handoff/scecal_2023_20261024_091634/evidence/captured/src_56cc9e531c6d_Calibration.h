#ifndef Calibration_h_
#define Calibration_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>

using namespace std;

class processCali
{
    // for draw some histrogram
    private:
    void Init();
    void Reset();
    private:
    static const int layerNu =32;
    static const int chipNu = 6;
    static const int channelNu = 36;
    static const int memoryNu = 16;
    static const int EBUchn = 210;
    // for fitting parameters
    TFile* fitfile;
    TTree* fNtuple1;
    std::vector<int>* _pedCellID;
    std::vector<float>* _pedMeanCharges;
    std::vector<float>* _pedSigmaCharges;
    std::vector<float>* _pedMeanTimes;
    std::vector<float>* _pedSigmaTimes;
    double pedMeanChargesChnLevel[layerNu][chipNu][channelNu]={0};
    double pedSigmaChargesChnLevel[layerNu][chipNu][channelNu]={0};
    double pedMeanTimesChnLevel[layerNu][chipNu][channelNu]={0};
    double pedSigmaTimesChnLevel[layerNu][chipNu][channelNu]={0};
    double interFactor[layerNu][chipNu][channelNu]={1};
    double interFactorIntercept[layerNu][chipNu][channelNu]={0};
    double mipFactor[layerNu][chipNu][channelNu]={1};
    bool mipFactorValid[layerNu][chipNu][channelNu]={false};
    bool interFactorValid[layerNu][chipNu][channelNu]={false};
    
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
    //std::vector<double> _hitTag_new;
    std::vector<double> _hitE;
    std::vector<double> _hit_hg_E;
    std::vector<double> _hit_lg_E;
    std::vector<double> _posx;
    std::vector<double> _posy;
    std::vector<double> _posz;
    //std::vector<double> _hit_time_new;
    //std::vector<double> _hg_cali;
    //std::vector<double> _lg_cali;
    std::vector<double> _temp;
    double _tot_energy;
    int _tot_hitlayer;
    //MIP Fitting
    
public:
    void PedestalExtract(string pedestal_name);
    void InterFactorExtract(string interfactor_name);
    void MIPFactorExtract(string mipfactor_name);
    void Calibration(string root_name,string pedestal_name,string interfactor_name, string mipfactor_name, string cali_name);
    void PedestalExtract_edit(string pedestal_name);
    void InterFactorExtract_edit(string interfactor_name);
    void MIPFactorExtract_edit(string mipfactor_name);
    void PedestalExtract_simulation(string pedestal_name);
    void InterFactorExtract_simulation(string interfactor_name);
    void MIPFactorExtract_simulation(string mipfactor_name);
    processCali(){};
    //int MIPFit(string input_name,string hist_name,string mip_name);

};

#endif
