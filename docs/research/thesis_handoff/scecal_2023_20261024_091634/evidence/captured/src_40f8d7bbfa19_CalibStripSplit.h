#ifndef TreeCalibStripSplit_h_
#define TreeCalibStripSplit_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include <map>

using namespace std;

class treeCalibStripSplit
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
    const double _xInterval = 5.3;
    const double _yInterval = 45.4;
    const int _splitNu = 9;
    const double _sInterval = 5.;
    
    // for Calibration
    TFile* califile;
    TTree* fNtuple2;
    // double _cycleID, _triggerID;
    // std::vector<int>* _cellID=nullptr;
    // std::vector<int>* _bcid=nullptr;
    // std::vector<int>* _hitTag=nullptr;
    // std::vector<int>* _gainTag=nullptr;
    // std::vector<double>* _charge=nullptr;
    // std::vector<double>* _time=nullptr;
    // std::vector<double>* _temp=nullptr;
    // std::vector<double>* _posX=nullptr;
    // std::vector<double>* _posY=nullptr;
    // std::vector<double>* _posZ=nullptr;
    // std::vector<double>* _low=nullptr;
    // std::vector<double>* _high=nullptr;
    int _run_num, _event_num, _event_time, _detectorID;
    std::vector<int> *_cellID=nullptr;
    //std::vector<double> _hitTag_new;
    std::vector<double> *_hitE=nullptr;
    std::vector<double> *_posX=nullptr;
    std::vector<double> *_posY=nullptr;
    std::vector<double> *_posZ=nullptr;
    //std::vector<double> _hit_time_new;
    //std::vector<double> _hg_cali;
    //std::vector<double> _lg_cali;
    std::vector<double> *_temp=nullptr;

    // for SSA
    TFile* ssafile;
    TTree* tree;
    int _run_num_new, _event_num_new, _event_time_new, _detectorID_new;
    // double _cycleID_cali, _triggerID_cali;
    std::vector<int> _cellID_cali;
    // std::vector<int> _bcid_cali;
    // std::vector<int> _hitTag_cali;
    // // std::vector<int> _gainTag_cali;
    // std::vector<double> _charge_cali;
    // std::vector<double> _time_cali;
    // std::vector<double> _temp_cali;
    // std::vector<double> _posX_cali;
    // std::vector<double> _posY_cali;
    // std::vector<double> _posZ_cali;
    // std::vector<double> _low_cali;
    // std::vector<double> _high_cali;
    std::vector<double> _hitE_cali;
    std::vector<double> _posX_cali;
    std::vector<double> _posY_cali;
    std::vector<double> _posZ_cali;
    std::vector<double> _temp_cali;
    std::vector<int> _ssaTag;
    //MIP Fitting
    
public:
//    ~treeSSA();
    int StripSplit(string txt_name,string cali_name,string ssa_name);
    treeCalibStripSplit(){};

};

#endif
