#ifndef TreeSelectRoot_h_
#define TreeSelectRoot_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH3.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TMinuit.h>
#include<vector>
#include<TF2.h>

using namespace std;

class trackFit
{

  public:
    static const int layerNu = 32;
    static const int chipNu = 6;
    static const int chnNu = 36;
    const double _xInterval = 5.3; // 300 um gap in width direction
    const double _yInterval = 45.4; // 400 um gap in length direction
    const int rowNu = 42;
    const int columnNu = 5;
    double EBUposX[16]={0.195964, -0.0818129, 0.00573527, 0.0400398, -0.13863, 0.0408884, -0.121727, -0.0522908, 0.0647114, -0.180934, 0.0555704, 0.00811854, -0.00357011, 0.0576592, -0.0627433, 0.138101};
    double EBUposY[16]={-0.00311279, 0.0613501, -0.0762173, 0.0746585, -0.134201, 0.0625489, -0.0402717, 0.189599, 0.067039, -0.154837, 0.0120084, -0.173652, -0.0551821, 0.124792, -0.046264, 0.0852524};
    double EBUposX2[16]={-0.0797351, 0.0219277, -0.0258218, -0.028357, 0.0272062, 0.0042528, 0.0177905, 0.033127, -0.00132277, 0.0465462, 0.0010297, -0.00559773, 0.0156754, -0.0344404, 0.0497271, -0.087874};
    double EBUposY2[16]={-0.00318017, -0.0295886, 0.0514142, -0.0447454, 0.052376, 0.00375967, -0.0176469, -0.0213683, -0.0178936, 0.0117386, -0.0314558, 0.0100134, 0.0444055, -0.0246971, 0.063703, -0.045664};

  private:
    vector<int> _hitCell;
    vector<double> _hitPosX;
    vector<double> _hitPosY;
    vector<double> _hitPosZ;
    vector<double> _hitCharge;
    double *_position;
    vector<int> _hitCellnew;
    vector<double> _hitPosXnew;
    vector<double> _hitPosYnew;
    vector<double> _hitPosZnew;
    vector<int> _hitCellnew_out;
    vector<double> _hitPosXnew_out;
    vector<double> _hitPosYnew_out;
    vector<double> _hitPosZnew_out;
    

    int _EventNum, _DetectorID;
    std::vector<int>* _CellID = nullptr;
    std::vector<double>* _Hit_Energy = nullptr;
    std::vector<double>* _HG_Charge = nullptr;
    std::vector<double>* _LG_Charge = nullptr;
    std::vector<double>* _Hit_Time = nullptr;
    std::vector<double>* _Hit_X = nullptr;
    std::vector<double>* _Hit_Y = nullptr;
    std::vector<double>* _Hit_Z = nullptr;
    std::vector<double>* _NewTemperature = nullptr;

    void Reset();
    void Clear();

  private:
    TFile *rootfile;
    TTree *fNtuple1;

    // track fitting
    TFile *fitfile;
    TTree *fTree;
    int _EventNum_new, _DetectorID_new;
    vector<int> _CellID_new;
    vector<double> _Hit_Energy_new;
    vector<double> _HG_Charge_new;
    vector<double> _LG_Charge_new;
    vector<double> _Hit_Time_new;
    vector<double> _Hit_X_new;
    vector<double> _Hit_Y_new;
    vector<double> _Hit_Z_new;
    vector<double> _temp;
    vector<double> _temp_new;
    vector<double> _trackFitPars;
    vector<double> _residualX;
    vector<double> _residualY;
    vector<double> _chargeCor;
    vector<double> _timeCor;
    vector<int> _layer_efficiency;
    vector<int> _trackTag;
    double _realLength;
    double _ndf;

    TMinuit* gMinuitX;
    TMinuit* gMinuitY;

  public:
    trackFit(string root_name,string map_name,string event);

    double* TrackFit(int slience);
};

#endif
