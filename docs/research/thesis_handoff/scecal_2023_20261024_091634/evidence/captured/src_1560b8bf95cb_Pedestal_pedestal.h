#ifndef TreePedestal_pedestal_h_
#define TreePedestal_pedestal_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>


using namespace std;

class treePedestalPed
{

    // for draw some histrogram
private:
    void pedInit();
    void Reset();
    TH2F* hh1;
    TH2F* hh2;
    TH1F *hPed[6];
    vector<int>* _chan1;
    vector<int>* _chan2;
private:
    static const int layerNu =32;
    static const int chipNu = 6;
    static const int channelNu = 36;
    static const int memoryNu = 16;
    // for fitting parameters
    std::vector<int> _pedCellID;
    std::vector<float> _pedMean;
    std::vector<float> _pedSigma;
    std::vector<float> _pedMeanCharges;
    std::vector<float> _pedSigmaCharges;
    std::vector<float> _pedMeanTimes;
    std::vector<float> _pedSigmaTimes;
//    std::vector<float> _pedTempSlope;  // for temperature dependence
//    std::vector<float> _pedTempIntercept;  // for temperature dependence
    std::vector<float> _pedChiCharges;
    std::vector<float> _pedChiTimes;

    int _leftC;
    int _righC;
    int _leftT;
    int _righT;

public:
    // for rawRoot
    TFile* rawfile;
    TTree* tree;
    double _cycleID;
    std::vector<int>* _cellID;
    std::vector<int>* _bcid;
    std::vector<int>* _hitTag;
    std::vector<int>* _gainTag;
    std::vector<double>* _charge;
    std::vector<double>* _time;
    
    ~treePedestalPed();
    int Pedestal(string root_name,string pedestal_name,string fitdat_name);
//    void TemperatureDependence(string root_name,string pedestal_name,string fitdat_name);  // for temperature dependence

    private:
    TF1* chargefun;
    // TF1* chargefun_double;
    TF1* timefun;
    TF1* timefun_d;
    TFile *fitfile;
    TFile *selectfile;
    // Memory Cells Level
    TTree *fNtuple1;
    TH1F *hCharge1[layerNu][chipNu][channelNu][memoryNu];
    TH1F *hTime1[layerNu][chipNu][channelNu][memoryNu];
    TObjArray *Hlist1;
    // Channels Level
    TTree *fNtuple2;
    TH1F *hCharge2[layerNu][chipNu][channelNu];
    TH1F *hTime2[layerNu][chipNu][channelNu];
    TObjArray *Hlist2;
    // Chips Level
    TTree *fNtuple3;
    TH1F *hCharge3[layerNu][chipNu];
    TH1F *hTime3[layerNu][chipNu];
    TObjArray *Hlist3;

};

#endif
