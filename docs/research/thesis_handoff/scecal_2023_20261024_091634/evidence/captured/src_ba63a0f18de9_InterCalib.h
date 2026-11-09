#ifndef TreeSelectRoot_h_
#define TreeSelectRoot_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include<TGraph.h>

using namespace std;

class interCalib
{

    public:
    void Init();
    interCalib(string root_name,string sps_name);

    void caliSPS(int chnID,vector<double> meansps, vector<double> sigmasps);

    protected:
    static const int chipNu = 6;
    static const int channelNu = 36;
    static const int memoryNu = 2;
    static const int LayerNo=32;
    const string Elink[LayerNo] = {"-040","-041","-042","-043","-044","-045","-046","-047","-048","-049","-04A","-04B","-04C","-04D","-04E","-04F",
                                   "-050","-051","-052","-053","-054","-055","-056","-057","-058","-059","-05A","-05B","-05C","-05D","-05E","-05F"};
    // Channels Level
    TObjArray *Hlist;
    TObjArray *Glist = new TObjArray(210);
    TGraph *gp[LayerNo][chipNu][channelNu];
	TH2F* hisg[LayerNo][chipNu][channelNu];
    int eventNo[LayerNo][chipNu][channelNu];
	int Xmax[LayerNo][chipNu][channelNu];

    int _channelID, _Npeaks, _DAC;
    double _Interval;
    vector<double> _spsPar;
    double _slope, _slopeError;
    double _intercept, _interceptError;
    double _chiSquare,_NDF;
    TF1 *fit;

    vector<int> _selectID;
    int _BlindChannel;
   
   	std::vector<int>* _cellID=nullptr;
   	std::vector<int>* _hitTag=nullptr;
    std::vector<double>* _lg_charge=nullptr;
    std::vector<double>* _hg_charge=nullptr;

    private:
    TFile* rawfile;

};

#endif
