#ifndef Treew_h_
#define Treew_h_

//#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<fstream>
#include<iostream>
#include <TCanvas.h>
#include <TH2.h>
//#include "Event.h"

using namespace std;

class raw2Root
{
    private:
    static const int channel_FEE = 73;//(36charges+36times + BCIDs )*16column+ ChipID
    static const int cell_SP = 16;
    static const int chip_No = 6; // Yuzhi: Changed to ScECAL
    static const int channel_No = 36;
    static const int Layer_No = 32; // Yuzhi: Changed to ScECAL
    static const int temp_No = 16;
    const double _Pos_X[channel_No]={100.2411,100.2411,100.2411,59.94146,59.94146,59.94146,19.64182,19.64182,19.64182,19.64182,59.94146,100.2411,100.2411,59.94146,19.64182,100.2411,59.94146,19.64182,-20.65782,-60.95746,-101.2571,-20.65782,-60.95746,-101.2571,-101.2571,-60.95746,-20.65782,-20.65782,-20.65782,-20.65782,-60.95746,-60.95746,-60.95746,-101.2571,-101.2571,-101.2571};
    const double _Pos_Y[channel_No]={141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,261.94766,261.94766,261.94766,302.2473,302.2473,302.2473,342.54694,342.54694,342.54694,342.54694,342.54694,342.54694,302.2473,302.2473,302.2473,261.94766,261.94766,261.94766,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874};
    const double chip_dis_X=239.3;
    const double chip_dis_Y=241.8;
    const double HBU_X=239.3;
    const double HBU_Y=725.4;    
    int LayerNo;
    vector<int> chipNo;

    public:
    int Decode(string binary_name,string raw_name);
    //int BER(string binary_name);
    //int RMFelixTag(string str_datalist,string outputDir);
    int CatchSPIROCBag(ifstream &f_in,vector<int> &buffer_v, int &layer_id,int &cycleID,int &triggerID);
    int FillChipBuffer(vector<int> &buffer_v,int cycleID,int triggerID,int layer_id);
    void SetTreeBranch(TTree *tree);
    //void ReadTreeBranch(TTree *tree);
    // int FillTreeBranch(TTree *tree,vector<int> &buffer_v,int cycleID,int triggerID,int layer_id);
    int DecodeAEvent(vector<int> &chip_v,vector<int> &tmpr_v,int layer_ID,int Memo_ID);
    //int AnaSPS(string str_in,string outputDir);    
    //int AnaPed(string str_in,string outputDir);    
    //int AnaCos(string str_in,string outputDir);    
    //void Reset();
    void BranchClear();
    //void deleteLocator();
    //void triggerBuild(vector<int> bufferAsic);
    double cycleIDbuffer;
    int prevCycleID;
    int tluOverflow;
    double Pos_X(int channel_ID,int chip_ID,int HBU_ID=1);
    double Pos_Y(int channel_ID,int chip_ID,int HBU_ID=1);
    int Chipbuffer_empty();
    private:
//    Event* EventClass;
   ifstream _f_in;
   int  _layer_id=0;
   vector< int > _buffer_v;
   vector< int > _chip_v[Layer_No][chip_No];
   vector< int > _tmpr_v[Layer_No];
   //vector< int > _hitTagOri; ///< Debug: original hit tag
   vector< double > _time;
   int _temp_checkfirst[Layer_No];
   
   Int_t _run_no;
   Int_t _event_time;
   Int_t _cycleID;
   Int_t _triggerID;
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
