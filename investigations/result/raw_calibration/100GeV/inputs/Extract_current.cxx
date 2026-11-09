#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TPaveText.h>
#include <TRandom3.h>
#include <TKey.h>
#include "Extract.h"

using namespace std; 

void Extract::PedestalExtract(string inPedFileName)
{
    cout<<"=========== PedestalExtract ============"<<endl;
	fill_n(&pedMeanHigh[0][0][0], 	LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&pedSigmaHigh[0][0][0], 	LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&pedMeanLow[0][0][0], 	LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&pedSigmaLow[0][0][0], 	LAYERNO*CHIPNO*CHANNO, 1);
    TFile* pedFile = new TFile(TString(inPedFileName),"READ");
    TTree* pedTree = (TTree*)pedFile->Get("ChnLevel"); 
    if(!pedTree)	{cout<<"!!! GET PED FILE FAILED !!!"<<endl;}
	vector<int>* _pedCellID = nullptr;
	vector<double>* _pedMeanHigh = nullptr;
	vector<double>* _pedSigmaHigh = nullptr;
	vector<double>* _pedMeanLow = nullptr;
	vector<double>* _pedSigmaLow = nullptr;
    pedTree ->SetBranchAddress("CellID",&_pedCellID);
    pedTree ->SetBranchAddress("PedHighMean",&_pedMeanHigh);
    pedTree ->SetBranchAddress("PedHighSig",&_pedSigmaHigh);
    pedTree ->SetBranchAddress("PedLowMean",&_pedMeanLow);
    pedTree ->SetBranchAddress("PedLowSig",&_pedSigmaLow);

	double meanPedHigh15  = 0 ; double meanPedHigh10  = 0;
	double sigmaPedHigh15 = 0 ; double sigmaPedHigh10 = 0;
    double chnNoPedHigh15 = 0 ; double chnNoPedHigh10 = 0;
	double meanPedLow15   = 0 ; double meanPedLow10   = 0;
	double sigmaPedLow15  = 0 ; double sigmaPedLow10  = 0;
    double chnNoPedLow15  = 0 ; double chnNoPedLow10  = 0;
	for(int entry=0; entry!=pedTree->GetEntries(); ++entry)
    {
        pedTree->GetEntry(entry);
        //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(size_t ihit=0; ihit!=_pedCellID->size(); ++ihit)
        {
            int layerID = _pedCellID->at(ihit)/1e5;
			if(layerID>=LAYERNO)	continue;// 
            int chipID = _pedCellID->at(ihit)%100000/1e4;
            int channelID = _pedCellID->at(ihit)%100;        
			
			if (layerID>3 && layerID<28) {
			   // if(_ChiSqr/(double)_NDF>2) continue;
				chnNoPedHigh10++;
				chnNoPedLow10++;
				meanPedHigh10 += _pedMeanHigh->at(ihit);
				meanPedLow10  += _pedMeanLow ->at(ihit);
				sigmaPedHigh10 += _pedSigmaHigh->at(ihit);
				sigmaPedLow10  += _pedSigmaLow ->at(ihit);
			} else {
			   // if(_ChiSqr/(double)_NDF>1.7) continue;
				chnNoPedHigh15++;
				chnNoPedLow15++;
				meanPedHigh15 += _pedMeanHigh->at(ihit);
				meanPedLow15  += _pedMeanLow ->at(ihit);
				sigmaPedHigh15 += _pedSigmaHigh->at(ihit);
				sigmaPedLow15  += _pedSigmaLow ->at(ihit);
			}
            pedMeanHigh[layerID][chipID][channelID] = _pedMeanHigh->at(ihit);
            pedSigmaHigh[layerID][chipID][channelID] = _pedSigmaHigh->at(ihit);
            pedMeanLow[layerID][chipID][channelID] = _pedMeanLow->at(ihit);
            pedSigmaLow[layerID][chipID][channelID] = _pedSigmaLow->at(ihit);
        }
    }
	meanPedHigh15 /= chnNoPedHigh15;
	meanPedLow15  /= chnNoPedLow15;
	meanPedHigh10 /= chnNoPedHigh10;
	meanPedLow10  /= chnNoPedLow10;
	sigmaPedHigh15 /= chnNoPedHigh15;
	sigmaPedLow15  /= chnNoPedLow15;
	sigmaPedHigh10 /= chnNoPedHigh10;
	sigmaPedLow10  /= chnNoPedLow10;
	cout << " HIGH GAIN PEDESTAL ... ..."  << endl;
    cout << " channel count(15um:10um) : " << chnNoPedHigh15 << " , " << chnNoPedHigh10 << endl;
    cout << " mean of pedestal High : "    << meanPedHigh15  << " , " << meanPedHigh10  << endl;
    cout << " sigma of pedestal High: "    << sigmaPedHigh15 << " , " << sigmaPedHigh10 << endl;
	cout << " LOW GAIN PEDESTAL ... ..."   << endl;
    cout << " channel count(15um:10um) : " << chnNoPedLow15  << " , " << chnNoPedLow10  << endl;
    cout << " mean of pedestal Low : "     << meanPedLow15   << " , " << meanPedLow10   << endl;
    cout << " sigma of pedestal Low: "     << sigmaPedLow15  << " , " << sigmaPedLow10  << endl;
	
	for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
                if (il>3 && il<28) {
					pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh10 : pedMeanHigh[il][ip][ic];
					pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1) 	 ? meanPedLow10 : pedMeanLow[il][ip][ic];
					pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh10 : pedSigmaHigh[il][ip][ic];
					pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow10 : pedSigmaLow[il][ip][ic];
                } else {
					pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh15 : pedMeanHigh[il][ip][ic];
					pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1) 	 ? meanPedLow15 : pedMeanLow[il][ip][ic];
					pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh15 : pedSigmaHigh[il][ip][ic];
					pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow15 : pedSigmaLow[il][ip][ic];
                }
            }
		}
	}
	pedFile->Close();
}

void Extract::HLGainRatioExtract(string inHLGainRatioFileName)
{
	cout<<"=========== HLGainRatioExtract ============"<<endl;
	fill_n(&HLRatio[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&HLIntercept[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 0);
	fill_n(&HighADCSwitch[0][0][0], LAYERNO*CHIPNO*CHANNO, 1);
    TFile* HLRatioFile = new TFile(TString(inHLGainRatioFileName),"READ");
    TTree* HLRatioTree = (TTree*)HLRatioFile->Get("InterCalib");
	if(!HLRatioTree)	cout<<"!!! GET HLRatio FILE FAILED !!!"<<endl;
	Int_t _HLRatioCellID;
	Double_t _slope, _slopeError, _intercept, _interceptError, _XMax;
    HLRatioTree->SetBranchAddress("CellID",&_HLRatioCellID);
    HLRatioTree->SetBranchAddress("Slope",&_slope);
    HLRatioTree->SetBranchAddress("SlopeError",&_slopeError);
    HLRatioTree->SetBranchAddress("Intercept",&_intercept);
    HLRatioTree->SetBranchAddress("InterceptError",&_interceptError);
    HLRatioTree->SetBranchAddress("XMax",&_XMax);
    
	double meanSlope  	 = 0; 
	double meanIntercept  	 = 0; 
    double chnNoSlope 	 = 0;
	double minHighSwitch = 4000;
    for(int entry=0; entry!=HLRatioTree->GetEntries(); entry++)
    {
        HLRatioTree->GetEntry(entry);
        int layerID = _HLRatioCellID/1e5;
		if(layerID>=LAYERNO)	continue;	
        int chipID  = (_HLRatioCellID%100000) /1e4;
        int chanID   = _HLRatioCellID%100;
		
		meanSlope += 1. / _slope;
        meanIntercept += _intercept;
		chnNoSlope++;
		minHighSwitch = ((_XMax-600) < minHighSwitch) ? (_XMax-600) : minHighSwitch;

        HLRatio[layerID][chipID][chanID] = 1. / _slope;
        HLIntercept[layerID][chipID][chanID] = _intercept;
		HighADCSwitch[layerID][chipID][chanID] = _XMax-600;
    }
    meanSlope /= chnNoSlope;
    meanIntercept /= chnNoSlope;
    cout << " channel count       : " << chnNoSlope << endl;
    cout << " mean of HLRatio     : " << meanSlope  << endl;
    cout << " mean of HLIntercept : " << meanIntercept  << endl;
    
	for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
					HLRatio[il][ip][ic] = (HLRatio[il][ip][ic]==1) ? meanSlope : HLRatio[il][ip][ic];
					// HLIntercept[il][ip][ic] = (HLIntercept[il][ip][ic]==0) ? meanIntercept : HLIntercept[il][ip][ic];
					HighADCSwitch[il][ip][ic] = (HighADCSwitch[il][ip][ic]==1) ? minHighSwitch : HighADCSwitch[il][ip][ic];
            }
		}
	}
	HLRatioFile->Close(); // jiaxuan
}

void Extract::MIPPeakValueExtract(string inMIPFileName)
{
    cout<<"=========== MIPExtract ============"<<endl;
	fill_n(&MIPPeakValue[0][0][0], 	LAYERNO*CHIPNO*CHANNO, 1);
    TFile *MIPFile = new TFile(TString(inMIPFileName),"READ");
    TTree* MIPTree = (TTree*)MIPFile->Get("MIP_Fit");
	if(!MIPTree)	{cout<<"!!! GET MIP FILE FAILED !!!"<<endl;}
	Int_t _mipCellID, _NDF;
	Double_t _landauMPV, _ChiSqr;
    MIPTree->SetBranchAddress("CellID",&_mipCellID);
    MIPTree->SetBranchAddress("LandauMPV",&_landauMPV);
    MIPTree->SetBranchAddress("ChiSquare",&_ChiSqr);
    MIPTree->SetBranchAddress("NDF",&_NDF);
    
	double meanMIP15  = 0 ; double meanMIP10  = 0;
    double chnNoMIP15 = 0 ; double chnNoMIP10 = 0;
    for(int entry=0; entry!=MIPTree->GetEntries(); entry++)
    {
        MIPTree->GetEntry(entry);
        int layerID = _mipCellID/1e5;
		if(layerID>=LAYERNO)	continue;	
        int chipID  = (_mipCellID%100000) /1e4;
        int chanID   = _mipCellID%100;
        if(_landauMPV<0) continue;
        if (layerID>3 && layerID<28) {
           // if(_ChiSqr/(double)_NDF>2) continue;
            chnNoMIP10++;
            meanMIP10 +=_landauMPV;
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
            chnNoMIP15++;
            meanMIP15 +=_landauMPV;
        }
        MIPPeakValue[layerID][chipID][chanID] = _landauMPV;
    }
    meanMIP10 /= chnNoMIP10;
    meanMIP15 /= chnNoMIP15;
    cout<<" channel count(15um:10um) : "<<chnNoMIP15<<" , "<<chnNoMIP10<<endl;
    cout<<" mean of MIPPeakValue : "<<meanMIP15<<" , "<<meanMIP10<<endl;
    
	for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
                if (il>3 && il<28) {
					MIPPeakValue[il][ip][ic] = (MIPPeakValue[il][ip][ic]==1) ? meanMIP10 : MIPPeakValue[il][ip][ic];
                } else {
					MIPPeakValue[il][ip][ic] = (MIPPeakValue[il][ip][ic]==1) ? meanMIP15 : MIPPeakValue[il][ip][ic];
                }
            }
		}
	}
	MIPFile->Close();
}

void Extract::SiPMGainExtract(string inSiPMGainFileName)
{
    cout<<"=========== SiPMGainExtract ============"<<endl;
	fill_n(&SiPMGain[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&SiPMGainError[0][0][0], LAYERNO*CHIPNO*CHANNO, 1);
    TFile *SiPMGainFile = new TFile(TString(inSiPMGainFileName),"READ");
    TTree* SiPMGainTree = (TTree*)SiPMGainFile->Get("T_Event");
	if(!SiPMGainTree)	{cout<<"!!! GET SiPM Gain FILE FAILED !!!"<<endl;}
	Int_t _gainCellID;
	Double_t _gain, _gainError;
    SiPMGainTree->SetBranchAddress("cellIDs",&_gainCellID);
    SiPMGainTree->SetBranchAddress("gain",&_gain);
    SiPMGainTree->SetBranchAddress("gainError",&_gainError);
    
	double meanGain15=0; double meanGain10=0;
	double meanGainError15=0; double meanGainError10=0;
    double chnNoGain15=0; double chnNoGain10=0;
	for(int entry=0; entry!=SiPMGainTree->GetEntries(); entry++)
    {
        SiPMGainTree->GetEntry(entry);
		// jiaxuan : notice tatsuki version layerID is cellID/1e6
        int layerID = _gainCellID/1e6;
		if(layerID>=LAYERNO)	continue;	
        int chipID  = (_gainCellID%1000000) /1e4;
        int chanID   = _gainCellID%100;
        if (layerID>3 && layerID<28) {
           // if(_ChiSqr/(double)_NDF>2) continue;
            chnNoGain10++;
            meanGain10 += _gain;
			meanGainError10 += _gainError;
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
            chnNoGain15++;
            meanGain15 += _gain;
			meanGainError15 += _gainError;
        }
        SiPMGain[layerID][chipID][chanID] = _gain;
		SiPMGainError[layerID][chipID][chanID] = _gainError;
    }
    meanGain10 /= chnNoGain10;
    meanGain15 /= chnNoGain15;
    meanGainError10 /= chnNoGain10;
    meanGainError15 /= chnNoGain15;
    cout << " channel count(15um:10um) : " << chnNoGain15 << " , " << chnNoGain10 << endl;
    cout << " mean of SiPM Gain : " << meanGain15 << " , " << meanGain10 << endl;
    cout << " mean of SiPM Gain Error: " << meanGainError15 << " , " << meanGainError10 << endl;
    
	for (int il=0; il<LAYERNO; il++) {
        for (int ip=0; ip<CHIPNO; ip++) {
            for (int ic=0; ic<CHANNO; ic++) {
                if (il>3 && il<28) {
					SiPMGain[il][ip][ic] = (SiPMGain[il][ip][ic]==1) ? meanGain10 : SiPMGain[il][ip][ic];
					SiPMGainError[il][ip][ic] = (SiPMGainError[il][ip][ic]==1) ? meanGainError10 : SiPMGainError[il][ip][ic];
                } else {
					SiPMGain[il][ip][ic] = (SiPMGain[il][ip][ic]==1) ? meanGain15 : SiPMGain[il][ip][ic];
					SiPMGainError[il][ip][ic] = (SiPMGainError[il][ip][ic]==1) ? meanGainError15 : SiPMGainError[il][ip][ic];
                } 
            }
		}
	}
	SiPMGainFile->Close();
}


void Extract::TresholdExtract(string inThresholdFileName, int thre)
{
    cout<<"=========== Threshold Extract ============"<<endl;
	fill_n(&thresholdValue[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 1);
	fill_n(&thresholdValue_with_pedestal[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 0);
	fill_n(&thresholdSigmaValue[0][0][0], 		LAYERNO*CHIPNO*CHANNO, 0);
	fill_n(&thresholdValue_chip[0][0], 		LAYERNO*CHIPNO, 1);
	fill_n(&thresholdSigmaValue_chip[0][0], 		LAYERNO*CHIPNO, 0);
	fill_n(&thresholdValue_chipCount[0][0], 		LAYERNO*CHIPNO, 0);
    TFile *ThresholdFile = new TFile(TString(inThresholdFileName),"READ");
    // TTree* ThresholdTree = (TTree*)ThresholdFile->Get("threshold");
    TTree* ThresholdTree = (TTree*)ThresholdFile->Get("expErfThre");
	if(!ThresholdTree){
        cout<<"!!! GET Threshold FILE FAILED !!!"<<endl;
        
        for(int layerID=0;layerID<LAYERNO;layerID++){
            int thres = layerID>3 && layerID<28 ? thre : 0;
            for(int chipID=0;chipID<CHIPNO;chipID++){
                for(int chanID=0;chanID<CHANNO;chanID++){
                    thresholdValue[layerID][chipID][chanID] = thres;
                    thresholdValue_with_pedestal[layerID][chipID][chanID] = pedMeanHigh[layerID][chipID][chanID] + thres;
                }
            }
        }
        cout << " channel count(15um:10um) : " << 0 << " , " << thre << endl;
    } else {
        Int_t _thresholdCellID;
    	Double_t _threshold, _thresholdSigma, _thoresold_with_pedestal;
        ThresholdTree->SetBranchAddress("CellID",&_thresholdCellID);
        ThresholdTree->SetBranchAddress("threshold",&_threshold);
        ThresholdTree->SetBranchAddress("sigma",&_thresholdSigma);
        // ThresholdTree->SetBranchAddress("thresold_with_pedestal",&_thoresold_with_pedestal);

        // int thres = thre;

    	double meanThreshold15=0; double meanThreshold10=0;
    	double sigmaThreshold15=0; double sigmaThreshold10=0;
        double chnNoThreshold15=0; double chnNoThreshold10=0;
    	for(int entry=0; entry!=ThresholdTree->GetEntries(); entry++)
        {
            ThresholdTree->GetEntry(entry);
    		// jiaxuan : notice tatsuki version layerID is cellID/1e6
            int layerID = _thresholdCellID/1e5;
    		if(layerID>=LAYERNO)	continue;	
            int chipID  = (_thresholdCellID%100000) /1e4;
            int chanID   = _thresholdCellID%100;
            // thres = layerID==14 ? 35 : thre;
            if (layerID>3 && layerID<28) {
               // if(_ChiSqr/(double)_NDF>2) continue;
                chnNoThreshold10++;
                meanThreshold10  += _threshold;
                sigmaThreshold10 += _thresholdSigma;
                // meanThreshold10 += _thoresold_with_pedestal - thres;
                // meanThreshold10 += _thoresold_with_pedestal -_threshold + thres;
            } else {
               // if(_ChiSqr/(double)_NDF>1.7) continue;
                chnNoThreshold15++;
                meanThreshold15  += _threshold;
                sigmaThreshold15 += _thresholdSigma;
                // meanThreshold15 += _thoresold_with_pedestal - thres;
                // meanThreshold15 += _thoresold_with_pedestal -_threshold + thres;
            }
            // if (_threshold>50) {
            thresholdValue_chip[layerID][chipID] += _threshold;
            thresholdSigmaValue_chip[layerID][chipID] += _thresholdSigma;
            thresholdValue_chipCount[layerID][chipID]++;
            // }
            thresholdValue[layerID][chipID][chanID] = _threshold;
            thresholdSigmaValue[layerID][chipID][chanID] = _thresholdSigma;
            // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal - thres;
            // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal -_threshold + thres;
        }
        meanThreshold10 /= chnNoThreshold10;
        meanThreshold15 /= chnNoThreshold15;
        sigmaThreshold10 /= chnNoThreshold10;
        sigmaThreshold15 /= chnNoThreshold15;
        cout << " channel count(15um:10um) : " << chnNoThreshold15 << " , " << chnNoThreshold10 << endl;
        cout << " mean of Threshold : " << meanThreshold15 << " , " << meanThreshold10 << endl;
        cout << " sigma of Threshold : " << sigmaThreshold15 << " , " << sigmaThreshold10 << endl;

    	for (int il=0; il<LAYERNO; il++) {
            for (int ip=0; ip<CHIPNO; ip++) {
                for (int ic=0; ic<CHANNO; ic++) {
                    if (il>3 && il<28) {
    					// thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? meanThreshold10 : thresholdValue[il][ip][ic];
    					thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? thresholdValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdValue[il][ip][ic];
    					thresholdSigmaValue[il][ip][ic] = (thresholdSigmaValue[il][ip][ic]==0) ? thresholdSigmaValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdSigmaValue[il][ip][ic];
    					// thresholdValue_with_pedestal[il][ip][ic] = (thresholdValue_with_pedestal[il][ip][ic]==1) ? meanThreshold10 : thresholdValue_with_pedestal[il][ip][ic];
                    } else {
    					// thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? meanThreshold15 : thresholdValue[il][ip][ic];
    					// thresholdValue[il][ip][ic] = 0;
                        thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? thresholdValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdValue[il][ip][ic];
    					thresholdSigmaValue[il][ip][ic] = (thresholdSigmaValue[il][ip][ic]==0) ? thresholdSigmaValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdSigmaValue[il][ip][ic];
                    } 
                }
    		}
    	}
    	ThresholdFile->Close();
    }
}

