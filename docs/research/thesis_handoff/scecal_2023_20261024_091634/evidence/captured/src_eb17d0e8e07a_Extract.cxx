#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>
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

namespace {
bool connectedChannel(int layer, int chip, int channel)
{
    return layer >= 0 && layer < Extract::LAYERNO && chip >= 0 && chip < Extract::CHIPNO &&
           channel >= 0 && channel < Extract::CHANNO && !(chip == 5 && channel > 29);
}

int calibrationGroup(int layer)
{
    return (layer > 3 && layer < 28) ? 1 : 0;
}

std::unique_ptr<TFile> openFactorFile(const string& name)
{
    std::unique_ptr<TFile> file(TFile::Open(name.c_str(), "READ"));
    if (!file || file->IsZombie())
        throw std::runtime_error("Cannot open calibration factor file: " + name);
    return file;
}

template<typename T>
void bindFactorBranch(TTree* tree, const char* name, T* address)
{
    if (tree->SetBranchAddress(name, address) < 0)
        throw std::runtime_error(string("Cannot read calibration branch: ") + name);
}
}

void Extract::PedestalExtract(string inPedFileName)
{
    cout << "=========== PedestalExtract ============" << endl;
    // Match processCali: missing entries stay zero; vector<float> reading also
    // matches its rounding when the source ROOT stores vector<double>.
    fill_n(&pedMeanHigh[0][0][0], LAYERNO*CHIPNO*CHANNO, 0.);
    fill_n(&pedSigmaHigh[0][0][0], LAYERNO*CHIPNO*CHANNO, 0.);
    fill_n(&pedMeanLow[0][0][0], LAYERNO*CHIPNO*CHANNO, 0.);
    fill_n(&pedSigmaLow[0][0][0], LAYERNO*CHIPNO*CHANNO, 0.);
    auto file = openFactorFile(inPedFileName);
    TTree* tree = nullptr;
    file->GetObject("ChnLevel", tree);
    if (!tree) throw std::runtime_error("Missing ChnLevel tree: " + inPedFileName);
    vector<int> idValues;
    vector<float> highValues, lowValues, highSigmaValues, lowSigmaValues;
    auto* ids = &idValues;
    auto* high = &highValues;
    auto* low = &lowValues;
    auto* highSigma = &highSigmaValues;
    auto* lowSigma = &lowSigmaValues;
    bindFactorBranch(tree, "CellID", &ids);
    bindFactorBranch(tree, "PedHighMean", &high);
    bindFactorBranch(tree, "PedLowMean", &low);
    bindFactorBranch(tree, "PedHighSig", &highSigma);
    bindFactorBranch(tree, "PedLowSig", &lowSigma);
    int loaded = 0;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        if (tree->GetEntry(entry) <= 0) throw std::runtime_error("Cannot read pedestal entry");
        if (ids->size() != high->size() || ids->size() != low->size() ||
            ids->size() != highSigma->size() || ids->size() != lowSigma->size())
            throw std::runtime_error("Pedestal vector size mismatch");
        for (size_t hit = 0; hit < ids->size(); ++hit) {
            int id = ids->at(hit), l = id / 100000, p = (id % 100000) / 10000, c = id % 100;
            if (id < 0 || !connectedChannel(l, p, c)) continue;
            if (!std::isfinite(high->at(hit)) || !std::isfinite(low->at(hit)) ||
                !std::isfinite(highSigma->at(hit)) || !std::isfinite(lowSigma->at(hit)) ||
                highSigma->at(hit) < 0. || lowSigma->at(hit) < 0.)
                throw std::runtime_error("Invalid pedestal for CellID " + std::to_string(id));
            pedMeanHigh[l][p][c] = high->at(hit);
            pedMeanLow[l][p][c] = low->at(hit);
            pedSigmaHigh[l][p][c] = highSigma->at(hit);
            pedSigmaLow[l][p][c] = lowSigma->at(hit);
            ++loaded;
        }
    }
    tree->ResetBranchAddresses();
    cout << " pedestal records loaded: " << loaded << "; missing entries use zero" << endl;
}

void Extract::HLGainRatioExtract(string inHLGainRatioFileName)
{
    cout << "=========== HLGainRatioExtract ============" << endl;
    auto file = openFactorFile(inHLGainRatioFileName);
    TTree* tree = nullptr;
    file->GetObject("InterCalib", tree);
    if (!tree) throw std::runtime_error("Missing InterCalib tree: " + inHLGainRatioFileName);
    int cellID = 0;
    double slope = 0., intercept = 0., xMax = 0.;
    bindFactorBranch(tree, "CellID", &cellID);
    bindFactorBranch(tree, "Slope", &slope);
    bindFactorBranch(tree, "Intercept", &intercept);
    bindFactorBranch(tree, "XMax", &xMax);
    bool valid[LAYERNO][CHIPNO][CHANNO] = {};
    bool switchValid[LAYERNO][CHIPNO][CHANNO] = {};
    double offsets[LAYERNO][CHIPNO][CHANNO] = {};
    fill_n(&HLRatio[0][0][0], LAYERNO*CHIPNO*CHANNO, 1.);
    fill_n(&HLIntercept[0][0][0], LAYERNO*CHIPNO*CHANNO, 0.);
    fill_n(&HighADCSwitch[0][0][0], LAYERNO*CHIPNO*CHANNO, 1.);
    int rejected = 0;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        if (tree->GetEntry(entry) <= 0) throw std::runtime_error("Cannot read HL factor entry");
        int l = cellID / 100000, p = (cellID % 100000) / 10000, c = cellID % 100;
        if (cellID < 0 || !connectedChannel(l, p, c) ||
            !std::isfinite(slope) || slope < 0.02 || slope > 0.05 || !std::isfinite(intercept)) {
            ++rejected;
            continue;
        }
        double factor = 1. / slope, offset = -intercept / slope;
        if (!std::isfinite(factor) || !std::isfinite(offset)) { ++rejected; continue; }
        HLRatio[l][p][c] = factor;
        HLIntercept[l][p][c] = intercept;
        offsets[l][p][c] = offset;
        valid[l][p][c] = true;
        switchValid[l][p][c] = std::isfinite(xMax) && xMax > 600.;
        if (switchValid[l][p][c]) HighADCSwitch[l][p][c] = xMax - 600.;
    }
    // Same connected-channel pool, unique-channel weighting and transformed
    // coefficients as InterFactorExtract_simulation in Calibration.cxx.
    long double sum[2] = {0., 0.}, offsetSum[2] = {0., 0.};
    int count[2] = {0, 0};
    double minHighSwitch = 4000.;
    for (int l = 0; l < LAYERNO; ++l)
        for (int p = 0; p < CHIPNO; ++p)
            for (int c = 0; c < CHANNO; ++c)
                if (connectedChannel(l, p, c) && valid[l][p][c]) {
                    int group = calibrationGroup(l);
                    sum[group] += HLRatio[l][p][c];
                    offsetSum[group] += offsets[l][p][c];
                    ++count[group];
                    if (switchValid[l][p][c]) minHighSwitch = std::min(minHighSwitch, HighADCSwitch[l][p][c]);
                }
    double mean[2], meanOffset[2];
    for (int group = 0; group < 2; ++group) {
        if (count[group] == 0) throw std::runtime_error("No valid HL factors in group " + std::to_string(group));
        mean[group] = static_cast<double>(sum[group] / count[group]);
        meanOffset[group] = static_cast<double>(offsetSum[group] / count[group]);
        if (!std::isfinite(mean[group]) || mean[group] <= 0. || !std::isfinite(meanOffset[group]))
            throw std::runtime_error("Invalid mean HL factors");
    }
    int replaced = 0;
    for (int l = 0; l < LAYERNO; ++l)
        for (int p = 0; p < CHIPNO; ++p)
            for (int c = 0; c < CHANNO; ++c) {
                if (!connectedChannel(l, p, c)) continue;
                if (!valid[l][p][c]) {
                    int group = calibrationGroup(l);
                    HLRatio[l][p][c] = mean[group];
                    // Generation uses LG = HG / gain + b; reconstruction uses
                    // HG = LG * gain + offset. Hence b = -offset / gain.
                    HLIntercept[l][p][c] = -meanOffset[group] / mean[group];
                    ++replaced;
                }
                if (!switchValid[l][p][c]) HighADCSwitch[l][p][c] = minHighSwitch;
            }
    cout << " mean of HLRatio (EBU15, EBU10): " << mean[0] << " , " << mean[1] << endl;
    cout << " mean of calibration HG offset (EBU15, EBU10): " << meanOffset[0] << " , " << meanOffset[1]
         << "; rejected records: " << rejected << "; supplemented channels: " << replaced << endl;
}

void Extract::MIPPeakValueExtract(string inMIPFileName)
{
    cout << "=========== MIPExtract ============" << endl;
    auto file = openFactorFile(inMIPFileName);
    TTree* tree = nullptr;
    file->GetObject("MIP_Fit", tree);
    if (!tree) throw std::runtime_error("Missing MIP_Fit tree: " + inMIPFileName);
    int cellID = 0, ndf = 0;
    double mpv = 0., chiSquare = 0.;
    bindFactorBranch(tree, "CellID", &cellID);
    bindFactorBranch(tree, "LandauMPV", &mpv);
    bindFactorBranch(tree, "ChiSquare", &chiSquare);
    bindFactorBranch(tree, "NDF", &ndf);
    bool valid[LAYERNO][CHIPNO][CHANNO] = {};
    fill_n(&MIPPeakValue[0][0][0], LAYERNO*CHIPNO*CHANNO, 1.);
    int rejected = 0;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        if (tree->GetEntry(entry) <= 0) throw std::runtime_error("Cannot read MIP factor entry");
        int l = cellID / 100000, p = (cellID % 100000) / 10000, c = cellID % 100;
        if (cellID < 0 || !connectedChannel(l, p, c)) { ++rejected; continue; }
        double limit = calibrationGroup(l) == 1 ? 2. : 1.7;
        if (!std::isfinite(mpv) || mpv <= 0. || ndf <= 0 ||
            !std::isfinite(chiSquare) || chiSquare < 0. || chiSquare / ndf > limit) {
            ++rejected;
            continue;
        }
        MIPPeakValue[l][p][c] = mpv;
        valid[l][p][c] = true;
    }
    long double sum[2] = {0., 0.};
    int count[2] = {0, 0};
    for (int l = 0; l < LAYERNO; ++l)
        for (int p = 0; p < CHIPNO; ++p)
            for (int c = 0; c < CHANNO; ++c)
                if (connectedChannel(l, p, c) && valid[l][p][c]) {
                    int group = calibrationGroup(l);
                    sum[group] += MIPPeakValue[l][p][c];
                    ++count[group];
                }
    double mean[2];
    for (int group = 0; group < 2; ++group) {
        if (count[group] == 0) throw std::runtime_error("No valid MIP factors in group " + std::to_string(group));
        mean[group] = static_cast<double>(sum[group] / count[group]);
        if (!std::isfinite(mean[group]) || mean[group] <= 0.)
            throw std::runtime_error("Invalid mean MIP factor");
    }
    int replaced = 0;
    for (int l = 0; l < LAYERNO; ++l)
        for (int p = 0; p < CHIPNO; ++p)
            for (int c = 0; c < CHANNO; ++c)
                if (connectedChannel(l, p, c) && !valid[l][p][c]) {
                    MIPPeakValue[l][p][c] = mean[calibrationGroup(l)];
                    ++replaced;
                }
    cout << " mean of MIPPeakValue (EBU15, EBU10): " << mean[0] << " , " << mean[1]
         << "; accepted channels: " << count[0] << " , " << count[1]
         << "; rejected records: " << rejected << "; supplemented channels: " << replaced << endl;
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
	fill_n(&thresholdValue_chip[0][0], 		LAYERNO*CHIPNO, 0);
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

