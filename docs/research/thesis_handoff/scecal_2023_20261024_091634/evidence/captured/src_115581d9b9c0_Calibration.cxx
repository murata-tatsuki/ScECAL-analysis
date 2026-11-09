#include<fstream>
#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TRandom.h>
#include<cmath>
#include<memory>
#include<stdexcept>
#include"Calibration.h"
#include"EBUdecode.h"

namespace {
bool connectedChannel(int layer, int chip, int channel)
{
    return layer >= 0 && layer < 32 && chip >= 0 && chip < 6 &&
           channel >= 0 && channel < 36 && !(chip == 5 && channel > 29);
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

void processCali::PedestalExtract(string pedestal_name)
{
    cout<<"=========== PedestalExtract ============"<<endl;
    // initial the branch of pedestal fitting
    _pedCellID=0;
    _pedMeanCharges=0;
    _pedSigmaCharges=0;

    // *********************************************************************************
    //          extract pedestal fit parameters
    // *********************************************************************************
    fitfile = new TFile(pedestal_name.c_str());
    fitfile->GetObject("ChnLevel",fNtuple1); // channel level fit parameter
    //fitfile->GetObject("ChipLevel",fNtuple1); // chip level fit parameter
    //fitfile->GetObject("CellLevel",fNtuple1); // cell level fit parameter
    if(!fitfile){ cout<<" Get pedestal fit file failed !!!"<<endl;return;}
    fNtuple1 ->SetBranchAddress("cellIDs",&_pedCellID);
    fNtuple1 ->SetBranchAddress("pedMeanCharges",&_pedMeanCharges);
    fNtuple1 ->SetBranchAddress("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple1 ->SetBranchAddress("pedMeanTimes",&_pedMeanTimes);
    fNtuple1 ->SetBranchAddress("pedSigmaTimes",&_pedSigmaTimes);

    for(int entry=0; entry!=fNtuple1->GetEntries(); ++entry)
    {
        fNtuple1->GetEntry(entry);
        //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(unsigned ihit=0; ihit!=_pedCellID->size(); ++ihit)
        {
            int layerID = _pedCellID->at(ihit)/1e6;
            int chipIDs = _pedCellID->at(ihit)%1000000/1e4;
          //  int memoryCells = _pedCellID->at(ihit)%10000/1e2;
            int channelIDs = _pedCellID->at(ihit)%100;
            pedMeanChargesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanCharges->at(ihit);
            pedSigmaChargesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaCharges->at(ihit);
            pedMeanTimesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanTimes->at(ihit);
            pedSigmaTimesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaTimes->at(ihit);
        }
    }
}

void processCali::MIPFactorExtract(string interfactor_name)
{
    cout<<"=========== MIPFactorExtract ============"<<endl;

    // ******************************************************************************
    // *********************************************************************************
    int _cellID, _ndf;
    double _landauMPV, _landauWidth,  _gausSigma, _ChiSqr;
    TFile *infile = new TFile(interfactor_name.c_str());
    TTree* intree = (TTree*)infile->Get("T_Event");
    intree->SetBranchAddress("cellIDs",&_cellID);
    intree->SetBranchAddress("landauMPV",&_landauMPV);
    intree->SetBranchAddress("landauWidth",&_landauWidth);
    intree->SetBranchAddress("GauSigma",&_gausSigma);
    intree->SetBranchAddress("ChiSqr",&_ChiSqr);
    intree->SetBranchAddress("ndf",&_ndf);

    double meanEBU15=0; double meanEBU10=0;
    double chnEBU15=0; double chnEBU10=0;
    for(int entry=0; entry!=intree->GetEntries(); entry++)
    {
        intree->GetEntry(entry);
        int _layerID = _cellID/1e6;
        int _chipID = _cellID%1000000/1e4;
        int _chnID = _cellID%10000%100;
        if(_landauMPV<0) continue;
        
        if(_layerID>3 && _layerID<28) {
            if(_ChiSqr/(double)_ndf>2) continue;
            chnEBU10++;
            meanEBU10 +=_landauMPV;
        }
        else if(_layerID<30) {
            if(_ChiSqr/(double)_ndf>1.7) continue;
            chnEBU15++;
            meanEBU15 +=_landauMPV;
        }
        mipFactor[_layerID][_chipID][_chnID] = _landauMPV;
    }

    meanEBU10 = meanEBU10/chnEBU10;
    meanEBU15 = meanEBU15/chnEBU15;

    cout<<" mean of mipFactor : "<<meanEBU15<<" , "<<meanEBU10<<endl;

    for(int il=0; il<layerNu; il++){
        for(int ip=0; ip<chipNu; ip++){
            for(int ic=0; ic<channelNu; ic++){
                if(il>3 && il<28) {
                    if(mipFactor[il][ip][ic]==1) mipFactor[il][ip][ic] = meanEBU10;
                }
                else{
                    if(mipFactor[il][ip][ic]==1) mipFactor[il][ip][ic] = meanEBU15;
                }
            }}}
}

void processCali::InterFactorExtract(string interfactor_name)
{

    cout<<"=========== InterFactorExtract ============"<<endl;
    // *********************************************************************************
    // *********************************************************************************
    int _cellID;
    double _interFactor;
    TFile *infile = new TFile(interfactor_name.c_str());
    TTree* intree = (TTree*)infile->Get("T");
    intree->SetBranchAddress("cellID",&_cellID);
    intree->SetBranchAddress("slope",&_interFactor);

    double meanEBU15=0; double meanEBU10=0;
    double chnEBU15=0; double chnEBU10=0;
    for(int entry=0; entry!=intree->GetEntries(); entry++)
    {
        intree->GetEntry(entry);
        int _layerID = entry/EBUchn;
        int _chipID = entry%EBUchn/36;
        int _chnID = entry%EBUchn%36;
       if(_interFactor<0.02 || _interFactor>0.05) continue;

       interFactor[_layerID][_chipID][_chnID] = 1./_interFactor;
        
        if(_layerID>3 && _layerID<28) {
            chnEBU10++;
            meanEBU10 += 1./_interFactor;
        }
        else if(_layerID<30) {
            chnEBU15++;
            meanEBU15 += 1./_interFactor;
        }
    }

    meanEBU10 = meanEBU10/chnEBU10;
    meanEBU15 = meanEBU15/chnEBU15;

    cout<<" mean of interFactor : "<<meanEBU15<<" , "<<meanEBU10<<endl;

    for(int il=0; il<layerNu; il++){
        for(int ip=0; ip<chipNu; ip++){
            for(int ic=0; ic<channelNu; ic++){
                if(il>3 && il<28) {
                    if(interFactor[il][ip][ic]==1) interFactor[il][ip][ic] = meanEBU10;
                }
                else{
                    if(interFactor[il][ip][ic]==1) interFactor[il][ip][ic] = meanEBU15;
                }
            }}}

}

void processCali::PedestalExtract_edit(string pedestal_name)
{
    cout<<"=========== PedestalExtract ============"<<endl;
    // initial the branch of pedestal fitting
    _pedCellID=0;
    _pedMeanCharges=0;
    _pedSigmaCharges=0;

    // *********************************************************************************
    //          extract pedestal fit parameters
    // *********************************************************************************
    fitfile = new TFile(pedestal_name.c_str());
    fitfile->GetObject("ChnLevel",fNtuple1); // channel level fit parameter
    //fitfile->GetObject("ChipLevel",fNtuple1); // chip level fit parameter
    //fitfile->GetObject("CellLevel",fNtuple1); // cell level fit parameter
    if(!fitfile){ cout<<" Get pedestal fit file failed !!!"<<endl;return;}
    fNtuple1 ->SetBranchAddress("CellID",&_pedCellID);
    fNtuple1 ->SetBranchAddress("PedLowMean",&_pedMeanCharges);
    fNtuple1 ->SetBranchAddress("PedLowSig",&_pedSigmaCharges);
    fNtuple1 ->SetBranchAddress("PedHighMean",&_pedMeanTimes);
    fNtuple1 ->SetBranchAddress("PedHighSig",&_pedSigmaTimes);

    for(int entry=0; entry!=fNtuple1->GetEntries(); ++entry)
    {
        fNtuple1->GetEntry(entry);
        //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(unsigned ihit=0; ihit!=_pedCellID->size(); ++ihit)
        {
            int layerID = _pedCellID->at(ihit)/1e5;
            int chipIDs = _pedCellID->at(ihit)%100000/1e4;
          //  int memoryCells = _pedCellID->at(ihit)%10000/1e2;
            int channelIDs = _pedCellID->at(ihit)%100;
            pedMeanChargesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanCharges->at(ihit);
            pedSigmaChargesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaCharges->at(ihit);
            pedMeanTimesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanTimes->at(ihit);
            pedSigmaTimesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaTimes->at(ihit);
        }
    }
}

void processCali::MIPFactorExtract_edit(string interfactor_name)
{
    cout<<"=========== MIPFactorExtract ============"<<endl;

    // ******************************************************************************
    // *********************************************************************************
    int _cellID, _ndf;
    double _landauMPV, _landauWidth,  _gausSigma, _ChiSqr, _HG_pedMean;
    TFile *infile = new TFile(interfactor_name.c_str());
    TTree* intree = (TTree*)infile->Get("T_Event");
    intree->SetBranchAddress("cellIDs",&_cellID);
    intree->SetBranchAddress("landauMPV",&_landauMPV);
    // intree->SetBranchAddress("landauMPV_jiaxuan",&_landauMPV);
    intree->SetBranchAddress("landauWidth",&_landauWidth);
    intree->SetBranchAddress("HG_pedMean",&_HG_pedMean);    // edition
    intree->SetBranchAddress("gausSigma",&_gausSigma);
    intree->SetBranchAddress("ChiSqr",&_ChiSqr);
    intree->SetBranchAddress("ndf",&_ndf);

    double meanEBU15=0; double meanEBU10=0;
    double chnEBU15=0; double chnEBU10=0;
    for(int entry=0; entry!=intree->GetEntries(); entry++)
    {
        intree->GetEntry(entry);
        int _layerID = _cellID/1e6;
        int _chipID = _cellID%1000000/1e4;
        int _chnID = _cellID%10000%100;

        // double __landauMPV = _landauMPV;
        double __landauMPV = _landauMPV - _HG_pedMean;
        if(__landauMPV<50) continue;
        
        if(_layerID>3 && _layerID<28) {
            if(_ChiSqr/(double)_ndf>2) continue;
            chnEBU10++;
            meanEBU10 += __landauMPV;
        }
        else if(_layerID<30) {
            if(_ChiSqr/(double)_ndf>1.7) continue;
            chnEBU15++;
            meanEBU15 += __landauMPV;
        }
        mipFactor[_layerID][_chipID][_chnID] = __landauMPV;
    }

    meanEBU10 = meanEBU10/chnEBU10;
    meanEBU15 = meanEBU15/chnEBU15;

    cout<<" mean of mipFactor : "<<meanEBU15<<" , "<<meanEBU10<<endl;

    for(int il=0; il<layerNu; il++){
        for(int ip=0; ip<chipNu; ip++){
            for(int ic=0; ic<channelNu; ic++){
                if(il>3 && il<28) {
                    if(mipFactor[il][ip][ic]==1) mipFactor[il][ip][ic] = meanEBU10;
                }
                else{
                    if(mipFactor[il][ip][ic]==1) mipFactor[il][ip][ic] = meanEBU15;
                }
            }}}
}

void processCali::InterFactorExtract_edit(string interfactor_name)
{

    cout<<"=========== InterFactorExtract ============"<<endl;
    // *********************************************************************************
    // *********************************************************************************
    int _cellID;
    double _interFactor, _intercept;
    TFile *infile = new TFile(interfactor_name.c_str());
    TTree* intree = (TTree*)infile->Get("InterCalib");
    intree->SetBranchAddress("CellID",&_cellID);
    intree->SetBranchAddress("Slope",&_interFactor);
    intree->SetBranchAddress("Intercept",&_intercept);

    double meanEBU15=0; double meanEBU10=0;
    double meanInterceptEBU15=0; double meanInterceptEBU10=0;
    double chnEBU15=0; double chnEBU10=0;
    for(int entry=0; entry!=intree->GetEntries(); entry++)
    {
        intree->GetEntry(entry);
        int _layerID = _cellID/1e5;
        int _chipID = _cellID%1000000/1e4;
        int _chnID = _cellID%10000%100;
       if(_interFactor<0.02 || _interFactor>0.05) continue;

       interFactor[_layerID][_chipID][_chnID] = 1./_interFactor;
       interFactorIntercept[_layerID][_chipID][_chnID] = -_intercept/_interFactor;
        
        if(_layerID>3 && _layerID<28) {
            chnEBU10++;
            meanEBU10 += 1./_interFactor;
            meanInterceptEBU10 += -_intercept/_interFactor;
        }
        else if(_layerID<30) {
            chnEBU15++;
            meanEBU15 += 1./_interFactor;
            meanInterceptEBU15 += -_intercept/_interFactor;
        }
    }

    meanEBU10 = meanEBU10/chnEBU10;
    meanEBU15 = meanEBU15/chnEBU15;
    meanInterceptEBU10 = meanInterceptEBU10/chnEBU10;
    meanInterceptEBU15 = meanInterceptEBU15/chnEBU15;

    cout<<" mean of interFactor : "<<meanEBU15<<" , "<<meanEBU10<<endl;
    cout<<" mean of intercept   : "<<meanInterceptEBU15<<" , "<<meanInterceptEBU10<<endl;

    for(int il=0; il<layerNu; il++){
        for(int ip=0; ip<chipNu; ip++){
            for(int ic=0; ic<channelNu; ic++){
                if(il>3 && il<28) {
                    if(interFactor[il][ip][ic]==1) interFactor[il][ip][ic] = meanEBU10;
                    if(interFactorIntercept[il][ip][ic]==1) interFactorIntercept[il][ip][ic] = meanInterceptEBU10;
                }
                else{
                    if(interFactor[il][ip][ic]==1) interFactor[il][ip][ic] = meanEBU15;
                    if(interFactorIntercept[il][ip][ic]==1) interFactorIntercept[il][ip][ic] = meanInterceptEBU15;
                }
            }}}

}

void processCali::PedestalExtract_simulation(string pedestal_name)
{
    cout<<"=========== PedestalExtract ============"<<endl;
    // initial the branch of pedestal fitting
    _pedCellID=0;
    _pedMeanCharges=0;
    _pedSigmaCharges=0;
    _pedMeanTimes=0;
    _pedSigmaTimes=0;

    // *********************************************************************************
    //          extract pedestal fit parameters
    // *********************************************************************************
    fitfile = new TFile(pedestal_name.c_str());
    fitfile->GetObject("ChnLevel",fNtuple1); // channel level fit parameter
    //fitfile->GetObject("ChipLevel",fNtuple1); // chip level fit parameter
    //fitfile->GetObject("CellLevel",fNtuple1); // cell level fit parameter
    if(!fitfile){ cout<<" Get pedestal fit file failed !!!"<<endl;return;}
    fNtuple1 ->SetBranchAddress("CellID",&_pedCellID);
    fNtuple1 ->SetBranchAddress("PedLowMean",&_pedMeanCharges);
    fNtuple1 ->SetBranchAddress("PedLowSig",&_pedSigmaCharges);
    fNtuple1 ->SetBranchAddress("PedHighMean",&_pedMeanTimes);
    fNtuple1 ->SetBranchAddress("PedHighSig",&_pedSigmaTimes);

    for(int entry=0; entry!=fNtuple1->GetEntries(); ++entry)
    {
        fNtuple1->GetEntry(entry);
        // cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(unsigned ihit=0; ihit!=_pedCellID->size(); ++ihit)
        {
            int layerID = _pedCellID->at(ihit)/1e5;
            int chipIDs = _pedCellID->at(ihit)%100000/1e4;
          //  int memoryCells = _pedCellID->at(ihit)%10000/1e2;
            int channelIDs = _pedCellID->at(ihit)%100;
            pedMeanChargesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanCharges->at(ihit);
            pedSigmaChargesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaCharges->at(ihit);
            pedMeanTimesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanTimes->at(ihit);
            pedSigmaTimesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaTimes->at(ihit);
        }
    }
}

void processCali::MIPFactorExtract_simulation(string interfactor_name)
{
    cout << "=========== MIPFactorExtract ============" << endl;
    auto infile = openFactorFile(interfactor_name);
    TTree* intree = nullptr;
    infile->GetObject("MIP_Fit", intree);
    if (!intree) throw std::runtime_error("Missing MIP_Fit tree: " + interfactor_name);
    int cellID = 0, ndf = 0;
    double mpv = 0., chiSquare = 0.;
    bindFactorBranch(intree, "CellID", &cellID);
    bindFactorBranch(intree, "LandauMPV", &mpv);
    bindFactorBranch(intree, "ChiSquare", &chiSquare);
    bindFactorBranch(intree, "NDF", &ndf);

    for (int l = 0; l < layerNu; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c) mipFactorValid[l][p][c] = false;

    int rejected = 0;
    for (Long64_t entry = 0; entry < intree->GetEntries(); ++entry) {
        if (intree->GetEntry(entry) <= 0) throw std::runtime_error("Cannot read MIP factor entry");
        int l = cellID / 100000, p = (cellID % 100000) / 10000, c = cellID % 100;
        if (cellID < 0 || !connectedChannel(l, p, c)) { ++rejected; continue; }
        double limit = calibrationGroup(l) == 1 ? 2. : 1.7;
        if (!std::isfinite(mpv) || mpv <= 0. || ndf <= 0 ||
            !std::isfinite(chiSquare) || chiSquare < 0. || chiSquare / ndf > limit) {
            ++rejected;
            continue;
        }
        mipFactor[l][p][c] = mpv;
        mipFactorValid[l][p][c] = true;
    }

    // Use each accepted physical channel once; layers 30/31 do not set the means.
    long double sum[2] = {0., 0.};
    int count[2] = {0, 0};
    for (int l = 0; l < 30; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c)
                if (connectedChannel(l, p, c) && mipFactorValid[l][p][c]) {
                    int group = calibrationGroup(l);
                    sum[group] += mipFactor[l][p][c];
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
    for (int l = 0; l < layerNu; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c)
                if (connectedChannel(l, p, c) && !mipFactorValid[l][p][c]) {
                    mipFactor[l][p][c] = mean[calibrationGroup(l)];
                    ++replaced;
                }
    cout << " mean of mipFactor (EBU15, EBU10): " << mean[0] << " , " << mean[1]
         << "; accepted channels: " << count[0] << " , " << count[1]
         << "; rejected records: " << rejected << "; supplemented channels: " << replaced << endl;
}

void processCali::InterFactorExtract_simulation(string interfactor_name)
{
    cout << "=========== InterFactorExtract ============" << endl;
    auto infile = openFactorFile(interfactor_name);
    TTree* intree = nullptr;
    infile->GetObject("InterCalib", intree);
    if (!intree) throw std::runtime_error("Missing InterCalib tree: " + interfactor_name);
    int cellID = 0;
    double slope = 0., intercept = 0.;
    bindFactorBranch(intree, "CellID", &cellID);
    bindFactorBranch(intree, "Slope", &slope);
    bindFactorBranch(intree, "Intercept", &intercept);

    for (int l = 0; l < layerNu; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c) interFactorValid[l][p][c] = false;

    int rejected = 0;
    for (Long64_t entry = 0; entry < intree->GetEntries(); ++entry) {
        if (intree->GetEntry(entry) <= 0) throw std::runtime_error("Cannot read HL factor entry");
        int l = cellID / 100000, p = (cellID % 100000) / 10000, c = cellID % 100;
        if (cellID < 0 || !connectedChannel(l, p, c) ||
            !std::isfinite(slope) || slope < 0.02 || slope > 0.05 || !std::isfinite(intercept)) {
            ++rejected;
            continue;
        }
        double factor = 1. / slope, offset = -intercept / slope;
        if (!std::isfinite(factor) || !std::isfinite(offset)) { ++rejected; continue; }
        interFactor[l][p][c] = factor;
        interFactorIntercept[l][p][c] = offset;
        interFactorValid[l][p][c] = true;
    }

    long double sum[2] = {0., 0.}, offsetSum[2] = {0., 0.};
    int count[2] = {0, 0};
    for (int l = 0; l < 30; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c)
                if (connectedChannel(l, p, c) && interFactorValid[l][p][c]) {
                    int group = calibrationGroup(l);
                    sum[group] += interFactor[l][p][c];
                    offsetSum[group] += interFactorIntercept[l][p][c];
                    ++count[group];
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
    for (int l = 0; l < layerNu; ++l)
        for (int p = 0; p < chipNu; ++p)
            for (int c = 0; c < channelNu; ++c)
                if (connectedChannel(l, p, c) && !interFactorValid[l][p][c]) {
                    int group = calibrationGroup(l);
                    interFactor[l][p][c] = mean[group];
                    interFactorIntercept[l][p][c] = meanOffset[group];
                    ++replaced;
                }
    cout << " mean of interFactor (EBU15, EBU10): " << mean[0] << " , " << mean[1] << endl;
    cout << " mean of intercept (EBU15, EBU10): " << meanOffset[0] << " , " << meanOffset[1]
         << "; rejected records: " << rejected << "; supplemented channels: " << replaced << endl;
}

void processCali::Init()
{
    for(int il=0; il<layerNu; il++){
        for(int ip=0; ip<chipNu; ip++){
            for(int ic=0; ic<channelNu; ic++){
                mipFactor[il][ip][ic]=1;
                interFactor[il][ip][ic]=1;
                interFactorIntercept[il][ip][ic]=0;
                mipFactorValid[il][ip][ic]=false;
                interFactorValid[il][ip][ic]=false;
            }}}
}

void processCali::Reset()
{	
	_cellID_new.clear();
	//_hit_time_newIclear();
    //_lg_cali.clear();
    //_hg_cali.clear();
	_hitE.clear();
	_hit_hg_E.clear();
	_hit_lg_E.clear();
    _temp.clear();
	_posx.clear();
	_posy.clear();
	_posz.clear();
}

void processCali::Calibration(string root_name,string pedestal_name,string interfactor_name, string mipfactor_name, string cali_name)
{
    //string txtname(root_name);
    //txtname.erase(0,txtname.find_last_of("/")+1);
    string caliname(cali_name);
    //caliname.insert(caliname.length(), txtname);
    //caliname.insert(caliname.length(), ".root");

    Init();
    //pedestal_name.insert(pedestal_name.length(),".root");
    //cout<<pedestal_name<<endl;

    // PedestalExtract(pedestal_name);
    // InterFactorExtract(interfactor_name);
    // MIPFactorExtract(mipfactor_name);

    // PedestalExtract_edit(pedestal_name);
    // InterFactorExtract_edit(interfactor_name);
    // MIPFactorExtract_edit(mipfactor_name);

    PedestalExtract_simulation(pedestal_name);
    InterFactorExtract_simulation(interfactor_name);
    MIPFactorExtract_simulation(mipfactor_name);

    
    // *********************************************************************************
    string dataname(root_name);
    //dataname.insert(dataname.length(), ".root");
    cout << dataname << endl;
    rawfile = TFile::Open(dataname.c_str(),"READ");
    rawfile->GetObject("Raw_Hit",fNtuple2);
    if(!fNtuple2){ cout<<" Get EventAlign file failed !!!"<<endl;return;}
    int Entries = fNtuple2->GetEntries();
    cout<<" [ Prepare to processCalib : ] Total Entries : "<<Entries<<" "<<endl;
    _cellID = nullptr;
    _hitTag = nullptr;
    _gainTag = nullptr;
    _hg_charge = nullptr;
    _lg_charge = nullptr;
    _hit_time = nullptr;
 	fNtuple2 ->SetBranchAddress("Run_Num",&_run_num);
 	fNtuple2 ->SetBranchAddress("Event_Time",&_event_time);
   	fNtuple2 ->SetBranchAddress("TriggerID",&_triggerID);
    fNtuple2 ->SetBranchAddress("CellID",&_cellID);
   	//fNtuple2 ->SetBranchAddress("BCID",&_bcid);
    fNtuple2 ->SetBranchAddress("HitTag",&_hitTag);
   	fNtuple2 ->SetBranchAddress("GainTag",&_gainTag);
    fNtuple2 ->SetBranchAddress("HG_Charge",&_hg_charge);
    fNtuple2 ->SetBranchAddress("LG_Charge",&_lg_charge);
    fNtuple2 ->SetBranchAddress("Hit_Time",&_hit_time);
    fNtuple2 ->SetBranchAddress("Temperature",&_tempLayer);
 	
    // *********************************************************************************
    califile = new TFile(caliname.c_str(),"RECREATE");
    tree = new TTree("Calib_Hit","after calibration");	
	tree ->Branch("Run_Num",&_run_num_new);
	tree ->Branch("Event_Time",&_event_time_new);
	tree ->Branch("Event_Num",&_event_num);
	tree ->Branch("DetectorID",&_detectorID);
	tree ->Branch("CellID",&_cellID_new);
    tree ->Branch("Hit_Energy",&_hitE);
    tree ->Branch("Hit_HG_Energy",&_hit_hg_E);
    tree ->Branch("Hit_LG_Energy",&_hit_lg_E);
	tree ->Branch("Hit_X",&_posx);
	tree ->Branch("Hit_Y",&_posy);
	tree ->Branch("Hit_Z",&_posz);
    //tree ->Branch("Hit_Time",&_hit_time_new);
    //tree ->Branch("HG_Cali",&_hg_cali);
    //tree ->Branch("LG_Cali",&_lg_cali);
    tree ->Branch("NewTemperature",&_temp);
    tree ->Branch("TotalEnergyDep",&_tot_energy);
	//tree ->Branch("TotalHitMap",&_tot_hitlayer);
    
    //int Entries_1=0; 
    //int Entries_2=0;
    // ==============================================================================
    //                          Calibration signal 
    // ==============================================================================
	TH1F *Edep = new TH1F("Edep","Edep(MeV)",25000,0,10000);
	TH1I *HitLayerNo = new TH1I("HitLayerNo","HitLayerNo",32,0,32);
	TH2I *HitMap = new TH2I("HitMap","HitMap",32,0,32,100,0,100);
	for(int entry=0; entry<Entries; entry++)
	//for(int entry=0; entry<100; entry++)
    {
        fNtuple2->GetEntry(entry);
        if(entry%1000==0)   
            cout<<" [ Process : Calibration : ]  event : "<<entry<<" "<<endl;
		double EnergyDep = 0;
		int hitNuPerLayer[32]={0};//hit channel Number: on each layer
		int hitFlagForLayer[32]={0};//whether this layer has hit
		int hitLayerNoPerEvent = 0;
        Reset();
		_run_num_new = _run_num;
		_event_time_new = _event_time;
		_event_num = _triggerID;
		_detectorID = 0;
        for(unsigned ihit=0; ihit!=_hitTag->size(); ++ihit)
        {
            
            if(_hitTag->at(ihit)==0) continue; // only save hit=1 events, above threshold for auto-trigger
	//		cout <<"hitTag = " << _hitTag->at(ihit)  << endl;
	//		int _cellidd		= _cellID->at(ihit);
            int _layerID        = _cellID->at(ihit)/1e5;
			hitNuPerLayer[_layerID]++;
			hitFlagForLayer[_layerID] = 1;
            int _chipID         = (_cellID->at(ihit) - _layerID*1e5 )/1e4;
           // int memoryCells = _cellID->at(ihit)%10000/1e2;
            int _chnID          = _cellID->at(ihit)%100;
            double caliCharge   = _lg_charge->at(ihit); // for auto gain setting
            double caliTime     = _hg_charge->at(ihit); // for High&Low gain setting, High gain

            if(_chipID==5 && _chnID>29) continue;

            // //subtract the pedestal for each channel 
            caliCharge -= pedMeanChargesChnLevel[_layerID][_chipID][_chnID];
            caliTime   -= pedMeanTimesChnLevel[_layerID][_chipID][_chnID];

            // // reconstruction temperature distribution at each SiPM location
            double *_position = EBUdecode(_layerID,_chipID,_chnID);
            double SiPMtemp;
            if(_tempLayer->at(_layerID).size()!=0) SiPMtemp = tempReconstruction(_layerID, _position, _tempLayer->at(_layerID));
            else SiPMtemp=20;
			
            // // temperature depedence correction
            // double tempCoefficient = 0.03;
            double tempCoefficient = (_layerID<4 || _layerID>=28) ? 3.5/230 : 1.6/135;
            caliCharge = (1-(SiPMtemp-20)*tempCoefficient)*caliCharge;
            caliTime = (1-(SiPMtemp-20)*tempCoefficient)*caliTime;
            if(caliTime<10) continue; // exclude bug of SP2E

            // //inter factor calibration for Gain Mode (high gain & low gain)
            caliCharge = caliCharge*interFactor[_layerID][_chipID][_chnID]+interFactorIntercept[_layerID][_chipID][_chnID];
            // //MIP response calibration for Gain Mode
            caliCharge = caliCharge/mipFactor[_layerID][_chipID][_chnID];
            caliTime   = caliTime/mipFactor[_layerID][_chipID][_chnID];
		
			_cellID_new.push_back(_cellID->at(ihit));
			_hitE.push_back(_hg_charge->at(ihit)<2600? caliTime*0.305: caliCharge*0.305);//MeV
			_hit_hg_E.push_back(caliTime*0.305);//MeV
			_hit_lg_E.push_back(caliCharge*0.305);//MeV
			_posx.push_back(*_position);
			_posy.push_back(*(_position+1));
			_posz.push_back(*(_position+2));
            //_lg_cali.push_back(caliCharge);
            //_hg_cali.push_back(caliTime);
            _temp.push_back(SiPMtemp);
			//if(!(ihit%500)) cout <<" Hit : "<< ihit << endl;

			//if(_time->at(ihit)>2000) cout << "_time->at("<<ihit<<"):"<< _time->at(ihit);
			if(_layerID<30) EnergyDep += _hg_charge->at(ihit)<2600 ? caliTime*0.305 : caliCharge*0.305 ;
        }
		_tot_energy = EnergyDep;
		for(int i_layer=0; i_layer<32;i_layer++)	
		{	
			HitMap->Fill(i_layer,hitNuPerLayer[i_layer]);
			hitLayerNoPerEvent += hitFlagForLayer[i_layer];
		}
		//_tot_hitlayer = hitLayerNoPerEvent;
		HitLayerNo->Fill(hitLayerNoPerEvent);
    	tree->Fill();
		Edep->Fill(EnergyDep);
    }

    //cout<<"Total Entries : "<<Entries<<" , Layer>=22 : "<<Entries_1<<" , Strips<=64 : "<<Entries_2<<endl;
    // *********************************************************************************
    califile->cd();
	Edep->Write();
	HitMap->Write();
	HitLayerNo->Write();
    tree->Write("",TObject::kOverwrite);
    califile->Close();
    rawfile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;

   // **********************************************************************************
}

