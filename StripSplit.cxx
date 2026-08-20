#include<fstream>
#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TMath.h>
#include<numeric>
#include"StripSplit.h"

using namespace std;

void treeStripSplit::Reset()
{
    _cycleID_cali=0;
    _cellID_cali.clear();
    _bcid_cali.clear();
    _hitTag_cali.clear();
    _gainTag_cali.clear();
    _charge_cali.clear();
    _time_cali.clear();
    _temp_cali.clear();
    _posX_cali.clear();
    _posY_cali.clear();
    _posZ_cali.clear();
    _low_cali.clear();
    _high_cali.clear();
    _ssaTag.clear();
}

int treeStripSplit::StripSplit(string txt_name,string cali_name,string ssa_name)
{
    cali_name.insert(cali_name.length(), txt_name);
    ssa_name.insert(ssa_name.length(),txt_name);

cout << cali_name << endl;
cout << ssa_name << endl; 

    // *********************************************************************************
    califile = TFile::Open(cali_name.c_str(),"READ");
    califile->GetObject("T_Event",fNtuple2);
    if(!fNtuple2){ cout<<" Get Calibration file failed !!!"<<endl;return 0;}
    int Entries = fNtuple2->GetEntries();
    fNtuple2 ->SetBranchAddress("cycleID",&_cycleID);
    fNtuple2 ->SetBranchAddress("calicellID",&_cellID);
    fNtuple2 ->SetBranchAddress("calibcid",&_bcid);
    fNtuple2 ->SetBranchAddress("hitTags",&_hitTag);
    fNtuple2 ->SetBranchAddress("gainTags",&_gainTag);
    fNtuple2 ->SetBranchAddress("charges",&_charge);
    fNtuple2 ->SetBranchAddress("times",&_time);
    fNtuple2 ->SetBranchAddress("newtemp",&_temp);
    fNtuple2 ->SetBranchAddress("posX",&_posX);
    fNtuple2 ->SetBranchAddress("posY",&_posY);
    fNtuple2 ->SetBranchAddress("posZ",&_posZ);

    // *********************************************************************************
    ssafile = new TFile(ssa_name.c_str(),"RECREATE");
    tree = new TTree("T_Event","after ssa");
    //tree = fNtuple2->CloneTree(0);
    tree ->Branch("cycleID",&_cycleID_cali);
    tree ->Branch("cellIDs",&_cellID_cali);
    tree ->Branch("BCIDs",&_bcid_cali);
    tree ->Branch("hitTags",&_hitTag_cali);
    tree ->Branch("gainTags",&_gainTag_cali);
    tree ->Branch("charges",&_charge_cali);
    tree ->Branch("times",&_time_cali);
    tree ->Branch("temp",&_temp_cali);
    tree ->Branch("posX",&_posX_cali);
    tree ->Branch("posY",&_posY_cali);
    tree ->Branch("posZ",&_posZ_cali);
    tree ->Branch("ssaTag",&_ssaTag);



    
    // ==============================================================================
    //                         Strip Splitting Algorithm
    // ==============================================================================
    for(int entry=0; entry!=Entries; ++entry)
//    for(int entry=1; entry!=2; ++entry)
    {
        fNtuple2->GetEntry(entry);
        if(entry%10000==0)   
            cout<<" [ Process : StripSplit : ]  event : "<<entry<<" "<<endl;

        Reset();
        _cycleID_cali = _cycleID;
        vector<int> _hitLayer;
        for(int i_hit=0; i_hit!=(int)_bcid->size(); ++i_hit)
        {
            // only save hit=1 events, above threshold for auto-trigger
            if(_hitTag->at(i_hit)==0) continue;
            int layerID = _cellID->at(i_hit)/1e6;
            int chipIDs = _cellID->at(i_hit)%1000000/1e4;
            int channelIDs = _cellID->at(i_hit)%100;
            double posX = _posX->at(i_hit);
            double posY = _posY->at(i_hit);
            if(chipIDs==5&&channelIDs>29) continue;
            
/*            vector<double> highWeight(_splitNu);
            vector<double> lowWeight(_splitNu);
            vector<double> highWeightPos(_splitNu);
            vector<double> lowWeightPos(_splitNu);*/
            vector<double> highWeight;
            vector<double> lowWeight;
            vector<double> weightPosX;
            vector<double> weightPosY;
//            highWeight.clear();
//            lowWeight.clear();
//            weightPos.clear();
            double highWeightSum = 0.;
            double lowWeightSum = 0.;
            for(int j=0; j!=(int)_bcid->size(); ++j)
            {
                if (j==i_hit) continue;
                int layerIDis = _cellID->at(j)/1e6;
                double posXis = _posX->at(j);
                double posYis = _posY->at(j);
                double posYedge = 0.;
                double posXedge = 0.;
                double chargeis = _charge->at(j);
                double timeis = _time->at(j);
                if (layerIDis==layerID-1 || layerIDis==layerID+1) {
                     if (layerID%2==0) {
                        if (posYis>posY-0.5*_yInterval-0.5*_xInterval && posYis<posY+0.5*_yInterval+0.5*_xInterval) {
/*                            if (posYis<posY-0.5*_yInterval+0.5*_xInterval) {
                                chargeis = chargeis * ( (posYis+0.5*_yInterval) - (posY-0.5*_xInterval) )/_xInterval;
                                timeis = timeis * ( (posYis+0.5*_yInterval) - (posY-0.5*_xInterval) )/_xInterval;
                                posYedge = 0.5*( (posYis+0.5*_yInterval) - (posY-0.5*_xInterval) ) - 0.5*_xInterval;
                            } else if (posYis>posY+0.5*_yInterval-0.5*_xInterval) {
                                chargeis = chargeis * ( (posY+0.5*_xInterval) - (posYis-0.5*_yInterval) )/_xInterval;
                                timeis = timeis * ( (posY+0.5*_xInterval) - (posYis-0.5*_yInterval) )/_xInterval;
                                posYedge = -0.5*( (posY+0.5*_xInterval) - (posYis-0.5*_yInterval) ) + 0.5*_xInterval;
                            }
*/                            
                            if (posXis>posX-0.5*_yInterval-0.5*_xInterval && posXis<posX+0.5*_yInterval+0.5*_xInterval) {
/*                                if (posXis<posX-0.5*_yInterval+0.5*_xInterval) {
                                    chargeis = chargeis * ( (posXis+0.5*_xInterval) - (posX-0.5*_yInterval) )/_xInterval;
                                    timeis = timeis * ( (posXis+0.5*_xInterval) - (posX-0.5*_yInterval) )/_xInterval;
                                    posXedge = -0.5*( (posXis+0.5*_xInterval) - (posX-0.5*_yInterval) ) + 0.5*_xInterval;
                                } else if (posXis>posX+0.5*_yInterval-0.5*_xInterval) {
                                    chargeis = chargeis * ( (posX+0.5*_yInterval) - (posXis-0.5*_xInterval) )/_xInterval;
                                    timeis = timeis * ( (posX+0.5*_yInterval) - (posXis-0.5*_xInterval) )/_xInterval;
                                    posXedge = 0.5*( (posX+0.5*_yInterval) - (posXis-0.5*_xInterval) ) - 0.5*_xInterval;
                                }
*/
//                        if (posYis>posY-0.5*_yInterval && posYis<posY+0.5*_yInterval) {
//                            if (posXis>posX-0.5*_yInterval && posXis<posX+0.5*_yInterval) {

                                if (weightPosX.size()==0) {
                                    highWeight.push_back(timeis);
                                    lowWeight.push_back(chargeis);
                                    weightPosX.push_back(posXis + posXedge);
                                    weightPosY.push_back(posY + posYedge);
                                } else {
                                    bool weightPosBool = true;
                                    for (int k=0;k!=(int)weightPosX.size();k++) {
                                        if (weightPosX[k]>posXis+posXedge-0.5*_xInterval && weightPosX[k]<posXis+posXedge+0.5*_xInterval) {
                                            highWeight[k] += timeis;
                                            lowWeight[k] += chargeis;
                                            weightPosX[k] = (weightPosX[k] + posXis + posXedge) / 2.;
                                            weightPosBool = false;
                                        }
                                    }
                                    if (weightPosBool) {
                                        highWeight.push_back(timeis);
                                        lowWeight.push_back(chargeis);
                                        weightPosX.push_back(posXis + posXedge);
                                        weightPosY.push_back(posY + posYedge);
                                    }
                                }
                                highWeightSum += timeis;
                                lowWeightSum += chargeis;
                            }
                        }
                    }
                    if (layerID%2==1) {
                        if (posXis>posX-0.5*_yInterval-0.5*_xInterval && posXis<posX+0.5*_yInterval+0.5*_xInterval) {
/*                            if (posXis<posX-0.5*_yInterval+0.5*_xInterval) {
                                chargeis = chargeis * ( (posXis+0.5*_yInterval) - (posX-0.5*_xInterval) )/_xInterval;
                                timeis = timeis * ( (posXis+0.5*_yInterval) - (posX-0.5*_xInterval) )/_xInterval;
                                posXedge = 0.5*( (posXis+0.5*_yInterval) - (posX-0.5*_xInterval) ) - 0.5*_xInterval;
                            } else if (posXis>posX+0.5*_yInterval-0.5*_xInterval) {
                                chargeis = chargeis * ( (posX+0.5*_xInterval) - (posXis-0.5*_yInterval) )/_xInterval;
                                timeis = timeis * ( (posX+0.5*_xInterval) - (posXis-0.5*_yInterval) )/_xInterval;
                                posXedge = -0.5*( (posX+0.5*_xInterval) - (posXis-0.5*_yInterval) ) + 0.5*_xInterval;
                            }
*/
                            if (posYis>posY-0.5*_yInterval-0.5*_xInterval && posYis<posY+0.5*_yInterval+0.5*_xInterval) {
/*                                if (posYis<posY-0.5*_yInterval+0.5*_xInterval) {
                                    chargeis = chargeis * ( (posYis+0.5*_xInterval) - (posY-0.5*_yInterval) )/_xInterval;
                                    timeis = timeis * ( (posYis+0.5*_xInterval) - (posY-0.5*_yInterval) )/_xInterval;
                                    posYedge = -0.5*( (posYis+0.5*_xInterval) - (posY-0.5*_yInterval) ) + 0.5*_xInterval;
                                } else if (posYis>posY+0.5*_yInterval-0.5*_xInterval) {
                                    chargeis = chargeis * ( (posY+0.5*_yInterval) - (posYis-0.5*_xInterval) )/_xInterval;
                                    timeis = timeis * ( (posY+0.5*_yInterval) - (posYis-0.5*_xInterval) )/_xInterval;
                                    posYedge = 0.5*( (posY+0.5*_yInterval) - (posYis-0.5*_xInterval) ) - 0.5*_xInterval;
                                }
*/                                
//                        if (posXis>posX-0.5*_yInterval && posXis<posX+0.5*_yInterval) {
//                            if (posYis>posY-0.5*_yInterval && posYis<posY+0.5*_yInterval) {

                                if (weightPosY.size()==0) {
                                    highWeight.push_back(timeis);
                                    lowWeight.push_back(chargeis);
                                    weightPosX.push_back(posX + posXedge);
                                    weightPosY.push_back(posYis + posYedge);
                                } else {
                                    bool weightPosBool = true;
                                    for (int k=0;k!=(int)weightPosY.size();k++) {
                                        if (weightPosY[k]>posYis+posYedge-0.5*_xInterval && weightPosY[k]<posYis+posYedge+0.5*_xInterval) {
                                            highWeight[k] += timeis;
                                            lowWeight[k] += chargeis;
                                            weightPosY[k] = (weightPosY[k] + posYis + posYedge) / 2.;
                                            weightPosBool = false;
                                        }
                                    }
                                    if (weightPosBool) {
                                        highWeight.push_back(timeis);
                                        lowWeight.push_back(chargeis);
                                        weightPosX.push_back(posX + posXedge);
                                        weightPosY.push_back(posYis + posYedge);
                                    }
                                }
                                highWeightSum += timeis;
                                lowWeightSum += chargeis;
                            }
                        }
                    }
                }
            }
           
//            cout << "highWeightSum: " << highWeightSum << endl;
            if (highWeight.size()!=0 && highWeightSum>0.) {
                for(int k=0; k!=(int)highWeight.size(); ++k) {
                    if (highWeight[k]==0.) continue;
/*                    cout << "i_hit: " << i_hit << "  cellID: " << _cellID->at(i_hit) << endl;
                    cout << "highWeight[" << k << "]: " << highWeight[k] << endl;
                    cout << "weightPosX[" << k << "]: " << weightPosX[k] << endl;
                    cout << "weightPosY[" << k << "]: " << weightPosY[k] << endl;
                    cout << "posZ: " << _posZ->at(i_hit) << endl;*/
                    _cellID_cali.push_back(_cellID->at(i_hit));
                    _bcid_cali.push_back(_bcid->at(i_hit));
                    _hitTag_cali.push_back(_hitTag->at(i_hit));
                    _gainTag_cali.push_back(_gainTag->at(i_hit));
                    _temp_cali.push_back(_temp->at(i_hit));

                    double highWeightSplit = highWeight[k] / highWeightSum;
                    double lowWeightSplit = lowWeight[k] / lowWeightSum;
                    _charge_cali.push_back(_charge->at(i_hit) * lowWeightSplit);
                    _time_cali.push_back(_time->at(i_hit) * highWeightSplit);
                    _posX_cali.push_back(weightPosX[k]);
                    _posY_cali.push_back(weightPosY[k]);
                    _posZ_cali.push_back(_posZ->at(i_hit));
                    _ssaTag.push_back(1);
                }
            } else {
                _cellID_cali.push_back(_cellID->at(i_hit));
                _bcid_cali.push_back(_bcid->at(i_hit));
                _hitTag_cali.push_back(_hitTag->at(i_hit));
                _gainTag_cali.push_back(_gainTag->at(i_hit));
                _temp_cali.push_back(_temp->at(i_hit));
                _charge_cali.push_back(_charge->at(i_hit));
                _time_cali.push_back(_time->at(i_hit));
                _posX_cali.push_back(_posX->at(i_hit));
                _posY_cali.push_back(_posY->at(i_hit));
                _posZ_cali.push_back(_posZ->at(i_hit));
                _ssaTag.push_back(0);
            }
        }
        tree->Fill();
    }

    // *********************************************************************************
    ssafile->cd();
    tree->Write("",TObject::kOverwrite);
    ssafile->Close();
    califile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;

    return 0;
   // **********************************************************************************
}

