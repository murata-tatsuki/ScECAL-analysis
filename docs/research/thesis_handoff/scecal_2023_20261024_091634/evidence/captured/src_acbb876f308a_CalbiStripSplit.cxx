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
#include"CalibStripSplit.h"

using namespace std;

void treeCalibStripSplit::Reset()
{
    _run_num_new = 0;
    _event_time_new = 0;
    _detectorID_new = 0;
    _cellID_cali.clear();
    // _bcid_cali.clear();
    // _hitTag_cali.clear();
    // _gainTag_cali.clear();
    // _charge_cali.clear();
    // _time_cali.clear();
    // _temp_cali.clear();
    // _posX_cali.clear();
    // _posY_cali.clear();
    // _posZ_cali.clear();
    // _low_cali.clear();
    // _high_cali.clear();
    _hitE_cali.clear();
    _posX_cali.clear();
    _posY_cali.clear();
    _posZ_cali.clear();
    _temp_cali.clear();
    _ssaTag.clear();
}

int treeCalibStripSplit::StripSplit(string txt_name,string cali_name,string ssa_name)
{
    cali_name.insert(cali_name.length(), txt_name);
    ssa_name.insert(ssa_name.length(),txt_name);

    // *********************************************************************************
    califile = TFile::Open(cali_name.c_str(),"READ");
    califile->GetObject("Calib_Hit",fNtuple2);
    if(!fNtuple2){ cout<<" Get Calibration file failed !!!"<<endl;return 0;}
    int Entries = fNtuple2->GetEntries();
    fNtuple2 ->SetBranchAddress("Run_Num",&_run_num);
    fNtuple2 ->SetBranchAddress("Event_Time",&_event_time);
    fNtuple2 ->SetBranchAddress("Event_Num",&_event_num);
    fNtuple2 ->SetBranchAddress("DetectorID",&_detectorID);
    fNtuple2 ->SetBranchAddress("CellID",&_cellID);
    fNtuple2 ->SetBranchAddress("Hit_Energy",&_hitE);
    fNtuple2 ->SetBranchAddress("Hit_X",&_posX);
    fNtuple2 ->SetBranchAddress("Hit_Y",&_posY);
    fNtuple2 ->SetBranchAddress("Hit_Z",&_posZ);
    fNtuple2 ->SetBranchAddress("NewTemperature",&_temp);



    // *********************************************************************************
    ssafile = new TFile(ssa_name.c_str(),"RECREATE");
    tree = new TTree("SSA_Hit","after ssa");
    tree ->Branch("Run_Num",&_run_num_new);
    tree ->Branch("Event_Time",&_event_time_new);
    tree ->Branch("Event_Num",&_event_num_new);
    tree ->Branch("DetectorID",&_detectorID_new);
    tree ->Branch("CellID",&_cellID_cali);
    tree ->Branch("Hit_Energy",&_hitE_cali);
    tree ->Branch("Hit_X",&_posX_cali);
    tree ->Branch("Hit_Y",&_posY_cali);
    tree ->Branch("Hit_Z",&_posZ_cali);
    tree ->Branch("NewTemperature",&_temp_cali);
    tree ->Branch("ssaTag",&_ssaTag);
    
    // ==============================================================================
    //                         Strip Splitting Algorithm
    // ==============================================================================
    for(int entry=0; entry!=Entries; ++entry)
//    for(int entry=1; entry!=2; ++entry)
    {
        fNtuple2->GetEntry(entry);
        // if(entry%10000==0)   
        if(entry%1000==0)   
            cout<<" [ Process : StripSplit : ]  event : "<<entry<<" "<<endl;

        Reset();
        // _cycleID_cali = _cycleID;
        // _triggerID_cali = _triggerID;
        _run_num_new = _run_num;
        _event_num_new = _event_num;
        _event_time_new = _event_time;
        _detectorID_new = _detectorID;
        // 高速化1: レイヤーごとにヒットのインデックスを事前仕分け
        std::vector<int> hitsInLayer[40]; // 最大レイヤー数より少し大きめに確保
        for(int i_hit=0; i_hit!=(int)_cellID->size(); ++i_hit) {
            int lID = _cellID->at(i_hit)/1e5;
            if(lID >= 0 && lID < 40) {
                hitsInLayer[lID].push_back(i_hit);
            }
        }
        // 高速化2: 重複検索用のインデックスマップ
        std::map<int, std::vector<int>> cellToIndices;
        // 重複Hit防止およびエネルギー統合用ラムダ関数
        auto addHit = [&](int cell, double energy, double x, double y, double z, double temp, int tag) {
            // for (size_t idx = 0; idx < _cellID_cali.size(); ++idx) {
            for (int idx : cellToIndices[cell]) {
                // if (_cellID_cali[idx] == cell &&
                //     TMath::Abs(_posX_cali[idx] - x) < 1e-4 &&
                if (TMath::Abs(_posX_cali[idx] - x) < 1e-4 &&
                    TMath::Abs(_posY_cali[idx] - y) < 1e-4 &&
                    TMath::Abs(_posZ_cali[idx] - z) < 1e-4) {
                    _hitE_cali[idx] += energy; // 既に同じ位置が存在すればエネルギーを加算
                    return;
                }
            }
            int newIdx = _cellID_cali.size();
            _cellID_cali.push_back(cell);
            _hitE_cali.push_back(energy);
            _posX_cali.push_back(x);
            _posY_cali.push_back(y);
            _posZ_cali.push_back(z);
            _temp_cali.push_back(temp);
            _ssaTag.push_back(tag);
            cellToIndices[cell].push_back(newIdx);
        };
        for(int i_hit=0; i_hit!=(int)_cellID->size(); ++i_hit) {
            // only save hit=1 events, above threshold for auto-trigger
            // if(_hitTag->at(i_hit)==0) continue;

            int layerID = _cellID->at(i_hit)/1e5;
            int chipIDs = (_cellID->at(i_hit)/100000)%10;
            int channelIDs = _cellID->at(i_hit)%100;
            double posX = _posX->at(i_hit);
            double posY = _posY->at(i_hit);

            if(chipIDs==5&&channelIDs>29) continue;
            
/*            vector<double> highWeight(_splitNu);
            vector<double> lowWeight(_splitNu);
            vector<double> highWeightPos(_splitNu);
            vector<double> lowWeightPos(_splitNu);*/
            vector<double> energyWeight;
            vector<double> weightPosX;
            vector<double> weightPosY;
            vector<int> weightPosCount;
//            highWeight.clear();
//            lowWeight.clear();
//            weightPos.clear();
            double energyWeightSum = 0.;
            // for(int j=0; j!=(int)_cellID->size(); ++j){
            //     if (j==i_hit) continue;
            //     int layerIDis = _cellID->at(j)/1e5;
            //     double posXis = _posX->at(j);
            //     double posYis = _posY->at(j);
            //     double energyis = _hitE->at(j);
            //     if (layerIDis==layerID-1 || layerIDis==layerID+1) {
            for (int step = -1; step <= 1; step += 2) {
                int searchLayer = layerID + step;
                // 最大レイヤー数が40想定（必要に応じて上限値は変更してください）
                if (searchLayer < 0 || searchLayer >= 40) continue;
                
                for (int j : hitsInLayer[searchLayer]) {
                    int layerIDis = searchLayer;
                    double posXis = _posX->at(j);
                    double posYis = _posY->at(j);
                    double energyis = _hitE->at(j);
                    
                    // ストリップの半幅 (長さ 45mm, 幅 5mm)
                    double x_int_i = (layerID % 2 == 0) ? 45.0 / 2.0 : 5.0 / 2.0;
                    double y_int_i = (layerID % 2 == 0) ? 5.0 / 2.0 : 45.0 / 2.0;
                    double x_int_j = (layerIDis % 2 == 0) ? 45.0 / 2.0 : 5.0 / 2.0;
                    double y_int_j = (layerIDis % 2 == 0) ? 5.0 / 2.0 : 45.0 / 2.0;

                    // オーバーラップ面積の計算
                    double overX = TMath::Max(0., TMath::Min(posX + x_int_i, posXis + x_int_j) - TMath::Max(posX - x_int_i, posXis - x_int_j));
                    double overY = TMath::Max(0., TMath::Min(posY + y_int_i, posYis + y_int_j) - TMath::Max(posY - y_int_i, posYis - y_int_j));
                    double overlapArea = overX * overY;

                    if (overlapArea > 0.) {
                        double weight = energyis * overlapArea; // 面積とエネルギーで重み付け

                        double hitX = (layerID % 2 == 0) ? posXis : posX;
                        double hitY = (layerID % 2 == 0) ? posY : posYis;

                        // 最適クラスタの探索
                        int bestIndex = -1;
                        double minDistance = 2.5; // セル幅(5.0)の半分

                        for (int k = 0; k < (int)weightPosX.size(); k++) {
                            double currentX = weightPosX[k] / weightPosCount[k];
                            double currentY = weightPosY[k] / weightPosCount[k];
                            double dist = (layerID % 2 == 0) ? TMath::Abs(currentX - hitX) : TMath::Abs(currentY - hitY);
                            if (dist < minDistance) {
                                minDistance = dist;
                                bestIndex = k;
                            }
                        }

                        if (bestIndex != -1) {
                            energyWeight[bestIndex] += weight;
                            weightPosX[bestIndex] += hitX;
                            weightPosY[bestIndex] += hitY;
                            weightPosCount[bestIndex]++;
                        } else {
                            energyWeight.push_back(weight);
                            weightPosX.push_back(hitX);
                            weightPosY.push_back(hitY);
                            weightPosCount.push_back(1);
                        }
                        energyWeightSum += weight;
                    }
                }
            }
            if (energyWeight.size()!=0 && energyWeightSum>0.) {
                for(int k=0; k!=(int)energyWeight.size(); ++k) {
                    if (energyWeight[k]==0.) continue;
                    double energyWeightSplit = energyWeight[k] / energyWeightSum;
                    double avgX = weightPosX[k] / weightPosCount[k];
                    double avgY = weightPosY[k] / weightPosCount[k];
                    addHit(_cellID->at(i_hit), _hitE->at(i_hit) * energyWeightSplit, avgX, avgY, _posZ->at(i_hit), _temp->at(i_hit), 1);
                }
            } else {
                addHit(_cellID->at(i_hit), _hitE->at(i_hit), _posX->at(i_hit), _posY->at(i_hit), _posZ->at(i_hit), _temp->at(i_hit), 0);
            }
        }
        tree->Fill();
    }

    // *********************************************************************************
    ssafile->cd();
    tree->Write("",TObject::kOverwrite);
    ssafile->Close();
    delete ssafile;
    califile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;

    return 0;
   // **********************************************************************************
}

