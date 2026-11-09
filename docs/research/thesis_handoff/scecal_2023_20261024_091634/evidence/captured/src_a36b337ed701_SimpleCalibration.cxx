#include<fstream>
#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TMath.h>
#include"SimpleCalibration.h"
#include"EBUdecode.h"

using namespace std;
//char char_tmp[1000];
void treeCalibration::PedestalExtract(string pedestal_name)
{
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
        cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
        for(unsigned i_hit=0; i_hit!=_pedCellID->size(); ++i_hit)
        {
            int layerID = _pedCellID->at(i_hit)/1e6;
            int chipIDs = _pedCellID->at(i_hit)%1000000/1e4;
          //  int memoryCells = _pedCellID->at(i_hit)%10000/1e2;
            int channelIDs = _pedCellID->at(i_hit)%100;
            pedMeanChargesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanCharges->at(i_hit);
            pedSigmaChargesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaCharges->at(i_hit);
            pedMeanTimesChnLevel[layerID][chipIDs][channelIDs] = _pedMeanTimes->at(i_hit);
            pedSigmaTimesChnLevel[layerID][chipIDs][channelIDs] = _pedSigmaTimes->at(i_hit);
        }
    }
}

void treeCalibration::Reset()
{
    _run_num_new=0;
    _event_time_new=0;
    _event_num=0;
    _detectorID=0;
    _cellID_new.clear();
    _hg_charge_new.clear();
    _lg_charge_new.clear();
    _temp.clear();
    _posX_cali.clear();
    _posY_cali.clear();
    _posZ_cali.clear();
}
int treeCalibration::Calibration(string root_name,string pedestal_name,string cali_name)
{
    string txtname(root_name);
    txtname.erase(0,txtname.find_last_of("/")+1);
    string caliname(cali_name);
    caliname.insert(caliname.length(), txtname);

    // pedestal_name.insert(pedestal_name.length(),txtname);
    cout<<pedestal_name<<endl;
    PedestalExtract(pedestal_name);
    
    // *********************************************************************************
    rawfile = TFile::Open(root_name.c_str(),"READ");
    cout<<root_name<<endl;
    rawfile->GetObject("Raw_Hit",fNtuple2);
    if(!fNtuple2){ cout<<" Get EventAlign file failed !!!"<<endl;return 0;}
    int Entries = fNtuple2->GetEntries();
    cout<<" [ Prepare to PedestalExtraction : ] Total Entries : "<<Entries<<" "<<endl;
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
    //tree = fNtuple2->CloneTree(0);
    tree ->Branch("Run_Num",&_run_num_new);
	tree ->Branch("Event_Time",&_event_time_new);
	tree ->Branch("Event_Num",&_event_num);
	tree ->Branch("DetectorID",&_detectorID);
    tree ->Branch("CellID",&_cellID_new);
    // tree ->Branch("BCIDs",&_bcid_cali);
    // tree ->Branch("hitTags",&_hitTag_cali);
    // tree ->Branch("gainTags",&_gainTag_cali);
    tree ->Branch("HG_Charge",&_hg_charge_new);
    tree ->Branch("LG_Charge",&_lg_charge_new);
    tree ->Branch("Hit_Energy",&_hg_charge_new);
    tree ->Branch("NewTemperature",&_temp);
    tree ->Branch("Hit_X",&_posX_cali);
    tree ->Branch("Hit_Y",&_posY_cali);
    tree ->Branch("Hit_Z",&_posZ_cali);
    
    int Entries_1=0; 
    int Entries_2=0;
    // ==============================================================================
    //                          Calibration signal 
    // ==============================================================================
    for(int entry=0; entry!=Entries; ++entry)
    {
        fNtuple2->GetEntry(entry);
        if(entry%10000==0)   
            cout<<" [ Process : Calibration : ]  event : "<<entry<<" "<<endl;

        Reset();
        // _cycleID_cali = _cycleID;
        _run_num_new = _run_num;
		_event_time_new = _event_time;
		_event_num = _triggerID;
		_detectorID = 0;
        vector<int> _hitLayer;
        int _totalHitStrips=0;
        for(unsigned i_hit=0; i_hit!=_cellID->size(); ++i_hit)
        {
            // only save hit=1 events, above threshold for auto-trigger
            if(_hitTag->at(i_hit)==0) continue;

            int layerID = _cellID->at(i_hit)/1e5;
            int chipIDs = (_cellID->at(i_hit)/10000)%10;
           // int memoryCells = _cellID->at(i_hit)%10000/1e2;
            int channelIDs = _cellID->at(i_hit)%100;
            // int charge = _lg_charge->at(i_hit); // for auto gain setting
            // int time = _hg_charge->at(i_hit); // for High&Low gain setting, High gain

            if(chipIDs==5&&channelIDs>29) continue;

            // subtract the pedestal for each channel 
            // double caliCharge = charge - pedMeanChargesChnLevel[layerID][chipIDs][channelIDs];
            // double caliTime   = time   - pedMeanTimesChnLevel[layerID][chipIDs][channelIDs];
            // if(caliTime<5*pedSigmaTimesChnLevel[layerID][chipIDs][channelIDs]) continue;
            //cout<<entry<<" : "<<time<<", "<<caliTime<<", "<<pedMeanTimesChnLevel[layerID][chipIDs][channelIDs] <<endl;

            vector<int>::iterator result = find(_hitLayer.begin(),_hitLayer.end(),layerID);
            if(result==_hitLayer.end()) _hitLayer.push_back(layerID);

            // decode the scintillator strips ID and positions corresponding
            double *_position = EBUdecode(layerID,chipIDs,channelIDs);
            double _xPosition = *_position;
            double _yPosition = *(_position+1);
            double _zPosition = *(_position+2);

            _totalHitStrips++;

            _cellID_new.push_back(_cellID->at(i_hit));
            // _bcid_cali.push_back(_bcid->at(i_hit));
            // _hitTag_cali.push_back(_hitTag->at(i_hit));
            // _gainTag_cali.push_back(_gainTag->at(i_hit));
            _hg_charge_new.push_back(_hg_charge->at(i_hit));
            _lg_charge_new.push_back(_lg_charge->at(i_hit));
            _posX_cali.push_back(_xPosition);
            _posY_cali.push_back(_yPosition);
            _posZ_cali.push_back(_zPosition);

            //SiPM gain temperature dependence correction
            double SiPM_temp;
            if(_tempLayer->at(layerID).size()!=0) SiPM_temp = tempReconstruction(layerID, _position, _tempLayer->at(layerID));
            else SiPM_temp=20;
            _temp.push_back(SiPM_temp);
        }
		// for test beam
        // if(_hitLayer.size()<22) continue;
        // if(_hitLayer.size()<15) continue;
        Entries_1++;
		// For test beam, this is not needed.
        // if(_totalHitStrips>64) continue;
        Entries_2++;
        tree->Fill();
    }

    cout<<"Total Entries : "<<Entries<<" , Layer>=22 : "<<Entries_1<<" , Strips<=64 : "<<Entries_2<<endl;
    // *********************************************************************************
    califile->cd();
    tree->Write("",TObject::kOverwrite);
    // califile->Close();
    // rawfile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;

   // **********************************************************************************
}

