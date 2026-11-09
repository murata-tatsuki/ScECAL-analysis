#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include "TCanvas.h"
#include "langaus_newConvoluted.C"
#include "EBUdecode.h"
#include "EBUdecode.cxx"

using namespace std;


int _run_num_new, _event_time_new, _event_num, _detectorID;
std::vector<int> _cellID_new;
//std::vector<double> _hitTag_new;
std::vector<double> _hitE;
std::vector<double> _posx;
std::vector<double> _posy;
std::vector<double> _posz;
//std::vector<double> _hit_time_new;
//std::vector<double> _hg_cali;
//std::vector<double> _lg_cali;
std::vector<double> _temp;
double _tot_energy;
int _tot_hitlayer;


void Reset()
{	
	_cellID_new.clear();
	//_hit_time_newIclear();
    //_lg_cali.clear();
    //_hg_cali.clear();
	_hitE.clear();
	_temp.clear();
	_posx.clear();
	_posy.clear();
	_posz.clear();
}


int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./combine_all_data output.root input1.root input1.root  ..." << endl;
        return 1;
		}

    //cout << argv[1] << endl;
	cout << "=====>  " << argv[1] << endl;
  
	int rawfilenum = argc - 2;

  TFile *filein[rawfilenum];
	TTree *tree[rawfilenum];
	int entry_max[rawfilenum];
	int total_entry_max=0;
	
  for(int i=0; i<rawfilenum; i++){
		// filein[i] = new TFile(argv[filename[i].second]);
		filein[i] = new TFile(argv[i+2]);
		tree[i] = (TTree*) filein[i]->Get("Calib_Hit");
		entry_max[i] = tree[i]->GetEntries();
	}

	const int layerNu = 32;
	const int chipNu = 6;
	const int channelNu = 36;
	const int rowNu = 5;
	const int colNu = 42;
	const int doublelayerNu = 2;


	TFile fileout(argv[1],"RECREATE");

  TTree* outtree = new TTree("Calib_Hit","after calibration");	
	outtree ->Branch("Run_Num",&_run_num_new);
	outtree ->Branch("Event_Time",&_event_time_new);
	outtree ->Branch("Event_Num",&_event_num);
	outtree ->Branch("DetectorID",&_detectorID);
	outtree ->Branch("CellID",&_cellID_new);
  outtree ->Branch("Hit_Energy",&_hitE);
	outtree ->Branch("Hit_X",&_posx);
	outtree ->Branch("Hit_Y",&_posy);
	outtree ->Branch("Hit_Z",&_posz);
  //outtree ->Branch("Hit_Time",&_hit_time_new);
  //outtree ->Branch("HG_Cali",&_hg_cali);
  //outtree ->Branch("LG_Cali",&_lg_cali);
  outtree ->Branch("NewTemperature",&_temp);
  outtree ->Branch("TotalEnergyDep",&_tot_energy);
	//outtree ->Branch("TotalHitMap",&_tot_hitlayer);




	
	int _run_num_old=0, _event_time_old=0, _event_num_old=0, _detectorID_old=0;
  std::vector<int>* _cellID_old=nullptr;
  //std::vector<double> _hitTag_new;
  std::vector<double>* _hitE_old=nullptr;
  std::vector<double>* _posx_old=nullptr;
  std::vector<double>* _posy_old=nullptr;
  std::vector<double>* _posz_old=nullptr;
  //std::vector<double> _hit_time_new;
  //std::vector<double> _hg_cali;
  //std::vector<double> _lg_cali;
  std::vector<double>* _temp_old=nullptr;
  double _tot_energy_old=0;
  int _tot_hitlayer_old=0;






	// data をとってきてる
	for(int irawfile=0; irawfile<rawfilenum; irawfile++){
		if(rawfilenum>1) cout << irawfile << "/" << rawfilenum << endl;


		tree[irawfile] ->Branch("Run_Num",&_run_num_old);
		tree[irawfile] ->Branch("Event_Time",&_event_time_old);
		tree[irawfile] ->Branch("Event_Num",&_event_num_old);
		tree[irawfile] ->Branch("DetectorID",&_detectorID_old);
		tree[irawfile] ->Branch("CellID",&_cellID_old);
		tree[irawfile] ->Branch("Hit_Energy",&_hitE_old);
		tree[irawfile] ->Branch("Hit_X",&_posx_old);
		tree[irawfile] ->Branch("Hit_Y",&_posy_old);
		tree[irawfile] ->Branch("Hit_Z",&_posz_old);
    // tree[irawfile] ->Branch("Hit_Time",&_hit_time_new);
    // tree[irawfile] ->Branch("HG_Cali",&_hg_cali);
    // tree[irawfile] ->Branch("LG_Cali",&_lg_cali);
    tree[irawfile] ->Branch("NewTemperature",&_temp_old);
    tree[irawfile] ->Branch("TotalEnergyDep",&_tot_energy_old);
		//tree[irawfile] ->Branch("TotalHitMap",&_tot_hitlayer);

		for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
      // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
			tree[irawfile]->GetEntry(ientry);
			if(ientry==0) cout << _tot_energy_old << endl;
			Reset();


			_run_num_new = _run_num_old;
			_event_time_new = _event_time_old;
			_event_num = _event_num_old;
			_detectorID = _detectorID_old;
			_tot_energy = _tot_energy_old;

			for(unsigned ihit=0; ihit!=_cellID_old->size(); ++ihit) {
				_cellID_new.push_back(_cellID_old->at(ihit));
				_hitE.push_back(_hitE_old->at(ihit));
				_temp.push_back(_temp_old->at(ihit));
				_posx.push_back(_posx_old->at(ihit));
				_posy.push_back(_posy_old->at(ihit));
				_posz.push_back(_posz_old->at(ihit));
			}
			outtree->Fill();


		}
		filein[irawfile]->Close();

		_cellID_old = nullptr;
		_hitE_old = nullptr;
		_posx_old = nullptr;
		_posy_old = nullptr;
		_posz_old = nullptr;
		_temp_old = nullptr;
	}

	outtree->Write();


}