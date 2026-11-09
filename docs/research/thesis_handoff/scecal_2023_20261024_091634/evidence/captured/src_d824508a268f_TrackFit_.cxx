#include<fstream>
#include<iostream>
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include"TrackFit_.h"
#include"../include/EBUdecode.h"

using namespace std;

static vector<float> xPosition;
static vector<float> yPosition;
static vector<float> zxPosition;
static vector<float> zyPosition;
//----------------------------------------------------
void trackFit_::Clear()
{
    _trackFitPars.clear();
    _residualX.clear();
    _residualY.clear();
    _chargeCor.clear();
    _timeCor.clear();

    _hitCell.clear();
    _hitPosX.clear();
    _hitPosY.clear();
    _hitPosZ.clear();
    _hitCharge.clear();

    _hitCellnew.clear();
    _hitPosXnew.clear();
    _hitPosYnew.clear();
    _hitPosZnew.clear();
    _layer_efficiency.clear();

    xPosition.clear();
    yPosition.clear();
    zxPosition.clear();
    zyPosition.clear();
}
//----------------------------------------------------
void trackFit_::Reset()
{
    xPosition.clear();
    yPosition.clear();
    zxPosition.clear();
    zyPosition.clear();
}
//----------------------------------------------------
// definition of function of Track Fitting
//----------------------------------------------------
void fcnX_(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = xPosition.size();
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (int i=0;i<nbins; i++) {
        delta  = (xPosition[i]-(zxPosition[i]*par[0]+par[1]));
        chisq += delta*delta;
    }
    f = chisq;
}
void fcnY_(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    const Int_t nbins = yPosition.size();
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta;
    for (int i=0;i<nbins; i++) {
        delta  = (yPosition[i]-(zyPosition[i]*par[0]+par[1]));
        chisq += delta*delta;
    }
    f = chisq;
}
//----------------------------------------------------
double* trackFit_::TrackFit(int slience)
{
    Double_t arglist[10];
    Int_t ierflg = 0;
    static Double_t vstart[2];
    static Double_t step[2];
    static double fitpar[8]={0};
    if(slience) cout<<"===================================================================="<<endl;
    //initialize TMinuit with a maximum of 2 params for XZ
    gMinuitX = new TMinuit();  
    gMinuitX->SetFCN(fcnX_);
    gMinuitX->SetPrintLevel(-1);

    arglist[0] = 1;
    gMinuitX->mnexcm("SET ERR", arglist ,1,ierflg);

    Int_t nbins = xPosition.size();
    //for(int i=0; i!=nbins; ++i){ cout<<xPosition[i]<<" , "<<endl;}
    double zLeft=291.55; double zRight=1; 
    for(int i=0; i!=nbins; ++i){if(zxPosition[i]<zLeft) zLeft=zxPosition[i];if(zxPosition[i]>zRight) zRight=zxPosition[i];}
    int itXleft, itXright;
    for(int i=0; i!=nbins; ++i){if(zxPosition[i]==zLeft) itXleft=i;if(zxPosition[i]==zRight) itXright=i;}
    double xLeft = xPosition[itXleft];
    double xRight = xPosition[itXright];
    if((zRight-zLeft)==0) vstart[0] = 0;
    else vstart[0] = (xRight-xLeft)/(zRight-zLeft);
    //vstart[0] = (xPosition[nbins-1]-xPosition[0])/(zxPosition[nbins-1]-zxPosition[0]);
    vstart[1] = xPosition[0] - vstart[0]*zxPosition[0];
    step[0]=fabs(0.05*vstart[0]);
    step[1]=fabs(0.05*vstart[1]);
    //cout<<"Initial Parameter x : "<<vstart[0]<<" "<<vstart[1]<<" "<<step[0]<<" "<<step[1]<<endl;
    gMinuitX->mnparm(0, "Xslope", vstart[0], step[0], -11000,11000,ierflg);
    gMinuitX->mnparm(1, "Xintercept", vstart[1], step[1], -10000,10000,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuitX->mnexcm("MIGRAD", arglist ,0,ierflg);

    gMinuitX->GetParameter(0,fitpar[0],fitpar[1]);
    gMinuitX->GetParameter(1,fitpar[2],fitpar[3]);

   if(slience)  cout<<"===================================================================="<<endl;
    //initialize TMinuit with a maximum of 2 params for YZ
    gMinuitY = new TMinuit();  
    gMinuitY->SetFCN(fcnY_);
    gMinuitY->SetPrintLevel(-1);

    arglist[0] = 1;
    gMinuitY->mnexcm("SET ERR", arglist ,1,ierflg);

    nbins = yPosition.size();
    //for(int i=0; i!=nbins; ++i){ cout<<yPosition[i]<<" , "<<endl;}
    zLeft=291.55; zRight=1;
    for(int i=0; i!=nbins; ++i){if(zyPosition[i]<zLeft) zLeft=zyPosition[i];if(zyPosition[i]>zRight) zRight=zyPosition[i];}
    int itYleft, itYright;
    for(int i=0; i!=nbins; ++i){if(zyPosition[i]==zLeft) itYleft=i;if(zyPosition[i]==zRight) itYright=i;}
    double yLeft = yPosition[itYleft];
    double yRight = yPosition[itYright];
    if((zRight-zLeft)==0) vstart[0] = 0;
    else vstart[0] = (yRight-yLeft)/(zRight-zLeft);
    //vstart[0] = (yPosition[nbins-1]-yPosition[0])/(zyPosition[nbins-1]-zyPosition[0]);
    vstart[1] = yPosition[0] - vstart[0]*zyPosition[0];
    step[0]=fabs(0.05*vstart[0]);
    step[1]=fabs(0.05*vstart[1]);
    //cout<<"Initial Parameter y : "<<vstart[0]<<" "<<vstart[1]<<" "<<step[0]<<" "<<step[1]<<endl;
    gMinuitY->mnparm(0, "Yslope", vstart[0], step[0], -11000,11000,ierflg);
    gMinuitY->mnparm(1, "Yintercept", vstart[1], step[1], -10000,10000,ierflg);

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuitY->mnexcm("MIGRAD", arglist ,0,ierflg);

    gMinuitY->GetParameter(0,fitpar[4],fitpar[5]);
    gMinuitY->GetParameter(1,fitpar[6],fitpar[7]);

    if(slience){ cout<<"Fitting Pars XZ slope :"<< fitpar[0]<<" +- "<< fitpar[1]<<" intercept : "<< fitpar[2]<< " +- " << fitpar[3]<<endl;
    cout<<"Fitting Pars YZ slope :"<< fitpar[4]<<" +- "<< fitpar[5]<<" intercept : "<< fitpar[6]<< " +- " << fitpar[7]<<endl;
    cout<<"===================================================================="<<endl;}

    delete gMinuitX; gMinuitX=NULL;
    delete gMinuitY; gMinuitY=NULL;
    return fitpar;
}   
//----------------------------------------------------
trackFit_::trackFit_(string root_name,string track_name, string layer)
{
    cout<<" ----------- "<<root_name<<" ------------"<<endl;
    string txtname(root_name);
    string rootname(root_name);
    rootname.erase(0,rootname.find_last_of("/")+1);
    cout<<" ----------- "<<rootname<<" ------------"<<endl;
    //txtname.insert( txtname.length(), ".root");

    // string txtname(root_name);
    rootfile = new TFile(txtname.c_str());
    if(!rootfile) {cout<<" open align file falied"<<endl; return; }
    fNtuple1 = (TTree*)rootfile->Get("Calib_Hit");
    int Entries = fNtuple1->GetEntries();
    cout<<" ------------ Total Events = "<<Entries<<" ---------------"<<endl;

    fNtuple1 ->SetBranchAddress("Event_Num",&_EventNum);
    fNtuple1 ->SetBranchAddress("DetectorID",&_DetectorID);
    fNtuple1 ->SetBranchAddress("CellID",&_CellID);
    fNtuple1 ->SetBranchAddress("Hit_Energy",&_Hit_Energy);
    // fNtuple1 ->SetBranchAddress("Hit_Time",&_Hit_Time);
    fNtuple1 ->SetBranchAddress("Hit_X",&_Hit_X);
    fNtuple1 ->SetBranchAddress("Hit_Y",&_Hit_Y);
    fNtuple1 ->SetBranchAddress("Hit_Z",&_Hit_Z);

    string trackname(track_name);
    trackname.insert( trackname.length(), rootname);
    fitfile = new TFile(trackname.c_str(),"RECREATE");
    fTree = new TTree("T_Event","");
    // fTree = fNtuple1->CloneTree(0);
    fTree ->Branch("TriggerID",&_EventNum_new);
    fTree ->Branch("trackFitPars",&_trackFitPars);
    fTree ->Branch("residualX",&_residualX);
    fTree ->Branch("residualY",&_residualY);
    fTree ->Branch("energyDep",&_timeCor);
    fTree ->Branch("realLength",&_realLength);
    //fTree ->Branch("cellIDs",&_hitCellnew);
    fTree ->Branch("hitCellnew",&_hitCellnew);
    fTree ->Branch("hitPosXnew",&_hitPosXnew);
    fTree ->Branch("hitPosYnew",&_hitPosYnew);
    fTree ->Branch("hitPosZnew",&_hitPosZnew);
    fTree ->Branch("layer_efficiency", &_layer_efficiency);

    // ========================= Total hit events for each channel ===================

    int _totalEntries1=0; int _totalEntries2=0;
    int _totalEntries3=0; int _totalEntries4=0;
    int _selectedLayerID =0; // stoi(event);
    stringstream degree(layer);
    degree >> _selectedLayerID;
    for(int entry=0; entry!=fNtuple1->GetEntries();entry++)
    {
        fNtuple1->GetEntry(entry);
        Clear();
        if(entry%1000==0) cout<<" Process : Events = "<<entry<</*" : "<<_bcid->size()<<*/endl;
        _EventNum_new = _EventNum;
        for(unsigned hit=0; hit!=_CellID->size(); ++hit)
        {
            if(_DetectorID!=0) continue;
            int _layerID = _CellID->at(hit)/1e5;
            int _chipID = (_CellID->at(hit)/10000)%10;
            //  int _cell = _cellID->at(hit)%10000/1e2;
            int _chnID = _CellID->at(hit)%100;
            //  int gainTag = _gainTag->at(hit);
            double charge = _Hit_Energy->at(hit);
            //int time = _time->at(hit);
            if(_layerID>=16) continue;

            int scintillatorID = _chipID*chnNu+_chnID;

            _hitCell.push_back(_CellID->at(hit));
            _position = EBUdecode(_layerID,_chipID,_chnID);
            // alignment for each EBU based on the track fitting
            if(_layerID%2==1) *_position = *(_position) - EBUposX[_layerID/2] - EBUposX2[_layerID/2];
            if(_layerID%2==0) *(_position+1) = *(_position+1) - EBUposY[_layerID/2] - EBUposY2[_layerID/2];
            // _hitPosX.push_back(*_position);
            // _hitPosY.push_back(*(_position+1));
            // _hitPosZ.push_back(*(_position+2));
            _hitPosX.push_back(_Hit_X->at(hit));
            _hitPosY.push_back(_Hit_Y->at(hit));
            _hitPosZ.push_back(_Hit_Z->at(hit));
            _hitCharge.push_back((double)charge);
            //cout<<_layerID<<" : "<<hit<<" : "<<*_position<<" : "<<*(_position+1)<<" : "<<*(_position+2)<<endl;
            // if(entry<3) cout << "(" << _Hit_X->at(hit) << "," << _Hit_Y->at(hit) << "," << _Hit_Z->at(hit) << ")    " << "(" << *_position << "," << *(_position+1) << "," << *(_position+2) << ")  " << endl;
        }
        
        // ---------------------- iteration fitting 1 ---------------------
        double *_fitPar;
        for(unsigned _chanl=0; _chanl!=_hitCell.size(); ++_chanl)
        {
            int layerID = _hitCell.at(_chanl)/1e5;
            double _xPosition = _hitPosX.at(_chanl);
            double _yPosition = _hitPosY.at(_chanl);
            double _zPosition = _hitPosZ.at(_chanl);
            double _signal = _hitCharge.at(_chanl);
            if(layerID%2==0) {
                yPosition.push_back(_yPosition);
                zyPosition.push_back(_zPosition);
            }
            if(layerID%2==1) {
                xPosition.push_back(_xPosition);
                zxPosition.push_back(_zPosition);
            }
        }
        if(xPosition.size()<3||yPosition.size()<3) continue;
        _fitPar = TrackFit(0);
        _totalEntries1++;
        // ---------------------- iteration fitting 2 ---------------------
        const int itNo=9;
        double step = 5;
        double start = 47.5;
        for(int it=0; it!=itNo; it++){
            Reset();
            for(unsigned _chanl=0; _chanl!=_hitCell.size(); ++_chanl)
            {
                int layerID = _hitCell.at(_chanl)/1e5;
                double _xPosition = _hitPosX.at(_chanl);
                double _yPosition = _hitPosY.at(_chanl);
                double _zPosition = _hitPosZ.at(_chanl);
                double xSlope = *(_fitPar);
                double xIntercept = *(_fitPar+2);
                double ySlope = *(_fitPar+4);
                double yIntercept = *(_fitPar+6);
                double _xFitting = xSlope*_zPosition+xIntercept;
                double _yFitting = ySlope*_zPosition+yIntercept;
                if(layerID%2==0) {
                    if(fabs(_yFitting-_yPosition)>start-step*it) continue;
                    yPosition.push_back(_yPosition);
                    zyPosition.push_back(_zPosition);
                }
                if(layerID%2==1) {
                    if(fabs(_xFitting-_xPosition)>start-step*it) continue;
                    xPosition.push_back(_xPosition);
                    zxPosition.push_back(_zPosition);
                }
            }
            if(xPosition.size()<3||yPosition.size()<3) continue;
            _fitPar = TrackFit(0);
        }

        // ====================== selection of track fitting 1  ================================
        double xSlope = *(_fitPar);
        double xIntercept = *(_fitPar+2);
        double ySlope = *(_fitPar+4);
        double yIntercept = *(_fitPar+6);
        double xLeft = xIntercept+xSlope*(12.95);
        double yLeft = yIntercept+ySlope*(12.95);
        double xRitgh = xIntercept+xSlope*(12.95+28/2*19.9);
        double yRitgh = yIntercept+ySlope*(12.95+28/2*19.9);
        double xPhi = TMath::ATan(xSlope);
        double yPhi = TMath::ATan(ySlope);
        if(fabs(xRitgh)>114||fabs(xLeft)>114||fabs(yLeft)>114||fabs(yRitgh)>114||fabs(xPhi)>0.7||fabs(yPhi)>0.7) continue;
        _totalEntries2++;

       // cout<<" ============================= finish track fitting =================================="<<endl;
        vector< vector<double> > minHit; minHit.resize(layerNu);
        vector< vector<double> > minCell; minCell.resize(layerNu);
        vector< vector<double> > hitCellX; hitCellX.resize(layerNu);
        vector< vector<double> > hitCellY; hitCellY.resize(layerNu);
        for(unsigned hit=0; hit!=_hitCell.size(); ++hit)
        {
            int layerID = _hitCell.at(hit)/1e5;
            double _xPosition = _hitPosX.at(hit);
            double _yPosition = _hitPosY.at(hit);
            double _zPosition = _hitPosZ.at(hit);
            double _xFitting = xSlope*_zPosition+xIntercept;
            double _yFitting = ySlope*_zPosition+yIntercept;
            double _xResidual = fabs(_xPosition-_xFitting);
            double _yResidual = fabs(_yPosition-_yFitting);

            if(layerID%2==0){
                if(minHit[layerID].size()==0) {minHit[layerID].push_back(_yResidual); minCell[layerID].push_back(_hitCell.at(hit));}
                else if(_yResidual<minHit[layerID].at(0)) {minHit[layerID].at(0) = _yResidual; minCell[layerID].at(0) = _hitCell.at(hit);}
            }
            if(layerID%2==1){
                if(minHit[layerID].size()==0) {minHit[layerID].push_back(_xResidual); minCell[layerID].push_back(_hitCell.at(hit));}
                else if(_yResidual<minHit[layerID].at(0)) {minHit[layerID].at(0) = _xResidual; minCell[layerID].at(0) = _hitCell.at(hit);}
            }
        }
        Reset();

        int _realHitX=0; int _realHitY=0;
        for(unsigned _chanl=0; _chanl!=_hitCell.size(); ++_chanl)
        {
            int layerID = _hitCell.at(_chanl)/1e5;
            double _xPosition = _hitPosX.at(_chanl);
            double _yPosition = _hitPosY.at(_chanl);
            double _zPosition = _hitPosZ.at(_chanl);
            double _xFitting = xSlope*_zPosition+xIntercept;
            double _yFitting = ySlope*_zPosition+yIntercept;
            double _xResidual = fabs(_xPosition-_xFitting);
            double _yResidual = fabs(_yPosition-_yFitting);
            int trueHit=0;
            //if(layerID%2==0&&fabs(_yFitting-_yPosition)<=7.5) {trueHit++;_realHitY++;}
            //if(layerID%2==1&&fabs(_xFitting-_xPosition)<=7.5) {trueHit++;_realHitX++;}
            if(layerID%2==0&&_yResidual<=7.5&&minCell[layerID].at(0)==_hitCell.at(_chanl)) {trueHit++;_realHitY++;}
            if(layerID%2==1&&_xResidual<=7.5&&minCell[layerID].at(0)==_hitCell.at(_chanl)) {trueHit++;_realHitX++;}
            if(trueHit==0) continue;
            //_hitCellnew.push_back(layerID);
            _hitCellnew.push_back(_hitCell.at(_chanl));
            _hitPosXnew.push_back(_xPosition);
            _hitPosYnew.push_back(_yPosition);
            _hitPosZnew.push_back(_zPosition);
            if(layerID%2==0) {
                yPosition.push_back(_yPosition);
                zyPosition.push_back(_zPosition);
                hitCellY[layerID].push_back(_hitCell.at(_chanl));
            }
            if(layerID%2==1) {
                xPosition.push_back(_xPosition);
                zxPosition.push_back(_zPosition);
                hitCellX[layerID].push_back(_hitCell.at(_chanl));
            }
        }
        if(_realHitX<3||_realHitY<3) continue;
        _fitPar = TrackFit(0);
        for(int ii=0; ii!=8; ++ii){_trackFitPars.push_back(*(_fitPar+ii));}

        // ---------------------- coincidence of trigger layers ---------------------
        double _hitLayer[layerNu]={0};
        for(unsigned _chanl=0; _chanl!=_hitCellnew.size(); ++_chanl)
        {
            int layerID = _hitCellnew.at(_chanl)/1e5;
            _hitLayer[layerID]=1;
        }
        // if(_hitLayer[1]==0||_hitLayer[29]==0) continue;

        // ====================== selection of track fitting 2 ================================
        double xNDF =0; double yNDF =0;
        vector< int > totalHitLayerX;
        vector< int > totalHitLayerY;
        for(unsigned _chanl=0; _chanl!=_hitCellnew.size(); ++_chanl)
        {
            int layerID = _hitCellnew.at(_chanl)/1e5;
            double _xPosition = _hitPosXnew.at(_chanl);
            double _yPosition = _hitPosYnew.at(_chanl);
            double _zPosition = _hitPosZnew.at(_chanl);
            double _xFitting = xSlope*_zPosition+xIntercept;
            double _yFitting = ySlope*_zPosition+yIntercept;
            double _xResidual = _xPosition-_xFitting;
            double _yResidual = _yPosition-_yFitting;
            if(layerID%2==0) {yNDF += pow(_yResidual,2);}
            if(layerID%2==1) {xNDF += pow(_xResidual,2);}
            if(layerID%2==1) {vector<int>::iterator result = find(totalHitLayerX.begin(),totalHitLayerX.end(),layerID);
                if(result==totalHitLayerX.end()) totalHitLayerX.push_back(layerID);}
            else  {vector<int>::iterator result = find(totalHitLayerY.begin(),totalHitLayerY.end(),layerID);
                if(result==totalHitLayerY.end()) totalHitLayerY.push_back(layerID);}
        }
        xNDF = (xNDF)/(_hitCellnew.size()-2);
        yNDF = (yNDF)/(_hitCellnew.size()-2);
        //if(xNDF>=9.6 || yNDF>=9.6) continue;
        if(totalHitLayerX.size()<=2||totalHitLayerY.size()<=2) continue;
        // if(totalHitLayerX.size()<=6||totalHitLayerY.size()<=6) continue;
        _totalEntries3++;

        // ====================== ADC correction by track fitting ==============================
        const double thickness = 2; // thickness of scintillator strip 
        //double xSlope = *(_fitPar);
        //double ySlope = *(_fitPar+4);
        double length = thickness*sqrt(pow(xSlope,2)+pow(ySlope,2)+pow(1,2));
        for(unsigned _chanl=0; _chanl!=_hitCellnew.size(); ++_chanl)
        {
            double _signal = _hitCharge.at(_chanl);
            double _signalCor = _signal*thickness/length;
            _timeCor.push_back(_signalCor);
        }
        _realLength = length;
        if(length>5) cout<<" 111111111111 "<<entry<<" : "<<length<<endl;

        // ====================== position residual calculation ================================
        vector< vector<double> > _trackPar; _trackPar.resize(layerNu);
        double *_trackFit;
        for(int s_layer=0; s_layer!=layerNu; ++s_layer)
        {
            // cout<<" --------- _selectedLayerID : "<<s_layer<<endl;
            xPosition.clear();
            yPosition.clear();
            zxPosition.clear();
            zyPosition.clear();
            for(unsigned _chanl=0; _chanl!=_hitCellnew.size(); ++_chanl)
            {
                int layerID = _hitCellnew.at(_chanl)/1e5;
                if(layerID==s_layer) continue;
                double _xPosition = _hitPosXnew.at(_chanl);
                double _yPosition = _hitPosYnew.at(_chanl);
                double _zPosition = _hitPosZnew.at(_chanl);
                if(layerID%2==0) {
                    yPosition.push_back(_yPosition);
                    zyPosition.push_back(_zPosition);
                }
                if(layerID%2==1) {
                    xPosition.push_back(_xPosition);
                    zxPosition.push_back(_zPosition);
                }
            }
            if(xPosition.size()<3||yPosition.size()<3) {
                cout<<xPosition.size()<<" : "<<yPosition.size()<<" The number of real hits too few to fitting"<<endl;
                //continue;
            }
            _trackFit = TrackFit(0);
            for(int ii=0; ii!=8; ++ii) _trackPar[s_layer].push_back(*(_trackFit+ii));
            //cout<<" ====================== position residual calculation ================================"<<endl;
        }

        for(unsigned _chanl=0; _chanl!=_hitCellnew.size(); ++_chanl)
        {
            int i_layer = _hitCellnew.at(_chanl)/1e5;
            double _xPosition = _hitPosXnew.at(_chanl);
            double _yPosition = _hitPosYnew.at(_chanl);
            double _zPosition = _hitPosZnew.at(_chanl);
            double xSlope = _trackPar[i_layer].at(0);
            double xIntercept = _trackPar[i_layer].at(2);
            double ySlope = _trackPar[i_layer].at(4);
            double yIntercept = _trackPar[i_layer].at(6);
            double _xFitting = xSlope*_zPosition+xIntercept;
            double _yFitting = ySlope*_zPosition+yIntercept;
            double _xResidual = _xPosition-_xFitting;
            double _yResidual = _yPosition-_yFitting;
            _residualX.push_back(_xResidual);
            _residualY.push_back(_yResidual);
        }

        for (int s_layer = 0; s_layer != layerNu; ++s_layer) {
            int result = 0;

            double xSlope = _trackPar[s_layer].at(0);
            double xIntercept = _trackPar[s_layer].at(2);
            double ySlope = _trackPar[s_layer].at(4);
            double yIntercept = _trackPar[s_layer].at(6);

            double _zPosition = 0, x_size=0, y_size=0;
            if(s_layer%2==0) {
                _zPosition = 1+s_layer/2*19.9;
                x_size = 22.5;
                y_size = 2.5;
            } else {
                _zPosition = 12.2+(s_layer-1)/2*19.9;
                x_size = 2.5;
                y_size = 22.5;
            }

            double _xFitting = xSlope * _zPosition + xIntercept;
            double _yFitting = ySlope * _zPosition + yIntercept;

bool cout_CCheck = false;
if(entry<3 &&cout_CCheck){
  cout << "   layer: " << s_layer << endl;
  cout << "      xSlope:" << xSlope << ", xIntercept:" << xIntercept << ", ySlope:" << ySlope << ", yIntercept:" << yIntercept << endl;
}

            for (unsigned _chanl = 0; _chanl != _hitCellnew.size(); ++_chanl) {
                int i_layer = _hitCellnew.at(_chanl) / 1e5;
                if(i_layer != s_layer) continue;

                double _xPosition = _hitPosXnew.at(_chanl);
                double _yPosition = _hitPosYnew.at(_chanl);


if(entry<3 &&cout_CCheck){
  cout << "         x,y= " << _xPosition << "," << _yPosition << "      fitting " << _xFitting << "," << _yFitting << endl;
}



                if(s_layer%2==1 && _xPosition-x_size<=_xFitting && _xFitting<=_xPosition+x_size) result = 1;
                if(s_layer%2==0 && _yPosition-y_size<=_yFitting && _yFitting<=_yPosition+y_size) result = 1;
                //if(_xPosition-x_size<=_xFitting && _xFitting<=_xPosition+x_size && _yPosition-y_size<=_yFitting && _yFitting<=_yPosition+y_size) result = 1;
            }
            _layer_efficiency.push_back(result);
        }

        _totalEntries4++;
        fTree->Fill();
    }
    cout<<_totalEntries1<<" : "<<_totalEntries2<<" : "<<_totalEntries3<<" : "<<_totalEntries4<<endl;

    // *********************************************************************************
    // *********************************************************************************

    fitfile->cd();
    fTree->Write("",TObject::kOverwrite);
    fitfile->Close();
}
