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

void efficiency()
{

  ifstream filelist("./datalist_sps.txt");
  // ifstream filelist("./datalist_ps.txt");
  if(!filelist.good()){cout<<"open datalist file failed !"<<endl; exit(0);}


  TH1F *nentry = new TH1F("nentry","nentry",300,0,300);
  TGraph *detection_efficiency = new TGraph();
  detection_efficiency->GetXaxis()->SetTitle("layer");
  detection_efficiency->GetYaxis()->SetTitle("efficiency");
  double efficient[32] = {0}, total_entry = 0;

  for(int j=0; !filelist.eof(); j++){
    string dname;
    filelist>>dname;
    if(filelist.eof()) break;

    TFile* file = new TFile(dname.c_str());
    TTree* T_Event = (TTree*)file->Get("T_Event");
    if(!T_Event) cout<<"open RawtoRoot file failed "<<endl;

    cout<<j<<" : "<<dname<<" : "<<T_Event->GetEntries()<<endl;
    int cout_check=0;

    vector<int> *hitCellnew = nullptr;
    vector<int> *hitTags = nullptr;
    vector<double> *temp = nullptr;
    vector<double> *residualX = nullptr;
    vector<double> *residualY = nullptr;
    vector<double> *hitPosXnew = nullptr;
    vector<double> *hitPosYnew = nullptr;
    vector<double> *hitPosZnew = nullptr;
    // vector<vector<double>> *trackFitPars;
    vector<int> *layer_efficiency = nullptr;
    vector<double> *trackFitPars = nullptr;
    T_Event->SetBranchAddress("hitCellnew", &hitCellnew);
    T_Event->SetBranchAddress("Temperature", &temp);
    T_Event->SetBranchAddress("residualX", &residualX);
    T_Event->SetBranchAddress("residualY", &residualY);
    T_Event->SetBranchAddress("hitPosXnew", &hitPosXnew);
    T_Event->SetBranchAddress("hitPosYnew", &hitPosYnew);
    T_Event->SetBranchAddress("hitPosZnew", &hitPosZnew);
    T_Event->SetBranchAddress("trackFitPars", &trackFitPars);
    T_Event->SetBranchAddress("layer_efficiency", &layer_efficiency);
          double xNDFcut=9.6; double yNDFcut=9.6;
          double xNDFcutL=0.0; double yNDFcutL=0.0;

    for(int ientry=0;ientry<T_Event->GetEntries();ientry++){
      T_Event->GetEntry(ientry);
      nentry->Fill(hitCellnew->size());

      // if(hitCellnew->size()<22 || hitCellnew->size()>64) continue;
      vector<int> layers;
      for(int i=0;i<hitCellnew->size();i++){
        int _layer = hitCellnew->at(i)/1000000;
        auto result = find(layers.begin(), layers.end(), _layer);
        if (result == layers.end()) layers.push_back(_layer);
      }
      if(layers.size()<22) continue;
          double xSlope = trackFitPars->at(0);
          double xSlopeError = trackFitPars->at(1);
          double xIntercept = trackFitPars->at(2);
          double xInterceptError = trackFitPars->at(3);
          double ySlope = trackFitPars->at(4);
          double ySlopeError = trackFitPars->at(5);
          double yIntercept = trackFitPars->at(6);
          double yInterceptError = trackFitPars->at(7);
          double xLeft = xIntercept+xSlope*(12.95);
          double yLeft = yIntercept+ySlope*(12.95);
          double xRitgh = xIntercept+xSlope*(12.95+28/2*19.9);
          double yRitgh = yIntercept+ySlope*(12.95+28/2*19.9);
          // double ratio = (double)_cellID->size()/(double)_newCell->size();
          double xPhi = TMath::ATan(xSlope);
          double yPhi = TMath::ATan(ySlope);
          double theta = sqrt(pow(xSlope,2)+pow(ySlope,2));
          double angle = theta/TMath::Pi()*180;
          if(fabs(xRitgh)>114||fabs(xLeft)>114||fabs(yLeft)>114||fabs(yRitgh)>114||fabs(xPhi)>0.7||fabs(yPhi)>0.7) continue;

          vector<int> totalHitLayer;
          vector<int> totalHitLayerX;
          vector<int> totalHitLayerY;
          // int hitLayerNo[Layer]={0};
          double xNDF =0; double yNDF =0;
          vector<int> hitCell;
          for(unsigned hit=0; hit!=hitCellnew->size(); ++hit)
          {   
              hitCell.push_back(hitCellnew->at(hit));
              int _layerID = hitCellnew->at(hit)/1e5;
              double layerZ = hitPosZnew->at(hit);
              double _xFitting = xSlope*layerZ+xIntercept;
              double _yFitting = ySlope*layerZ+yIntercept;
              double _xResidual = hitPosXnew->at(hit)-_xFitting;
              double _yResidual = hitPosYnew->at(hit)-_yFitting;
              if(_layerID%2==0) yNDF += pow(_yResidual,2);
              if(_layerID%2==1) xNDF += pow(_xResidual,2);
              if(_layerID%2==1) {auto result = find(totalHitLayerX.begin(),totalHitLayerX.end(),_layerID);
                  if(result==totalHitLayerX.end()) totalHitLayerX.push_back(_layerID);}
              else  {auto result = find(totalHitLayerY.begin(),totalHitLayerY.end(),_layerID);
                  if(result==totalHitLayerY.end()) totalHitLayerY.push_back(_layerID);}
              auto result = find(totalHitLayer.begin(),totalHitLayer.end(),_layerID);
              if(result==totalHitLayer.end()) totalHitLayer.push_back(_layerID);
              // hitLayerNo[_layerID]++;
          }
          xNDF = (xNDF)/(hitCellnew->size()-2);
          yNDF = (yNDF)/(hitCellnew->size()-2);
          if(xNDF>=xNDFcut || yNDF>=yNDFcut) continue;
          if(xNDF<xNDFcutL || yNDF<yNDFcutL) continue;
          if(totalHitLayerX.size()<=6||totalHitLayerY.size()<=6) continue;
          // ******************************** event selection *********************************************


          if(hitCellnew->size()>40) continue;
          if(hitCellnew->size()>100) continue;
          total_entry++;
      for(int i=0; i<32;i++){
        if(layer_efficiency->at(i)==1) efficient[i]++;
      }
      // if(layer_hit.size()<22) continue;
      if(hitCellnew->size()>40) continue;
      if(hitCellnew->size()>100) continue;
    }

  }
  
  // TFile fileout("spsResult.root","RECREATE");


  for(int i=0; i<32;i++){
    double eff = total_entry!=0 ? efficient[i]/total_entry : 0;
    detection_efficiency->SetPoint(i,i,eff);
  }
  TCanvas *c = new TCanvas("c", "c", 2560, 1440);
  c->cd();
  detection_efficiency->SetMarkerStyle(4);
  detection_efficiency->Draw("AP");




}
