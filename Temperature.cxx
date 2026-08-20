#include<fstream>
#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include"Temperature.h"
#include"EBUdecode.h"
#include<TMath.h>

using namespace std;

void treeTemperature::Init() 
{
    char hChargeName[200];
    Hlist1 = new TObjArray(10);
    Hlist2 = new TObjArray(10);
    Hlist3 = new TObjArray(10);
    Hlist4 = new TObjArray(10);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                for(int k=0; k!=memoryNu; ++k)
                {
                   // sprintf(hChargeName,"Charges_layer%d_chip%d_chan%d_cell%d",l,i,j,k);
                   // hCharge1[l][i][j][k] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                   // Hlist1->Add(hCharge1[l][i][j][k]);
                   // sprintf(hChargeName,"Times_layer%d_chip%d_chan%d_cell%d",l,i,j,k);
                   // hTime1[l][i][j][k] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                   // Hlist1->Add(hTime1[l][i][j][k]);
                }
                sprintf(hChargeName,"Charges_layer%d_chip%d_chan%d",l,i,j);
                hCharge2[l][i][j] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                Hlist2->Add(hCharge2[l][i][j]);                
                sprintf(hChargeName,"Times_layer%d_chip%d_chan%d",l,i,j);
                hTime2[l][i][j] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                Hlist2->Add(hTime2[l][i][j]);

                // temperature
                sprintf(hChargeName,"Charges_layer%d_chip%d_chan%d_temp",l,i,j);
                hCharge4[l][i][j] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                Hlist4->Add(hCharge4[l][i][j]);                
                sprintf(hChargeName,"Times_layer%d_chip%d_chan%d_temp",l,i,j);
                hTime4[l][i][j] = new TH1F(hChargeName,hChargeName,4096,0,4095);
                Hlist4->Add(hTime4[l][i][j]);
                // graph
                sprintf(hChargeName,"Graph_Charges_layer%d_chip%d_chan%d",l,i,j);
                gCharge4[l][i][j] = new TGraph();
                gCharge4[l][i][j]->SetNameTitle(hChargeName,hChargeName);
                Hlist4->Add(gCharge4[l][i][j]);                
                sprintf(hChargeName,"Graph_Times_layer%d_chip%d_chan%d",l,i,j);
                gTime4[l][i][j] = new TGraph();
                gTime4[l][i][j]->SetNameTitle(hChargeName,hChargeName);
                Hlist4->Add(gTime4[l][i][j]);

                gCharge4Count[l][i][j]=0;
                gTime4Count[l][i][j]=0;

                sprintf(hChargeName,"Graph_Charges_layer%d_chip%d_chan%d_fit",l,i,j);
                gCharge4_fit[l][i][j] = new TGraph();
                gCharge4_fit[l][i][j]->SetNameTitle(hChargeName,hChargeName);
                sprintf(hChargeName,"Graph_Times_layer%d_chip%d_chan%d_fit",l,i,j);
                gTime4_fit[l][i][j] = new TGraph();
                gTime4_fit[l][i][j]->SetNameTitle(hChargeName,hChargeName);
            }
            //sprintf(hChargeName,"Charges_layer%d_chip%d",l,i);
            //hCharge3[l][i] = new TH1F(hChargeName,hChargeName,4096,0,4095);
            //Hlist3->Add(hCharge3[l][i]);
            //sprintf(hChargeName,"Times_layer%d_chip%d",l,i);
            //hTime3[l][i] = new TH1F(hChargeName,hChargeName,4096,0,4095);
            //Hlist3->Add(hTime3[l][i]);
        }
    }

    // initial the brach of EventAlign
    _cellID=0;
    _bcid=0;
    _hitTag=0;
    _gainTag=0;
    _charge=0;
    _time=0;

    // draw some plots for check
    hh1 = new TH2F("hh1","",4096,0,4096,4096,0,4096);
    hh2 = new TH2F("hh2","",5000,0,80096,4096,0,4096);

    for(int i=0; i!=chipNu; ++i)
    {
        char hname[50];
        sprintf(hname,"time_%d",i);
        hPed[i] = new TH1F(hname,"",1000,0,1000);
    }
}

void treeTemperature::Reset()
{
    _pedCellID.clear();
    _pedMeanCharges.clear();
    _pedSigmaCharges.clear();
    _pedMeanTimes.clear();
    _pedSigmaTimes.clear();
    _pedTempSlopeCharges.clear();
    _pedTempInterceptCharges.clear();
    _pedTempSlopeTimes.clear();
    _pedTempInterceptTimes.clear();
}

treeTemperature::~treeTemperature()
{
    //for(int l=0; l!=layerNu; ++l)
    //{
    //    for(int i=0; i!=chipNu; ++i)
    //    {
    //        for(int j=0; j!=channelNu; ++j)
    //        {
    //            for(int k=0; k!=memoryNu; ++k)
    //            {
    //                hCharge1[l][i][j][k]->Delete();
    //                hTime1[l][i][j][k] ->Delete();
    //            }
    //            hCharge2[l][i][j]->Delete();
    //            hTime2[l][i][j]->Delete();
    //        }
    //        hCharge3[l][i]->Delete();
    //        hTime3[l][i]->Delete();
    //    }
    //}
}

int treeTemperature::Pedestal(string root_name,string pedestal_name,string fitdata_name)
{
    Init();
    string txtname(root_name);
    txtname.erase(0,txtname.find_last_of("/")+1);

    // *********************************************************************************
    //TFile *f =  TFile::Open("/mnt2/ScECAL_CR/results/eventAlign/20201216_1924_cosmicRayTst.root","READ");
    //TFile *f =  TFile::Open("test.root","RECREATE");
    //rawfile =  TFile::Open(root_name.c_str(),"READ");
    rawfile =  TFile::Open(root_name.c_str(),"READ");
    if(!rawfile){ cout<<" Get open eventAlign file failed !!!"<<endl;return 0;}
    rawfile->GetObject("T_Event",tree);
    if(!tree){ cout<<" Get tree object failed !!!"<<endl;return 0;}
    int Entries = tree->GetEntries();
    cout<<" [ Prepare to PedestalExtraction : ] Total Entries : "<<Entries<<" "<<endl;
    tree ->SetBranchAddress("cycleID",&_cycleID);
    tree ->SetBranchAddress("cellIDs",&_cellID);
    tree ->SetBranchAddress("BCIDs",&_bcid);
    tree ->SetBranchAddress("hitTags",&_hitTag);
    tree ->SetBranchAddress("gainTags",&_gainTag);
    tree ->SetBranchAddress("charges",&_charge);
    tree ->SetBranchAddress("times",&_time);

    // *********************************************************************************
    pedestal_name.insert(pedestal_name.length(),txtname);
    fitdata_name.insert(fitdata_name.length(), txtname);

    // ==============================================================================
    //                          Extracting pedestal 
    // ==============================================================================
    for(int entry=0; entry!=Entries; ++entry)
    {
        tree->GetEntry(entry);
        if(entry%10000==0)   
            cout<<" [ Process : PedestalExtraction : ]  event : "<<entry<<" "<<endl;

        // int cycleIDs = _cycleID;
        for(unsigned i_hit=0; i_hit!=_bcid->size(); ++i_hit)
        {
            int layerID = _cellID->at(i_hit)/1e6;
            int chipIDs = _cellID->at(i_hit)%1000000/1e4;
            int memoryCells = _cellID->at(i_hit)%10000/1e2;
            int channelIDs = _cellID->at(i_hit)%100;

            int  charge = _charge->at(i_hit);
            int  time   = _time->at(i_hit);

            // only save hit=0 events, for pedestal by auto-trigger
            if(_hitTag->at(i_hit)!=0) continue;

            if(1)
            {
                // memoryCells Level
               // hCharge1[layerID][chipIDs][channelIDs][memoryCells] ->Fill(charge);
               // hTime1[layerID][chipIDs][channelIDs][memoryCells] ->Fill(time);

                // channel Level
                hCharge2[layerID][chipIDs][channelIDs] ->Fill(charge);
                hTime2[layerID][chipIDs][channelIDs] ->Fill(time);

                // chip Level
                //hCharge3[layerID][chipIDs] ->Fill(charge);
                //hTime3[layerID][chipIDs] ->Fill(time);
            }
        }
    }


    // ==============================================================================
    //                          Fitting pedestal 
    // ==============================================================================

    // *********************************************************************************

    chargefun = new TF1("chargefun","gaus");
    timefun = new TF1("timefun","gaus");
    cout<<fitdata_name<<endl;
    fitfile = new TFile(fitdata_name.c_str(),"RECREATE");
    // **************************** Memory Cells Level *****************************************************
    fNtuple1 = new TTree("CellLevel","fit pedestal at memoryCell Level");
    /*
    cout<<" ============== Begin memoryCell Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple1 ->Branch("cellIDs",&_pedCellID);
    fNtuple1 ->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple1 ->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple1 ->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple1 ->Branch("pedSigmaTimes",&_pedSigmaTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                for(int k=0; k!=memoryNu; ++k)
                    //for(int k=0; k!=1; ++k) // for CR test only one memory used
                {
                    if(hCharge1[l][i][j][k]->GetEntries()==0) continue;
                    _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+k*pow(10,2)+j);

                    float  meanCharges = hCharge1[l][i][j][k]->GetMean();
                    float  rmsCharges = hCharge1[l][i][j][k]->GetRMS();
                    hCharge1[l][i][j][k] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                    meanCharges = chargefun->GetParameter(1);
                    rmsCharges = chargefun->GetParameter(2);
                    hCharge1[l][i][j][k] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                    meanCharges = chargefun->GetParameter(1);
                    rmsCharges = chargefun->GetParameter(2);
                    _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                    _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                    hCharge1[l][i][j][k] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                    _pedMeanCharges.push_back(meanCharges);
                    _pedSigmaCharges.push_back(rmsCharges);

                    if(hTime1[l][i][j][k]->GetEntries()==0) continue;

                    float  meanTimes = hTime1[l][i][j][k]->GetMean();
                    float  rmsTimes = hTime1[l][i][j][k]->GetRMS();
                    hTime1[l][i][j][k] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                    meanTimes = chargefun->GetParameter(1);
                    rmsTimes = chargefun->GetParameter(2);
                    hTime1[l][i][j][k] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                    meanTimes = chargefun->GetParameter(1);
                    rmsTimes = chargefun->GetParameter(2);
                    _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                    _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                    hTime1[l][i][j][k] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                    _pedMeanTimes.push_back(meanTimes);
                    _pedSigmaTimes.push_back(rmsTimes);
                }
            }
        }
    }
    fNtuple1->Fill();
    cout<<" ============== Finish memoryCell Level pedestal fitting Process==============="<<endl;
    */

    // **************************** Channel Cells Level *****************************************************
    cout<<" ============== Begin Channel Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple2 = new TTree("ChnLevel","fit pedestal at memoryCell Level");
    fNtuple2->Branch("cellIDs",&_pedCellID);
    fNtuple2->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple2->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple2->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple2->Branch("pedSigmaTimes",&_pedSigmaTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                // Channel Level
                if(hCharge2[l][i][j]->GetEntries()==0) continue;
                if(hTime2[l][i][j]->GetEntries()==0) continue;
                _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+0*pow(10,2)+j);

				//  The ADC distribution of hits with hitTag == 0 shows several small peaks in larger region.
				//  So the fit based on the gaussian mean and RMS works not well for the beam test data.
				/*
                float  meanCharges = hCharge2[l][i][j]->GetMean();
                float  rmsCharges = hCharge2[l][i][j]->GetRMS();
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);

                float  meanTimes = hTime2[l][i][j]->GetMean();
                float  rmsTimes = hTime2[l][i][j]->GetRMS();
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                hTime2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);
				*/
				
				// To address the issue of small peaks, we only fit the left 80% histogram.

				double qtl_pers[3] = {0., 0.5, 0.80};
				double qtl[] = {380 ,500,1000};
				hCharge2[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				double rangeLX = qtl[0];
				double rangeRX = qtl[2];
				hCharge2[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                double meanCharges = chargefun->GetParameter(1);
                double rmsCharges = chargefun->GetParameter(2);
                _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                _pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);

				hTime2[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				rangeLX = qtl[0];
				rangeRX = qtl[2];
				hTime2[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                double meanTimes = chargefun->GetParameter(1);
                double rmsTimes = chargefun->GetParameter(2);
                _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                hTime2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                _pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);
				
            }
        }
    }
    fNtuple2->Fill();
    cout<<" ============== Finish Channel Level pedestal fitting Process==============="<<endl;
    

    // **************************** Chip Cells Level *****************************************************
    fNtuple3 = new TTree("ChipLevel","fit pedestal at memoryCell Level");
   /* cout<<" ============== Begin Chip Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple3 ->Branch("cellIDs",&_pedCellID);
    fNtuple3 ->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple3 ->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple3 ->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple3 ->Branch("pedSigmaTimes",&_pedSigmaTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            // Chip Level
            if(hCharge3[l][i]->GetEntries()==0) continue;
            _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+0*pow(10,2)+0);

            float  meanCharges = hCharge3[l][i]->GetMean();
            float  rmsCharges = hCharge3[l][i]->GetRMS();
            hCharge3[l][i] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
            meanCharges = chargefun->GetParameter(1);
            rmsCharges = chargefun->GetParameter(2);
            hCharge3[l][i] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
            meanCharges = chargefun->GetParameter(1);
            rmsCharges = chargefun->GetParameter(2);
            _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
            _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
            hCharge3[l][i] ->GetXaxis()->SetRangeUser(_leftC,_righC);

            _pedMeanCharges.push_back(meanCharges);
            _pedSigmaCharges.push_back(rmsCharges);
        }
    }
    fNtuple3->Fill();
    cout<<" ============== Finish Chip Level pedestal fitting Process==============="<<endl;
*/
    // *********************************************************************************
    fNtuple1->Write();
    fNtuple2->Write();
    fNtuple3->Write();
    fitfile->Close();

    //    TCanvas* c = new TCanvas("c","",10,10,1500,800);
    //    c->Divide(3,2);
    //    for(int tt=0; tt!=6; tt++)
    //    {
    //        c->cd(tt+1);
    //        hCharge3[0][tt]->Draw();
    //    }
    // **********************************************************************************

    TFile *pedestalfile = new TFile(pedestal_name.c_str(),"RECREATE");
    //pedestalfile ->mkdir("MemCellLevel");
    //pedestalfile ->cd("MemCellLevel");
    //Hlist1 ->Write();
    //char char_tmp[1000];
    //for(int l=0; l!=layerNu; ++l)
    //{   
    //    for(int i=0; i!=chipNu; ++i)
    //    {
    //        pedestalfile->cd("MemCellLevel");
    //        sprintf(char_tmp,"layer%d/chip%d",l,i);
    //        pedestalfile->mkdir(char_tmp);
    //        pedestalfile->cd(char_tmp);
    //        for(int j=0; j!=channelNu; ++j)
    //        {
    //            for(int k=0; k!=memoryNu; ++k)
    //            {
    //                hTime1[l][i][j][k]->Write();          
    //                hCharge1[l][i][j][k]->Write();          
    //            }
    //        }
    //    }
    //}
    pedestalfile ->mkdir("ChannelLevel");
    pedestalfile ->cd("ChannelLevel");
    Hlist2 ->Write();
   // pedestalfile ->mkdir("ChipLevel");
   // pedestalfile ->cd("ChipLevel");
   // Hlist3 ->Write();
    pedestalfile->Close();
    
    cout<<"Ending of the process!!!!"<<endl;
}

void treeTemperature::TemperatureDependence(string root_name,string temperature_name,string fitdata_name)
{
    Init();
    string txtname(root_name);
    txtname.erase(0,txtname.find_last_of("/")+1);

    // *********************************************************************************
    //TFile *f =  TFile::Open("/mnt2/ScECAL_CR/results/eventAlign/20201216_1924_cosmicRayTst.root","READ");
    //TFile *f =  TFile::Open("test.root","RECREATE");
    //rawfile =  TFile::Open(root_name.c_str(),"READ");
    rawfile =  TFile::Open(root_name.c_str(),"READ");
    if(!rawfile){ cout<<" Get open eventAlign file failed !!!"<<endl;return;}
    rawfile->GetObject("T_Event",tree);
    if(!tree){ cout<<" Get tree object failed !!!"<<endl;return;}
    int Entries = tree->GetEntries();
    cout<<" [ Prepare to PedestalExtraction : ] Total Entries : "<<Entries<<" "<<endl;
    tree ->SetBranchAddress("cycleID",&_cycleID);
    tree ->SetBranchAddress("cellIDs",&_cellID);
    tree ->SetBranchAddress("BCIDs",&_bcid);
    tree ->SetBranchAddress("hitTags",&_hitTag);
    tree ->SetBranchAddress("gainTags",&_gainTag);
    tree ->SetBranchAddress("charges",&_charge);
    tree ->SetBranchAddress("times",&_time);
    tree ->SetBranchAddress("temp",&_tempLayer);

    // *********************************************************************************
    temperature_name.insert(temperature_name.length(),txtname);
    fitdata_name.insert(fitdata_name.length(), txtname);

    // ==============================================================================
    //                          Extracting pedestal 
    // ==============================================================================
    for(int entry=0; entry!=Entries; ++entry)
    {
        tree->GetEntry(entry);
        if(entry%10000==0)   
            cout<<" [ Process : PedestalExtraction : ]  event : "<<entry<<" "<<endl;

        // int cycleIDs = _cycleID;
        for(unsigned i_hit=0; i_hit!=_bcid->size(); ++i_hit)
        {
            int layerID = _cellID->at(i_hit)/1e6;
            int chipIDs = _cellID->at(i_hit)%1000000/1e4;
            int memoryCells = _cellID->at(i_hit)%10000/1e2;
            int channelIDs = _cellID->at(i_hit)%100;

            int  charge = _charge->at(i_hit);
            int  time   = _time->at(i_hit);

            // only save hit=0 events, for pedestal by auto-trigger
            if(_hitTag->at(i_hit)!=0) continue;

            if(1)
            {
                // memoryCells Level
               // hCharge1[layerID][chipIDs][channelIDs][memoryCells] ->Fill(charge);
               // hTime1[layerID][chipIDs][channelIDs][memoryCells] ->Fill(time);

                // channel Level
                hCharge2[layerID][chipIDs][channelIDs] ->Fill(charge);
                hTime2[layerID][chipIDs][channelIDs] ->Fill(time);

                // chip Level
                //hCharge3[layerID][chipIDs] ->Fill(charge);
                //hTime3[layerID][chipIDs] ->Fill(time);

                // channel Level
                hCharge4[layerID][chipIDs][channelIDs] ->Fill(charge);
                hTime4[layerID][chipIDs][channelIDs] ->Fill(time);

                double *_position = EBUdecode(layerID,chipIDs,channelIDs);
                double SiPMtemp;
                if(_tempLayer->at(layerID).size()!=0) SiPMtemp = tempReconstruction(layerID, _position, _tempLayer->at(layerID));
                else SiPMtemp=20;

                gCharge4[layerID][chipIDs][channelIDs]->SetPoint(gCharge4Count[layerID][chipIDs][channelIDs]++,SiPMtemp,charge);
                gTime4[layerID][chipIDs][channelIDs]->SetPoint(gTime4Count[layerID][chipIDs][channelIDs]++,SiPMtemp,time);
            }
        }
    }


    // ==============================================================================
    //                          Fitting pedestal 
    // ==============================================================================

    // *********************************************************************************

    chargefun = new TF1("chargefun","gaus");
    timefun = new TF1("timefun","gaus");
    gchargefun = new TF1("gchargefun","[0]+[1]*x",10,35);
    gtimefun = new TF1("gtimefun","[0]+[1]*x",10,35);
    cout<<fitdata_name<<endl;
    fitfile = new TFile(fitdata_name.c_str(),"RECREATE");
    // **************************** Memory Cells Level *****************************************************
    fNtuple1 = new TTree("CellLevel","fit pedestal at memoryCell Level");
    /*
    cout<<" ============== Begin memoryCell Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple1 ->Branch("cellIDs",&_pedCellID);
    fNtuple1 ->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple1 ->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple1 ->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple1 ->Branch("pedSigmaTimes",&_pedSigmaTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                for(int k=0; k!=memoryNu; ++k)
                    //for(int k=0; k!=1; ++k) // for CR test only one memory used
                {
                    if(hCharge1[l][i][j][k]->GetEntries()==0) continue;
                    _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+k*pow(10,2)+j);

                    float  meanCharges = hCharge1[l][i][j][k]->GetMean();
                    float  rmsCharges = hCharge1[l][i][j][k]->GetRMS();
                    hCharge1[l][i][j][k] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                    meanCharges = chargefun->GetParameter(1);
                    rmsCharges = chargefun->GetParameter(2);
                    hCharge1[l][i][j][k] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                    meanCharges = chargefun->GetParameter(1);
                    rmsCharges = chargefun->GetParameter(2);
                    _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                    _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                    hCharge1[l][i][j][k] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                    _pedMeanCharges.push_back(meanCharges);
                    _pedSigmaCharges.push_back(rmsCharges);

                    if(hTime1[l][i][j][k]->GetEntries()==0) continue;

                    float  meanTimes = hTime1[l][i][j][k]->GetMean();
                    float  rmsTimes = hTime1[l][i][j][k]->GetRMS();
                    hTime1[l][i][j][k] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                    meanTimes = chargefun->GetParameter(1);
                    rmsTimes = chargefun->GetParameter(2);
                    hTime1[l][i][j][k] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                    meanTimes = chargefun->GetParameter(1);
                    rmsTimes = chargefun->GetParameter(2);
                    _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                    _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                    hTime1[l][i][j][k] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                    _pedMeanTimes.push_back(meanTimes);
                    _pedSigmaTimes.push_back(rmsTimes);
                }
            }
        }
    }
    fNtuple1->Fill();
    cout<<" ============== Finish memoryCell Level pedestal fitting Process==============="<<endl;
    */

    // **************************** Channel Cells Level *****************************************************
    cout<<" ============== Begin Channel Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple2 = new TTree("ChnLevel","fit pedestal at memoryCell Level");
    fNtuple2->Branch("cellIDs",&_pedCellID);
    fNtuple2->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple2->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple2->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple2->Branch("pedSigmaTimes",&_pedSigmaTimes);
    const int nsigma = 4;
    double chargeadcmax[layerNu][chipNu][channelNu],chargeadcmin[layerNu][chipNu][channelNu],timeadcmax[layerNu][chipNu][channelNu],timeadcmin[layerNu][chipNu][channelNu];
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                // Channel Level
                if(hCharge2[l][i][j]->GetEntries()==0) {
                    chargeadcmax[l][i][j]=4095, chargeadcmin[l][i][j]=0;
                    continue;
                }
                if(hTime2[l][i][j]->GetEntries()==0) {
                    timeadcmax[l][i][j]=4095, timeadcmin[l][i][j]=0;
                    continue;
                }
                _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+0*pow(10,2)+j);

				//  The ADC distribution of hits with hitTag == 0 shows several small peaks in larger region.
				//  So the fit based on the gaussian mean and RMS works not well for the beam test data.
				/*
                float  meanCharges = hCharge2[l][i][j]->GetMean();
                float  rmsCharges = hCharge2[l][i][j]->GetRMS();
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);

                float  meanTimes = hTime2[l][i][j]->GetMean();
                float  rmsTimes = hTime2[l][i][j]->GetRMS();
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                hTime2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);
				*/
				
				// To address the issue of small peaks, we only fit the left 80% histogram.

				double qtl_pers[3] = {0., 0.5, 0.80};
				double qtl[] = {380 ,500,1000};
				hCharge2[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				double rangeLX = qtl[0];
				double rangeRX = qtl[2];
				hCharge2[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                double meanCharges = chargefun->GetParameter(1);
                double rmsCharges = chargefun->GetParameter(2);
                _leftC = (meanCharges-nsigma*rmsCharges)>0?(meanCharges-nsigma*rmsCharges):0;
                _righC = (meanCharges+nsigma*rmsCharges)<4096?(meanCharges+nsigma*rmsCharges):4096;
                hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                _pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);
                chargeadcmax[l][i][j] = _righC;
                chargeadcmin[l][i][j] = _leftC;

				hTime2[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				rangeLX = qtl[0];
				rangeRX = qtl[2];
				hTime2[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                double meanTimes = chargefun->GetParameter(1);
                double rmsTimes = chargefun->GetParameter(2);
                _leftC = (meanTimes-nsigma*rmsTimes)>0?(meanTimes-nsigma*rmsTimes):0;
                _righC = (meanTimes+nsigma*rmsTimes)<4096?(meanTimes+nsigma*rmsTimes):4096;
                hTime2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                _pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);
                timeadcmax[l][i][j] = _righC;
                timeadcmin[l][i][j] = _leftC;
				
            }
        }
    }
    fNtuple2->Fill();
    cout<<" ============== Finish Channel Level pedestal fitting Process==============="<<endl;
    

    // **************************** Chip Cells Level *****************************************************
    fNtuple3 = new TTree("ChipLevel","fit pedestal at memoryCell Level");
   /* cout<<" ============== Begin Chip Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple3 ->Branch("cellIDs",&_pedCellID);
    fNtuple3 ->Branch("pedMeanCharges",&_pedMeanCharges);
    fNtuple3 ->Branch("pedSigmaCharges",&_pedSigmaCharges);
    fNtuple3 ->Branch("pedMeanTimes",&_pedMeanTimes);
    fNtuple3 ->Branch("pedSigmaTimes",&_pedSigmaTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            // Chip Level
            if(hCharge3[l][i]->GetEntries()==0) continue;
            _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+0*pow(10,2)+0);

            float  meanCharges = hCharge3[l][i]->GetMean();
            float  rmsCharges = hCharge3[l][i]->GetRMS();
            hCharge3[l][i] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
            meanCharges = chargefun->GetParameter(1);
            rmsCharges = chargefun->GetParameter(2);
            hCharge3[l][i] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
            meanCharges = chargefun->GetParameter(1);
            rmsCharges = chargefun->GetParameter(2);
            _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
            _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
            hCharge3[l][i] ->GetXaxis()->SetRangeUser(_leftC,_righC);

            _pedMeanCharges.push_back(meanCharges);
            _pedSigmaCharges.push_back(rmsCharges);
        }
    }
    fNtuple3->Fill();
    cout<<" ============== Finish Chip Level pedestal fitting Process==============="<<endl;
*/


    // **************************** extraxting pedestal only in fist peak for fit *****************************************************
    int gCharge4Count_fit[layerNu][chipNu][channelNu], gTime4Count_fit[layerNu][chipNu][channelNu];
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                gCharge4Count_fit[l][i][j]=0;
                gTime4Count_fit[l][i][j]=0;
            }
        }
    }
    for(int entry=0; entry!=Entries; ++entry)
    {
        tree->GetEntry(entry);
        if(entry%10000==0)   
            cout<<" [ Process : PedestalExtraction : ]  event : "<<entry<<" "<<endl;

        // int cycleIDs = _cycleID;
        for(unsigned i_hit=0; i_hit!=_bcid->size(); ++i_hit)
        {
            int layerID = _cellID->at(i_hit)/1e6;
            int chipIDs = _cellID->at(i_hit)%1000000/1e4;
            int memoryCells = _cellID->at(i_hit)%10000/1e2;
            int channelIDs = _cellID->at(i_hit)%100;

            int  charge = _charge->at(i_hit);
            int  time   = _time->at(i_hit);

            // only save hit=0 events, for pedestal by auto-trigger
            if(_hitTag->at(i_hit)!=0) continue;

            if(1)
            {
                double *_position = EBUdecode(layerID,chipIDs,channelIDs);
                double SiPMtemp;
                if(_tempLayer->at(layerID).size()!=0) SiPMtemp = tempReconstruction(layerID, _position, _tempLayer->at(layerID));
                else SiPMtemp=20;

                if(charge>=chargeadcmin[layerID][chipIDs][channelIDs] && charge<=chargeadcmax[layerID][chipIDs][channelIDs]) gCharge4_fit[layerID][chipIDs][channelIDs]->SetPoint(gCharge4Count_fit[layerID][chipIDs][channelIDs]++,SiPMtemp,charge);
                if(time>=timeadcmin[layerID][chipIDs][channelIDs] && time<=timeadcmax[layerID][chipIDs][channelIDs]) gTime4_fit[layerID][chipIDs][channelIDs]->SetPoint(gTime4Count_fit[layerID][chipIDs][channelIDs]++,SiPMtemp,time);
            }
        }
    }


    // **************************** Channel Cells Level *****************************************************
    // temperature
    cout<<" ============== Begin Channel Level temperature vs pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple4 = new TTree("ChnLevel","fit temperature vs pedestal at channel Level");
    fNtuple4->Branch("cellIDs",&_pedCellID);
    fNtuple4->Branch("pedTempSlopeCharges",&_pedTempSlopeCharges);
    fNtuple4->Branch("pedTempInterceptCharges",&_pedTempInterceptCharges);
    fNtuple4->Branch("pedTempSlopeTimes",&_pedTempSlopeTimes);
    fNtuple4->Branch("pedTempInterceptTimes",&_pedTempInterceptTimes);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                // Channel Level
                if(hCharge4[l][i][j]->GetEntries()==0) continue;
                if(hTime4[l][i][j]->GetEntries()==0) continue;
                _pedCellID.push_back(l*pow(10,6)+i*pow(10,4)+0*pow(10,2)+j);

				//  The ADC distribution of hits with hitTag == 0 shows several small peaks in larger region.
				//  So the fit based on the gaussian mean and RMS works not well for the beam test data.
				/*
                float  meanCharges = hCharge2[l][i][j]->GetMean();
                float  rmsCharges = hCharge2[l][i][j]->GetRMS();
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                hCharge2[l][i][j] ->Fit(chargefun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = chargefun->GetParameter(1);
                rmsCharges = chargefun->GetParameter(2);
                _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
                hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);

                float  meanTimes = hTime2[l][i][j]->GetMean();
                float  rmsTimes = hTime2[l][i][j]->GetRMS();
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                hTime2[l][i][j] ->Fit(chargefun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = chargefun->GetParameter(1);
                rmsTimes = chargefun->GetParameter(2);
                _leftC = (meanTimes-10*rmsTimes)>0?(meanTimes-10*rmsTimes):0;
                _righC = (meanTimes+10*rmsTimes)<4096?(meanTimes+10*rmsTimes):4096;
                hTime2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);

                _pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);
				*/
				
				// To address the issue of small peaks, we only fit the left 80% histogram.

				// double qtl_pers[3] = {0., 0.5, 0.80};
				// double qtl[] = {380 ,500,1000};
				// hCharge4[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				// double rangeLX = qtl[0];
				// double rangeRX = qtl[2];
				// hCharge4[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                // double meanCharges = chargefun->GetParameter(1);
                // double rmsCharges = chargefun->GetParameter(2);
                // _leftC = (meanCharges-nsigma*rmsCharges)>0?(meanCharges-nsigma*rmsCharges):0;
                // _righC = (meanCharges+nsigma*rmsCharges)<4096?(meanCharges+nsigma*rmsCharges):4096;
                // // hCharge4[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                // _pedMeanCharges.push_back(meanCharges);
                // _pedSigmaCharges.push_back(rmsCharges);
                // // gCharge4[l][i][j] ->GetYaxis()->SetRangeUser(_leftC,_righC);
                // gchargefun->SetRange(10,35,_leftC,_righC);
                gCharge4[l][i][j]->Fit(gchargefun, "Q", "", 10,35);
                double slopeCharges = gchargefun->GetParameter(1);
                double interceptCharges = gchargefun->GetParameter(0);
                _pedTempSlopeCharges.push_back(slopeCharges);
                _pedTempInterceptCharges.push_back(interceptCharges);


				// hTime4[l][i][j]->GetQuantiles(3, qtl, qtl_pers);
				// rangeLX = qtl[0];
				// rangeRX = qtl[2];
				// hTime4[l][i][j]->Fit(chargefun, "Q", "", rangeLX, rangeRX);
                // double meanTimes = chargefun->GetParameter(1);
                // double rmsTimes = chargefun->GetParameter(2);
                // _leftC = (meanTimes-nsigma*rmsTimes)>0?(meanTimes-nsigma*rmsTimes):0;
                // _righC = (meanTimes+nsigma*rmsTimes)<4096?(meanTimes+nsigma*rmsTimes):4096;
                // // hTime4[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                // _pedMeanTimes.push_back(meanTimes);
                // _pedSigmaTimes.push_back(rmsTimes);
                // // gTime4[l][i][j] ->GetYaxis()->SetRangeUser(_leftC,_righC);
                // gtimefun->SetRange(10,35,_leftC,_righC);
                gTime4[l][i][j]->Fit(gtimefun, "Q", "", 10,35);
                double slopeTimes = gtimefun->GetParameter(1);
                double interceptTimes = gtimefun->GetParameter(0);
                _pedTempSlopeTimes.push_back(slopeTimes);
                _pedTempInterceptTimes.push_back(interceptTimes);
            }
        }
    }
    fNtuple4->Fill();
    cout<<" ============== Finish Channel Level temperature vs pedestal fitting  Process==============="<<endl;


    // *********************************************************************************

    
    // fNtuple1->Write();
    // fNtuple2->Write();
    // fNtuple3->Write();
    fNtuple4->Write();
    fitfile->Close();

    //    TCanvas* c = new TCanvas("c","",10,10,1500,800);
    //    c->Divide(3,2);
    //    for(int tt=0; tt!=6; tt++)
    //    {
    //        c->cd(tt+1);
    //        hCharge3[0][tt]->Draw();
    //    }
    // **********************************************************************************

    TFile *temperaturefile = new TFile(temperature_name.c_str(),"RECREATE");
    //temperaturefile ->mkdir("MemCellLevel");
    //temperaturefile ->cd("MemCellLevel");
    //Hlist1 ->Write();
    //char char_tmp[1000];
    //for(int l=0; l!=layerNu; ++l)
    //{   
    //    for(int i=0; i!=chipNu; ++i)
    //    {
    //        temperaturefile->cd("MemCellLevel");
    //        sprintf(char_tmp,"layer%d/chip%d",l,i);
    //        temperaturefile->mkdir(char_tmp);
    //        temperaturefile->cd(char_tmp);
    //        for(int j=0; j!=channelNu; ++j)
    //        {
    //            for(int k=0; k!=memoryNu; ++k)
    //            {
    //                hTime1[l][i][j][k]->Write();          
    //                hCharge1[l][i][j][k]->Write();          
    //            }
    //        }
    //    }
    //}
    temperaturefile ->mkdir("pedestalChannelLevel");
    temperaturefile ->cd("pedestalChannelLevel");
    Hlist2 ->Write();
    // temperaturefile ->mkdir("ChipLevel");
    // temperaturefile ->cd("ChipLevel");
    // Hlist3 ->Write();
    // temperaturefile->cd("..");
    temperaturefile ->mkdir("tempChannelLevel");
    temperaturefile ->cd("tempChannelLevel");
    Hlist4 ->Write();
    temperaturefile->Close();
    
    cout<<"Ending of the process!!!!"<<endl;
}
