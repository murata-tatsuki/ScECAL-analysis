/////////////////////////////////////////////////////
//////////////////				/////////////////////
//	save TCanvas instead of TH2F for each chn
//	TH2F has no fit curve in single-channel hist
//	but has fit curve in chip-level canvas
//////////////////				/////////////////////
/////////////////////////////////////////////////////
#include<fstream>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TGraph.h>
#include<TGraph.h>
#include<TRandom.h>
#include<TSpectrum.h>
#include<TVirtualFitter.h>
#include<TH2F.h>
#include"InterCalib.h"


using namespace std;

void interCalib::Init()
{
	char fname[50];
	for (int i_layer = 0; i_layer < LayerNo; ++i_layer)
    {
        for (int i_chip = 0; i_chip < chipNu; ++i_chip)
        {
            for (int i_chn = 0; i_chn < channelNu; ++i_chn)
            {
				sprintf(fname,"Layer%d_Chip%d_Chn%d",i_layer,i_chip,i_chn);
				hisg[i_layer][i_chip][i_chn] = new TH2F("hl_ratio","hl_ratio",1000,200,4000,200,300,700);
				hisg[i_layer][i_chip][i_chn]->SetName(fname);
				hisg[i_layer][i_chip][i_chn]->SetTitle(fname);
				hisg[i_layer][i_chip][i_chn]->GetXaxis()->SetRangeUser(0,3800);
				hisg[i_layer][i_chip][i_chn]->GetXaxis()->SetTitle("HighGain [ADC]");
				hisg[i_layer][i_chip][i_chn]->GetYaxis()->SetRangeUser(350,650);
				hisg[i_layer][i_chip][i_chn]->GetYaxis()->SetTitle("LowGain [ADC]");
				hisg[i_layer][i_chip][i_chn]->SetOption("colz");
				eventNo[i_layer][i_chip][i_chn]=0;
				Xmax[i_layer][i_chip][i_chn]=0;
			}
		}
	}

}

interCalib::interCalib(string root_name,string sps_name)
{
    Init();
    double start = clock();
    // ================ Read DATA =======================    
    //string dataname(root_name);
    TFile* infile = TFile::Open(root_name.c_str(),"READ");
    if(!infile)cout<<"read file failed"<<endl;
    TTree* mipTreeRead = (TTree*)infile->Get("Raw_Hit");
    if(!mipTreeRead)cout<<"Get tree failed"<<endl;
    //double _realLength;
    //mipTreeRead->SetBranchAddress("realLength",&_realLength);
    mipTreeRead->SetBranchAddress("CellID",&_cellID);
    mipTreeRead->SetBranchAddress("HitTag",&_hitTag);
    mipTreeRead->SetBranchAddress("LG_Charge",&_lg_charge);
    mipTreeRead->SetBranchAddress("HG_Charge",&_hg_charge);

    for(int entry=0; entry<mipTreeRead->GetEntries(); ++entry)
    {
        mipTreeRead->GetEntry(entry);
        if(entry%10000==0)
            cout<<" [ Process ]  event : "<< entry << endl;
        for(unsigned i_hit=0; i_hit!=_cellID->size(); ++i_hit)
        {
			if(_hitTag->at(i_hit)==0) continue;
            int layerID = _cellID->at(i_hit)/1e5;
            int chipID = _cellID->at(i_hit)%100000/1e4;
            int chnID = _cellID->at(i_hit)%100;
			//if(layerID==0 && chipID==0 && chnID==0) cout << _hitTag->at(i_hit) << endl;
            double lg_charge = _lg_charge->at(i_hit);
            double hg_charge = _hg_charge->at(i_hit);
			hisg[layerID][chipID][chnID]->Fill(hg_charge,lg_charge);
			eventNo[layerID][chipID][chnID]++;
			Xmax[layerID][chipID][chnID] = Xmax[layerID][chipID][chnID] > hg_charge ? Xmax[layerID][chipID][chnID] : hg_charge;
        }
    }
    infile->Close();

    // ================ LINEAR FITTING =======================    
    //string txtname(root_name);
    //txtname.erase(txtname.find('.'));
    //txtname.insert(txtname.length(),"_InterCalib.root");

    TFile* spsfile = new TFile(sps_name.c_str(),"RECREATE");
    TTree* fNtuple = new TTree("InterCalib","fit SPS with multigaus");
    fNtuple->Branch("cellID",&_channelID); //i_chn + i_chip*36 + i_layer*210
    fNtuple->Branch("slope",&_slope);
    fNtuple->Branch("intercept",&_intercept);
    fNtuple->Branch("slopeError",&_slopeError);
    fNtuple->Branch("interceptError",&_interceptError);
    fNtuple->Branch("chiSquare",&_chiSquare);
    fNtuple->Branch("NDF",&_NDF);

    char fname[50];
    _BlindChannel=0;
    gStyle->SetOptFit(1111);
	int Xmax_Fit;
    for (int i_layer = 0; i_layer < LayerNo; ++i_layer)
    //for (int i_layer = 27; i_layer < 28; ++i_layer)
    {
        cout<<" -------- "<<i_layer<<" -------- "<<endl;
        for (int i_chip = 0; i_chip < chipNu; ++i_chip)
        {
    		TCanvas* c1 = new TCanvas(Form("c_%d_%d",i_layer,i_chip),Form("c_%d_%d",i_layer,i_chip),10,10,1500,1200);
			c1->Divide(6,6);
			int pad_no = 0;
            spsfile->cd();
            sprintf(fname,"Layer%d/Chip%d",i_layer,i_chip);
            spsfile->mkdir(fname);
            spsfile->cd(fname);
            for (int i_chn = 0; i_chn < channelNu; ++i_chn)
            {
    			TCanvas* c2 = new TCanvas(Form("c_%d_%d_%d",i_layer,i_chip,i_chn),Form("c_%d_%d_%d",i_layer,i_chip,i_chn),800,600);
				if(Xmax[i_layer][i_chip][i_chn] < 1000 && Xmax[i_layer][i_chip][i_chn] != 0)  
					cout << "Xmax : " << Xmax[i_layer][i_chip][i_chn] <<" layer:" << i_layer <<" chip:" << i_chip <<" chn:"<< i_chn << endl;
				pad_no++;
				if((i_chip==5&&i_chn>29) && pad_no != 36) continue;
				if((i_chip==5&&i_chn>29) && pad_no == 36) {c1->Write(); continue;}
				if(hisg[i_layer][i_chip][i_chn]->GetEntries()<25) 
				{
					cout << "layer:chip:chan	"<< i_layer << " "<< i_chip << " "<< i_chn << endl;
					_BlindChannel++;
					if(pad_no!=36) continue;
					if(pad_no==36) {c1->Write(); continue;}
				}
				//cout<<"--- "<<i_layer<<" : "<<i_chip<<" : "<<i_chn<<" : "<<eventNo[i_layer][i_chip][i_chn]<<endl;
				TF1* fun = new TF1("fun","[0]*x+[1]",500,3000);
				fun->SetParameter(0,0.03);
				fun->SetParameter(1,370);
                fun->SetParLimits(0,0.023,0.037);//27~43.4
                fun->SetParLimits(1,300,500);
                fun->SetLineColor(2);
                fun->SetLineWidth(2);
				Xmax_Fit = Xmax[i_layer][i_chip][i_chn] > 2000 ? Xmax[i_layer][i_chip][i_chn]-500 : 1500 ;
                if(i_layer>3&&i_layer<28)
				{
                	hisg[i_layer][i_chip][i_chn]->Fit(fun,"qn","",600,Xmax_Fit);
                    //fun->SetRange(500,2500);
                }
                else
				{
                    hisg[i_layer][i_chip][i_chn]->Fit(fun,"qn","",600,Xmax_Fit);
                 	//fun->SetRange(500,2500);
                }
                // ==============================================================
                
                c1->cd(pad_no);
                hisg[i_layer][i_chip][i_chn]->Draw();
                fun->Draw("same");
					
				c2->cd();
				hisg[i_layer][i_chip][i_chn]->Draw();
                fun->Draw("same");
				c2->Write();
			
			
			//	if((i_layer==9&&i_chip==1&&i_chn==21)||(i_layer==13&&i_chip==1&&i_chn==32)||(i_layer==25&&i_chip==3&&i_chn==35)||(i_layer==24&&i_chip==1&&i_chn==21)||(i_layer==27&&i_chip==5&&i_chn==22))
			//	{
			//		c2->cd();
            //    	hisg[i_layer][i_chip][i_chn]->Draw();
            //    	fun->Draw("same");
			//		c2->Print(Form("test%d.png",i_layer));
			//	}
				
				if(pad_no==36) {c1->Write();}
                // ==============================================================
                _channelID = i_chn + i_chip*1e4 + i_layer*1e5;
				_slope = fun->GetParameter(0);
                _intercept = fun->GetParameter(1);
                _slopeError = fun->GetParError(0);
                _interceptError = fun->GetParError(1);
                _chiSquare = fun->GetChisquare();
                _NDF = fun->GetNDF();
				
                //cout<<" ---------- "<<chn<<" ------------ "<<_slope<<endl;

                fNtuple->Fill();
            }
        }
    }

    double end = clock();
    printf(" --- total select %d blind chns : Times %.1f s ************* \n",_BlindChannel,(end-start)/CLOCKS_PER_SEC);
    // *********************************************************************************
    spsfile->cd();
    //fNtuple->Write("",TObject::kOverwrite);
    fNtuple->Write();
    spsfile->Close();
}

