
void effEBU()
{
    const int LayerNo=30;
    TGraph* gEff[LayerNo];
    TGraph* gTemp[LayerNo];
    for(int i=0; i!=LayerNo; ++i)
    {
        gEff[i] = new TGraph();
        gTemp[i] = new TGraph();
    }
    // hit chn=0
   // double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,
   //                           2,1,4,2,1,2,3,3,3,1,
   //                           2,4,1,1,0,0,0,0,9,3};
    // hit chn<=10
    double blindChn[LayerNo]={3,1,4,1,1,3,1,1,1,1,
                              2,1,4,2,1,2,3,3,3,1,
                              2,4,1,2,0,1,0,0,9,3};

    TH1F* hEff0 = new TH1F("hEff0","",30,0,30);
    ifstream filelist("./datalist_sps.txt");
    // ifstream filelist("./datalist_ps.txt");
    if(!filelist.good()){cout<<"open datalist file failed !"<<endl; exit(0);}
    vector<string> datalist;
    string dname;
    while(filelist >> dname) datalist.push_back(dname);
    const int nfile = (int)datalist.size();
    if(nfile==0){cout<<"datalist is empty !"<<endl; return;}
    int nbin=nfile;
    double effmean = 0.;
    // ==============================================================
    //for(int j=0; j!=100; j++){
    for(int j=0; j<nfile; j++){
        const string& dname = datalist[j];
        std::vector<int>* _newCell = nullptr;
        std::vector<double>* _charge = nullptr;
        std::vector<double>* _trackFitPars = nullptr;
        std::vector<double>* _residualX = nullptr;
        std::vector<double>* _residualY = nullptr;
        std::vector<double>* _newPosX = nullptr;
        std::vector<double>* _newPosY = nullptr;
        std::vector<double>* _newPosZ = nullptr;
        std::vector<double>* _temp = nullptr;    

        TFile* file = new TFile(dname.c_str());
        TTree* fNtuple1 = (TTree*)file->Get("T_Event");
        int _totalEntries = fNtuple1->GetEntries();
        if(!fNtuple1) cout<<"open RawtoRoot file failed "<<endl;
        fNtuple1 ->SetBranchAddress("hitCellnew",&_newCell);
        fNtuple1 ->SetBranchAddress("energyDep",&_charge);
        fNtuple1 ->SetBranchAddress("trackFitPars",&_trackFitPars);
        fNtuple1 ->SetBranchAddress("residualX",&_residualX);
        fNtuple1 ->SetBranchAddress("residualY",&_residualY);
        fNtuple1 ->SetBranchAddress("hitPosXnew",&_newPosX);
        fNtuple1 ->SetBranchAddress("hitPosYnew",&_newPosY);
        fNtuple1 ->SetBranchAddress("hitPosZnew",&_newPosZ);
        fNtuple1 ->SetBranchAddress("Temperature",&_temp);

        TH1F* htemp[LayerNo];
        for(int i=0; i<LayerNo; i++){
            char hname[50];
            sprintf(hname,"layer_%d",i);
            htemp[i] = new TH1F(hname,hname,100,12,30);
        }

        cout<<j<<" : "<<dname<<" : "<<_totalEntries<<endl;
        int totalHitLayer[LayerNo]={0};
        for(int entry=0; entry!=_totalEntries; ++entry)
        {
            //if(entry%10000==0) 
            //cout<<" Event : "<<entry<<endl;
            fNtuple1->GetEntry(entry);

            vector< int > hitLayer;
            for(unsigned hit=0; hit!=_newCell->size(); ++hit)
            {
                if (_newCell->at(hit) < 0) continue;
                int _layerID = _newCell->at(hit)/1e5;
                if(_layerID>29) continue;

                vector<int>::iterator result = find(hitLayer.begin(),hitLayer.end(),_layerID);
                if(result==hitLayer.end()) {
                    totalHitLayer[_layerID]++;
                    hitLayer.push_back(_layerID);
                    htemp[_layerID]->Fill(_temp->at(hit));
                }
            }
        }
        for(int i=0; i!=LayerNo; ++i)
        {
            double eff = (double)totalHitLayer[i]/(double)_totalEntries;
            // if(i!=1&&i!=29) eff += blindChn[i]*0.5/100.;
            gEff[i]->SetPoint(j,j,eff);
            if(j==nfile-1) {
                hEff0->SetBinContent(i+1,eff);
                if(i!=1&&i!=29) effmean+=eff;
            }
        }
        for(int i=0; i<LayerNo; i++){
            double temp = htemp[i]->GetMean();
            // double temp = htemp[i]->GetMean()/20.;
            gTemp[i]->SetPoint(j,j,temp);
        }
    }
    // ==============================================================
    TCanvas* c = new TCanvas("c","",10,10,800,600);
//    TLine* line = new TLine(0,0.935,30,0.935);
    TLine* line = new TLine(0,effmean/28,30,effmean/28);
    line->SetLineWidth(2);
    line->SetLineColor(2);
    hEff0->GetYaxis()->SetTitle("Efficiency");
    hEff0->GetXaxis()->SetTitle("LayerID");
    hEff0->GetXaxis()->CenterTitle();
    hEff0->GetYaxis()->CenterTitle();
    hEff0->GetYaxis()->SetRangeUser(0,1.1);
    hEff0->SetLineColor(4);
    hEff0->SetLineWidth(2);
    hEff0->SetFillStyle(3004);
    hEff0->SetFillColor(4);
    hEff0->Draw("HIST");
    line->Draw("same");
    TLatex* latex = new TLatex(.4,.7,"Preliminary");
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.06);
    latex->SetTextAlign(33);
    latex->SetTextColor(15);
    latex->Draw();
   
    // ==============================================================
    TCanvas* cc = new TCanvas("cc","",10,10,1800,800);
    gStyle->SetOptStat(false);
    TH2D* hTop = new TH2D("hTop","",10,0,nbin,10,0.5,1.01);
    hTop->GetYaxis()->SetTitle("Efficiency");
    hTop->GetXaxis()->SetTitle("Date");
    hTop->GetXaxis()->CenterTitle();
    hTop->GetYaxis()->CenterTitle();
    hTop->GetXaxis()->SetLabelSize(0);
    hTop->GetYaxis()->SetTitleSize(0.06);
    hTop->GetYaxis()->SetTitleOffset(.41);
    cc->cd();
    TPad* pad1 = new TPad("pad1","pad1",.01,.55,1,1);
    pad1->SetMargin(.07,.1,.03,.1);
    pad1->Draw();
    pad1->cd();
    hTop->Draw();
    TLegend* lg = new TLegend(.9,.03,.98,.9);
    char lname[50];
    for(int i=0; i<LayerNo; i++)
    {
        if(i<=3 || i >=28){
            gEff[i]->SetMarkerStyle(21);
            gEff[i]->SetMarkerColor(1+i);
            gEff[i]->SetLineColor(1+i);
            gEff[i]->SetMarkerSize(1);
            gEff[i]->SetLineWidth(1);
            gEff[i]->SetFillColor(false);
            gEff[i]->Draw("sameLP");
            sprintf(lname,"Layer_%d",i);
            lg->AddEntry(gEff[i],lname);
        }
    }
    lg->Draw("same");
    // -----------------------------------------------------------
    TH2D* hBottom = new TH2D("hBottom","",9,0,nbin,10,0.5,1.01);
    hBottom->GetYaxis()->SetTitle("Efficiency");
    hBottom->GetXaxis()->SetTitle("Date");
    hBottom->GetXaxis()->CenterTitle();
    hBottom->GetYaxis()->CenterTitle();
    hBottom->GetXaxis()->SetLabelSize(0.06);
    hBottom->GetXaxis()->SetTitleSize(0.06);
    hBottom->GetYaxis()->SetTitleSize(0.06);
    hBottom->GetXaxis()->SetTitleOffset(1.2);
    hBottom->GetYaxis()->SetTitleOffset(.41);
    for(int bin=0; bin<9; bin++){
        char hname[50];
        sprintf(hname,"%d",(bin+1)*10);
        hBottom->GetXaxis()->SetBinLabel(bin+1,hname);
    }
    cc->cd();
    TPad* pad2 = new TPad("pad2","pad2",.01,.05,1,.55);
    pad2->SetMargin(.07,.1,.15,.05);
    pad2->Draw();
    pad2->cd();
    hBottom->Draw();
    TLegend* lg1 = new TLegend(.9,.15,.98,.95);
    for(int i=0; i<LayerNo; i++)
    {
        if(i>3 &&i <28){
            gEff[i]->SetMarkerStyle(20);
            gEff[i]->SetMarkerColor(1+i);
            gEff[i]->SetLineColor(1+i);
            gEff[i]->SetMarkerSize(1);
            gEff[i]->SetLineWidth(1);
            gEff[i]->SetFillColor(false);
            gEff[i]->Draw("sameLP");
            sprintf(lname,"Layer_%d",i);
            lg1->AddEntry(gEff[i],lname);
        }
    }
    lg1->Draw("same");
    
    // -----------------------------------------------------------
    TCanvas* ctemp = new TCanvas("ctemp","",10,10,1800,600);
    int _MaxTemp=30;
    int _MinTemp=12;
    TH2D* h = new TH2D("h","",10,0,nbin,10,_MinTemp,_MaxTemp);
    // TH2D* h = new TH2D("h","",10,0,nbin,10,0.6,1.5);
    h->Draw();
    h->GetYaxis()->SetTitle("Temperature [ {}^{o}C ]");
    h->GetYaxis()->SetTitleOffset(.55);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitle("Date");
    h->GetXaxis()->SetTitleOffset(1.15);
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
   
    TLegend* lg2 = new TLegend(.9,.1,.97,.9);
    for(int i=0; i<LayerNo; i++){
        gTemp[i]->SetMarkerColor(1+i);
        gTemp[i]->SetMarkerSize(1.5);
        gTemp[i]->SetMarkerStyle(20);
        gTemp[i]->SetFillColor(0);
        gTemp[i]->Draw("sameP");
        char lname[50];
        sprintf(lname,"Layer_%d",i);
        lg2->AddEntry(gTemp[i],lname);
    }
    lg2->SetFillColor(0);
    lg2->SetLineWidth(0);
    lg2->Draw();

    // -----------------------------------------------------------
    TCanvas* cet = new TCanvas("cet","",10,10,1800,600);
    gPad->SetGridy();
    const int layer = 14;
    const int layer2 = 15;
    TH2D* hTempEff = new TH2D("hTempEff","",9,0,nbin,10,0.75,1.3);
    hTempEff->GetYaxis()->SetTitle("Efficiency");
    hTempEff->GetXaxis()->SetTitle("Date");
    hTempEff->GetXaxis()->CenterTitle();
    hTempEff->GetYaxis()->CenterTitle();
    hTempEff->GetXaxis()->SetLabelSize(0.04);
    hTempEff->GetXaxis()->SetTitleSize(0.04);
    hTempEff->GetYaxis()->SetLabelSize(0.04);
    hTempEff->GetYaxis()->SetTitleSize(0.04);
    hTempEff->GetXaxis()->SetTitleOffset(1.15);
    hTempEff->GetYaxis()->SetTitleOffset(.6);
    for(int bin=0; bin<9; bin++){
        char hname[50];
        sprintf(hname,"%d",(bin+1)*10);
        hTempEff->GetXaxis()->SetBinLabel(bin+1,hname);
    }
    hTempEff->Draw();

    //TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
    //                          gPad->GetUxmax(),gPad->GetUymax(),_MinTemp,_MaxTemp,510,"+L");
    TGaxis* axis = new TGaxis(nbin,.75,
                              nbin,1.3,_MinTemp,_MaxTemp,510,"+L");
    axis->SetLineColor(2);
    axis->SetLabelColor(2);
    axis->SetTitle("Temperature [ {}^{o}C ]");
    axis->SetTitleColor(2);
    axis->SetTitleOffset(.6);
    axis->SetLabelSize(.04);
    axis->CenterTitle();
    axis->Draw();

    char hname[50];
    sprintf(hname,"Layer_%d",layer);
    hTempEff->SetTitle(hname);
    TLegend* lg3 = new TLegend(.6,.75,.85,.9);
    gEff[layer]->SetMarkerColor(3);
    gEff[layer]->SetMarkerSize(1.5);
    gEff[layer]->SetMarkerStyle(20);
    gEff[layer]->SetLineWidth(2);
    gEff[layer]->SetLineColor(3);
    gEff[layer]->Draw("sameLP");
    gEff[layer2]->SetMarkerColor(4);
    gEff[layer2]->SetMarkerSize(1.5);
    gEff[layer2]->SetMarkerStyle(21);
    gEff[layer2]->SetLineWidth(2);
    gEff[layer2]->SetLineColor(4);
    gEff[layer2]->Draw("sameLP");
    gTemp[layer]->SetMarkerColor(2);
    gTemp[layer]->SetMarkerSize(1.5);
    gTemp[layer]->SetMarkerStyle(22);
    gTemp[layer]->SetLineWidth(2);
    gTemp[layer]->SetLineColor(2);
    gTemp[layer]->Draw("sameLP");
    gTemp[layer2]->SetMarkerColor(5);
    gTemp[layer2]->SetMarkerSize(1.5);
    gTemp[layer2]->SetMarkerStyle(23);
    gTemp[layer2]->SetLineWidth(2);
    gTemp[layer2]->SetLineColor(5);
    gTemp[layer2]->Draw("sameLP");
    sprintf(hname,"Layer_%d Detection Efficiency",layer);
    lg3->AddEntry(gEff[layer],hname);
    sprintf(hname,"Layer_%d Detection Efficiency",layer2);
    lg3->AddEntry(gEff[layer2],hname);
    sprintf(hname,"Layer_%d Average Temperature",layer);
    lg3->AddEntry(gTemp[layer],hname);
    sprintf(hname,"Layer_%d Average Temperature",layer2);
    lg3->AddEntry(gTemp[layer2],hname);
    lg3->Draw();
}
