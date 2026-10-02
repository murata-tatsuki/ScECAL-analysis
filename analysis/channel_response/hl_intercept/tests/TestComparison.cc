#define HL_NO_MAIN
#include "../CompareHLIntercept.cc"
using namespace hl;
void fixture(const fs::path& path,int sample,bool changed=false) {
    TFile f(path.c_str(),"RECREATE");
    TNamed("selection","same selection").Write();
    TNamed("input_pairs","same data manifest").Write();
    TParameter<Long64_t>("processed_events",1000).Write();
    for(int id: {
        0,1
    }) {
        auto*d=f.mkdir(("channel_"+std::to_string(id)).c_str());
        d->cd();
        TParameter<double>("pedestal_hg",500).Write();
        TParameter<double>("pedestal_lg",500).Write();
        for(std::string name: {
            "hg_lg","hg_lg_zoom","hg_energy","lg_energy","hg_selected_energy","lg_selected_energy"
        }) {
            TH2D h(name.c_str(),"",128,-50,3000,128,-5,130);
            h.SetDirectory(nullptr);
            TProfile p((name+"_profile").c_str(),"",128,-50,3000);
            p.SetDirectory(nullptr);
            for(int bin=1;bin<=128;++bin) {
                double x=h.GetXaxis()->GetBinCenter(bin);
                if(x<600||x>2500)continue;
                int n=id==1?4:40;
                for(int i=0;i<n;++i) {
                    double y=.03*x+(sample==1?0:-5)+(i%2?.5:-.5);
                    if(changed)y+=.1;
                    h.Fill(x,y);
                    p.Fill(x,y);
                }
            }
            h.Write();
            p.Write();
        }
        f.cd();
    }
    f.Close();
}
void reference(const fs::path& path) {
    TFile f(path.c_str(),"RECREATE");
    TTree t("InterCalib","fixture");
    int id,statistics=500;
    double slope=.03,b=-5,se=.001,be=.1,chi=1,ndf=20,xmax=4000;
    t.Branch("CellID",&id,"CellID/I");
    t.Branch("Slope",&slope,"Slope/D");
    t.Branch("Intercept",&b,"Intercept/D");
    t.Branch("SlopeError",&se,"SlopeError/D");
    t.Branch("InterceptError",&be,"InterceptError/D");
    t.Branch("ChiSquare",&chi,"ChiSquare/D");
    t.Branch("NDF",&ndf,"NDF/D");
    t.Branch("XMax",&xmax,"XMax/D");
    t.Branch("Statistics",&statistics,"Statistics/I");
    for(id=0;id<3;++id)t.Fill();
    t.Write();
}
int main(int argc,char**argv) {
    try {
        require(argc==2,"Need temporary directory");
        fs::path p=argv[1];
        Options o;
        o.output=(p/"result").string();
        o.reference=(p/"reference.root").string();
        o.channels= {
            0,1,2
        };
        o.top=1;
        reference(o.reference);
        for(int k=0;k<3;++k) {
            o.inputs[k]=(p/(labels[k]+".root")).string();
            fixture(o.inputs[k],k);
        }
        o.dataBefore=(p/"data_old.root").string();
        fixture(o.dataBefore,0);
        run(o);
        for(int k=0;k<3;++k) {
            auto f=open(p/"result"/("InterCalib_"+labels[k]+"_100GeV.root"));
            auto*t=get<TTree>(f.get(),"InterCalib");
            require(t->GetEntries()==1,"Only valid channel belongs in calibration tree");
            require(t->GetListOfBranches()->GetEntries()==9,"Schema branch count");
            require(std::string(t->GetLeaf("NDF")->GetTypeName())=="Double_t","NDF schema");
            require(std::string(t->GetLeaf("Statistics")->GetTypeName())=="Int_t","Statistics schema");
            double s,b,x;
            t->SetBranchAddress("Slope",&s);
            t->SetBranchAddress("Intercept",&b);
            t->SetBranchAddress("XMax",&x);
            t->GetEntry(0);
            require(std::abs(s-.03)<1e-10&&std::abs(b-(k==1?0:-5))<1e-8,"Known linear fit closure");
            require(x==4000,"XMax must be inherited, not set to fit upper bound");
            require(get<TTree>(f.get(),"FitDiagnostics")->GetEntries()==3,"All channels in diagnostic tree");
        }
        {
            auto f=open(p/"result/comparison/100GeV/summary.root");
            auto*g=get<TGraphErrors>(f.get(),"delta_intercept_vs_input");
            require(g->GetN()==1&&std::abs(g->GetPointY(0)+5)<1e-8,"Intercept-shift closure");
        }
        {
            auto data=open(o.inputs[0]),before=open(o.inputs[1]),after=open(o.inputs[2]);
            auto a=hist<TH2D>(before->GetDirectory("channel_0"),"hg_lg"),b=hist<TH2D>(after->GetDirectory("channel_0"),"hg_lg");
            require(projectionTV(*a,*b,true)==0,"HG distribution must remain unchanged in pure intercept fixture");
            require(projectionTV(*a,*b,false)>0,"LG projection must detect intercept shift");
            auto profile=hist<TProfile>(data->GetDirectory("channel_0"),"hg_selected_energy_profile");
            auto sw=switchCheck(*profile,500,o);
            require(!sw.leftStatus&&!sw.rightStatus&&std::abs(sw.jump)<1e-8,"Continuous energy has no artificial switch jump");
        }
        bool failed=false;
        try {
            run(o);
        }
        catch(const std::exception&) {
            failed=true;
        }
        require(failed,"Existing output must fail");
        fixture(o.dataBefore,0,true);
        o.output=(p/"bad_data").string();
        failed=false;
        try {
            run(o);
        }
        catch(const std::exception&) {
            failed=true;
        }
        require(failed,"Changed data must fail");
        TProfile profile("p","",10,0,1000);
        profile.SetDirectory(nullptr);
        for(int j=1;j<=10;++j)for(int i=0;i<40;++i)profile.Fill(profile.GetBinCenter(j),.03*profile.GetBinCenter(j)-5+(i%2?.5:-.5));
        Options limits;
        std::vector<int>bins= {
            1,2,3
        };
        require(fitProfile(profile,bins,limits).status==2,"Minimum bins guard");
        limits.minBins=3;
        require(fitProfile(profile,bins,limits).status==3,"Minimum span guard");
        limits.minSpan=100;
        require(fitProfile(profile,bins,limits).status==4,"Minimum hits guard");
        std::cout<<"PASS: fit closure, branch types, XMax provenance, sparse/missing channels, intercept shift, output protection, data integrity, fit guards"<<std::endl;
        return 0;
    }
    catch(const std::exception&e) {
        std::cerr<<"TEST FAILED: "<<e.what()<<std::endl;
        return 1;
    }
}
