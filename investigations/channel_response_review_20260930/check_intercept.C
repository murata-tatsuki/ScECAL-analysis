#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <iostream>
#include <iomanip>
void check_intercept(){
 const char*base="/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/";
 TFile hl((std::string(base)+"all_hl_electron2023_hlratio_v3.root").c_str());auto*t=hl.Get<TTree>("InterCalib");int id;double s,b; t->SetBranchAddress("CellID",&id);t->SetBranchAddress("Slope",&s);t->SetBranchAddress("Intercept",&b);
 double slope=0,intercept=0;for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(id==240008){slope=s;intercept=b;}}
 TFile mip((std::string(base)+"all_auto_muon_v4_trackfit.root").c_str());auto*m=mip.Get<TTree>("MIP_Fit");double mpv,chi;int ndf; m->SetBranchAddress("CellID",&id);m->SetBranchAddress("LandauMPV",&mpv);m->SetBranchAddress("ChiSquare",&chi);m->SetBranchAddress("NDF",&ndf);
 double value=0,total=0;int count=0;for(Long64_t i=0;i<m->GetEntries();++i){m->GetEntry(i);int layer=id/100000;if(layer>=30||(layer>3&&layer<28))continue;if(chi/ndf>1.7)continue;total+=mpv;++count;if(id==240008)value=mpv;}
 bool fallback=value==0;if(fallback)value=total/count;double k=.305/value,shift=-k*intercept/slope;
 std::cout<<std::setprecision(12)<<"cell=240008 slope="<<slope<<" intercept="<<intercept<<" mip="<<value<<" fallback="<<fallback<<" k="<<k<<" expected_shift_MeV="<<shift<<std::endl;
 TFile hist("/home/murata_t/ScECAL_BeamTest/investigations/result/channel_response_review_20260930/channel8_cell240008_calibrated_hg_lg.root");
 for(const char*name:{"data_p","mc_p"}){auto*p=hist.Get<TProfile>(name);TGraphErrors g;double n=0,sum=0;bool mc=std::string(name)=="mc_p";
 for(int j=1;j<=p->GetNbinsX();++j){double x=p->GetBinCenter(j);if(x<.3||x>=1.2)continue;double entries=p->GetBinEntries(j);if(entries==0)continue;double y=p->GetBinContent(j)-(mc?shift:0);int point=g.GetN();g.SetPoint(point,x,y);g.SetPointError(point,0,p->GetBinError(j));n+=entries;sum+=entries*y;}
 TF1 fit("fit","[0]+[1]*(x-0.75)",.3,1.2);g.Fit(&fit,"Q");std::cout<<name<<" corrected_mean="<<sum/n<<" n="<<n<<" residual_at_0.75="<<fit.GetParameter(0)<<" +- "<<fit.GetParError(0)<<" slope="<<fit.GetParameter(1)<<" +- "<<fit.GetParError(1)<<" chi2/ndf="<<fit.GetChisquare()<<"/"<<fit.GetNDF()<<std::endl;
 }
}
