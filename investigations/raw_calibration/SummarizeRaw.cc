#include "StudyCommon.hh"
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(const std::string&s){std::istringstream i(s);std::string x;std::vector<std::string>v;while(std::getline(i,x,'\t'))v.push_back(x);return v;}
std::vector<Row> read(const fs::path&p){std::ifstream f(p);require(bool(f),"Cannot read "+p.string());std::string s;std::getline(f,s);auto h=split(s);std::vector<Row> out;while(std::getline(f,s)){auto v=split(s);require(v.size()==h.size(),"Bad TSV "+p.string());Row r;for(size_t i=0;i<h.size();++i)r[h[i]]=v[i];out.push_back(r);}return out;}
double v(const Row&r,const char*k){return std::stod(r.at(k));}
double median(std::vector<double>x){if(x.empty())return 0;std::sort(x.begin(),x.end());size_t n=x.size();return (x[(n-1)/2]+x[n/2])/2;}
int main(int argc,char**argv)try{
 require(argc==2,"SummarizeRaw OUTPUT_DIRECTORY");fs::path out=argv[1];std::map<int,std::array<Row,2>> channel,resid;
 for(auto&r:read(out/"channels.tsv"))channel[int(v(r,"cell"))][r.at("sample")=="mc"]=r;
 for(auto&r:read(out/"robust_residuals.tsv"))resid[int(v(r,"cell"))][r.at("sample")=="mc"]=r;
 gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
 std::array<std::array<std::unique_ptr<TH1D>,2>,2> h;
 for(int metric=0;metric<2;++metric)for(int s=0;s<2;++s)h[metric][s]=std::make_unique<TH1D>((std::to_string(metric)+(s?"mc":"data")).c_str(),metric?";RMS / expected pedestal noise;Channels":";Central 68% width / expected pedestal noise;Channels",120,0,metric?12:3);
 std::ofstream summary(out/"comparison_summary.tsv");summary<<std::setprecision(12)<<"metric\tdata\tmc\tdefinition\n";
 std::vector<int> common;for(auto&[id,a]:resid)if(!a[0].empty()&&!a[1].empty()&&v(a[0],"n")>=100&&v(a[1],"n")>=100)common.push_back(id);
 summary<<"common_channels\t"<<common.size()<<'\t'<<common.size()<<"\tpositive MIP, common good mask, raw HG 800-2200, at least 100 hits/sample/channel\n";
 std::array<double,2> totalN{},tail{},pos{},neg{},weighted{};double sumWeight=0;
 for(int id:common){double w=std::min(v(resid[id][0],"n"),v(resid[id][1],"n"));sumWeight+=w;for(int s=0;s<2;++s){const auto&r=resid[id][s];totalN[s]+=v(r,"n");tail[s]+=v(r,"tail5sigma");pos[s]+=v(r,"tail_positive");neg[s]+=v(r,"tail_negative");weighted[s]+=w*v(r,"tail5sigma")/v(r,"n");h[0][s]->Fill(v(r,"central68_sigma_LG_ADC")/v(r,"predicted_noise_LG_ADC"));h[1][s]->Fill(v(r,"RMS_LG_ADC")/v(r,"predicted_noise_LG_ADC"));}}
 for(auto [key,a]:std::vector<std::pair<std::string,std::array<double,2>>>{{"overlap_hits",totalN},{"tail5sigma_hits",tail},{"positive_tail_hits",pos},{"negative_tail_hits",neg}})summary<<key<<'\t'<<a[0]<<'\t'<<a[1]<<"\tcommon channel set\n";
 summary<<"matched_weight_tail_fraction\t"<<weighted[0]/sumWeight<<'\t'<<weighted[1]/sumWeight<<"\tchannel weight=min(data hits,MC hits); tail relative to channel/sample median\n";
 for(auto key:{"central68_sigma_LG_ADC","RMS_LG_ADC","core3sigma_RMS"}){summary<<key;for(int s=0;s<2;++s){std::vector<double>x;for(int id:common)x.push_back(v(resid[id][s],key));summary<<'\t'<<median(x);}summary<<"\tmedian across same channels, ADC\n";}
 for(auto key:{"core_ped_HG_mean","core_ped_LG_mean","core_ped_HG_sigma","core_ped_LG_sigma"}){std::vector<double>x;for(auto&[id,a]:channel){auto&r=a[0];if(v(r,"bad")==0&&v(r,"mip")>0&&v(r,"core_ped_n")>=100)x.push_back(v(r,key));}summary<<key<<'\t'<<median(x)<<"\tnan\tdata HitTag=0 with both ADC within 5 pedestal sigma; median over "<<x.size()<<" channels, not an unbiased pedestal refit\n";}
 std::ofstream early(out/"early_switch_channels.tsv");early<<"cell\tgood_positive_mip\tHL_present\tHG_switch\tMC_early_hits\tdata_hits\tMC_hits\n";
 int aCount=0,earlyCount=0;for(auto&[id,a]:channel){auto&d=a[0];auto&m=a[1];bool good=v(m,"bad")==0&&v(m,"mip")>0;if(std::abs(v(m,"calib_HL_A")/v(m,"digi_HL_A")-1)>1e-6)++aCount;if(v(m,"MC_HG_switch")<2600){++earlyCount;early<<id<<'\t'<<good<<'\t'<<v(m,"HL_present")<<'\t'<<v(m,"MC_HG_switch")<<'\t'<<v(m,"HG4000_impliedHG_lt2600")<<'\t'<<v(d,"cal_hits")<<'\t'<<v(m,"cal_hits")<<'\n';}}
 summary<<"calib_digi_HL_slope_mismatch_channels\tnan\t"<<aCount<<"\tall connected channels\n"<<"configured_HG_switch_below2600_channels\tnan\t"<<earlyCount<<"\tall connected channels; configuration alone does not prove a populated effect\n";
 TFile f((out/"comparison_summary.root").c_str(),"RECREATE");TCanvas canvas("summary","",1300,550);canvas.Divide(2,1);for(int metric=0;metric<2;++metric){canvas.cd(metric+1);gPad->SetLogy();double max=0;for(int s=0;s<2;++s)max=std::max(max,h[metric][s]->GetMaximum());TLegend legend(.6,.73,.89,.89);for(int s=0;s<2;++s){auto&p=h[metric][s];p->SetLineColor(s?kBlue+1:kBlack);p->SetLineWidth(2);p->SetMaximum(max*2);p->SetMinimum(.5);p->Write();p->Draw(s?"HIST SAME":"HIST");legend.AddEntry(p.get(),s?"MC":"data","l");}legend.DrawClone();}canvas.SaveAs((out/"figures/residual_width_comparison.png").c_str());
 // Event-level coincidence of the low-LG tails, on the fixed shower+CoG20 selection.
 auto readout=read(out/"event_readout.tsv");std::array<std::array<Stats,3>,2> energy;std::array<std::array<long long,3>,2> negativeHits{};
 TH1D multiplicityData("tail_multiplicity_data",";Low-LG outliers / event (HG 800-2200);Events",131,-.5,130.5),multiplicityMC("tail_multiplicity_mc","",131,-.5,130.5);
 TH1D quiet("data_quiet",";Positive good event energy [GeV];Probability / bin",90,3,5.7),affected("data_affected","",90,3,5.7);
 for(auto&r:readout){int sample=r.at("sample")=="mc";bool selected=v(r,"shower_cog20")>0,flag=v(r,"overlap_negative_tail")>=3;if(selected){(sample?multiplicityMC:multiplicityData).Fill(v(r,"overlap_negative_tail"));if(!sample)(flag?affected:quiet).Fill(v(r,"good_energy_GeV"));for(int k=0;k<3;++k)if(k==0||(k==1&&flag)||(k==2&&!flag)){energy[sample][k].add(v(r,"good_energy_GeV"));negativeHits[sample][k]+=v(r,"LG_selected_negative");}}}
 std::ofstream pop(out/"event_readout_summary.tsv");pop<<std::setprecision(12)<<"sample\tselection\tevents\tmean_GeV\tRMS_GeV\tRMS_over_mean\tnegative_selected_LG_hits\n";for(int sample=0;sample<2;++sample)for(int k=0;k<3;++k){auto&a=energy[sample][k];pop<<(sample?"mc":"data")<<'\t'<<(k==0?"shower_cog20":k==1?"shower_cog20_lowLG_atleast3":"shower_cog20_lowLG_below3")<<'\t'<<a.n<<'\t'<<a.mean<<'\t'<<a.rms()<<'\t'<<(a.n?a.rms()/a.mean:0)<<'\t'<<negativeHits[sample][k]<<'\n';}
 f.cd();multiplicityData.Write();multiplicityMC.Write();quiet.Write();affected.Write();
 canvas.Clear();canvas.Divide(2,1);canvas.cd(1);gPad->SetLogy();multiplicityData.SetLineColor(kBlack);multiplicityMC.SetLineColor(kBlue+1);multiplicityData.SetMinimum(.5);multiplicityData.SetMaximum(10000);multiplicityData.Draw("HIST");multiplicityMC.Draw("HIST SAME");TLegend ml(.53,.73,.89,.89);ml.AddEntry(&multiplicityData,"data","l");ml.AddEntry(&multiplicityMC,"MC","l");ml.DrawClone();
 canvas.cd(2);if(quiet.Integral()>0)quiet.Scale(1./quiet.Integral());if(affected.Integral()>0)affected.Scale(1./affected.Integral());quiet.SetLineColor(kBlack);affected.SetLineColor(kRed+1);quiet.SetMaximum(std::max(quiet.GetMaximum(),affected.GetMaximum())*1.3);quiet.Draw("HIST");affected.Draw("HIST SAME");TLegend el(.46,.72,.89,.89);el.AddEntry(&quiet,"data: <3 low-LG outliers","l");el.AddEntry(&affected,"data: >=3 low-LG outliers","l");el.DrawClone();canvas.SaveAs((out/"figures/event_low_lg_coincidence.png").c_str());
 // Independent consistency checks between the two passes.
 auto ev=read(out/"event_summary.tsv"),tailEv=read(out/"tail_event_effect.tsv");require(ev.size()==tailEv.size(),"Selection count mismatch");for(size_t i=0;i<ev.size();++i){require(ev[i].at("sample")==tailEv[i].at("sample")&&ev[i].at("selection")==tailEv[i].at("selection"),"Selection mismatch");require(v(ev[i],"events")==v(tailEv[i],"events")&&std::abs(v(ev[i],"good_energy_GeV")-v(tailEv[i],"mean_GeV"))<1e-9,"Event sum mismatch between passes");}
 for(auto&r:read(out/"validation.tsv"))require(v(r,"formula_mismatch_hits")==0,"Saved calibration does not match expected formula");
 std::cout<<"Validated event samples/totals across passes and raw-to-calibration closure; common channels="<<common.size()<<std::endl;
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
