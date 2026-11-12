#include "Common.hh"
#include "../adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TBox.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
namespace {
using namespace de;
using Key=std::tuple<int,int,int,int,int>;
struct Point{double x,p,lo,hi,n;};
struct Plotter{
 int worker,workers,task=0;std::string out;std::ofstream manifest;
 Plotter(int w,int ws,std::string o):worker(w),workers(ws),out(std::move(o)),manifest(table(out+"/plots_worker"+std::to_string(w)+".tsv","path\tbytes")){}
 bool take(){return task++%workers==worker;}
 void save(TCanvas&c,const std::string&path){auto dest=out+"/"+path;fs::create_directories(fs::path(dest).parent_path());fresh(dest);c.Update();fastPng(c,dest.c_str());require(fs::exists(dest)&&fs::file_size(dest)>100,"PNG creation failed: "+dest);manifest<<path<<'\t'<<fs::file_size(dest)<<'\n';}
};
void label(const std::string&s){TLatex text;text.SetNDC();text.SetTextSize(.026);text.DrawLatex(.14,.94,s.c_str());}
TGraphAsymmErrors graph(const std::vector<Point>&pts,int color,bool denom=false){TGraphAsymmErrors g;for(auto&p:pts){if(!std::isfinite(p.p)&&!denom)continue;int n=g.GetN();g.SetPoint(n,p.x,denom?p.n:p.p);g.SetPointError(n,0,0,denom?0:std::max(0.,p.p-p.lo),denom?0:std::max(0.,p.hi-p.p));}g.SetLineColor(color);g.SetMarkerColor(color);g.SetMarkerStyle(color==kBlack?20:24);g.SetMarkerSize(.65);return g;}
void comparison(Plotter&plot,const std::string&path,const std::string&title,const std::string&axis,const std::array<std::vector<Point>,2>&p,double low,double high){
 if(!plot.take())return;
 TCanvas c("compare","",1500,600);c.Divide(3,1);std::array<TGraphAsymmErrors,2>ge={graph(p[0],kBlack),graph(p[1],kRed+1)},gn={graph(p[0],kBlack,true),graph(p[1],kRed+1,true)};
 c.cd(1);gPad->SetLeftMargin(.14);auto*frame=gPad->DrawFrame(low,0,high,1.06,(title+";"+axis+";Conditional acquisition fraction").c_str());(void)frame;for(auto&g:ge)g.Draw("P SAME");TLegend legend(.15,.16,.6,.29);legend.AddEntry(&ge[0],"data","p");legend.AddEntry(&ge[1],"MC","p");legend.Draw();
 std::map<double,Point>mc;for(auto&a:p[1])mc[a.x]=a;std::vector<Point>d;for(auto&a:p[0]){auto b=mc.find(a.x);if(b==mc.end())continue;auto delta=difference({a.p,a.lo,a.hi},{b->second.p,b->second.lo,b->second.hi});d.push_back({a.x,delta.p,delta.lo,delta.hi,0});}
 auto gd=graph(d,kBlue+1);c.cd(2);gPad->SetLeftMargin(.14);gPad->DrawFrame(low,-1,high,1,(";"+axis+";data - MC (confidence interval)").c_str());gd.Draw("P SAME");
 double maxN=1;for(auto&v:p)for(auto&a:v)maxN=std::max(maxN,a.n);c.cd(3);gPad->SetLeftMargin(.16);gPad->DrawFrame(low,0,high,maxN*1.12,(";"+axis+";Absolute predicted crossings").c_str());for(auto&g:gn)g.Draw("P SAME");plot.save(c,path);
}
}
int main(int argc,char**argv)try{
 using namespace de;require(argc==6,"Draw statistics.root configuration.tsv output_dir worker workers");gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetPalette(kViridis);gStyle->SetNumberContours(100);Config cfg(argv[2]);int worker=std::stoi(argv[4]),workers=std::stoi(argv[5]);require(workers>=1&&worker>=0&&worker<workers,"Invalid plot worker");Plotter plot(worker,workers,argv[3]);auto f=openFile(argv[1]);auto*savedCfg=dynamic_cast<TNamed*>(f->Get("configuration"));require(savedCfg&&std::string(savedCfg->GetTitle())==slurp(argv[2]),"Statistics configuration mismatch");
 auto*t=tree(*f,"Summary");Result r;r.bindTo(t);std::map<Key,Result>rows;std::set<int>samples,runs;for(Long64_t i=0;i<t->GetEntries();++i){read(t,i);rows[{r.sample,r.variant,r.matched,r.level,r.id}]=r;samples.insert(r.sample);if(r.level==3)runs.insert(r.id);}t->ResetBranchAddresses();auto get=[&](int s,int v,int m,int level,int id){auto it=rows.find({s,v,m,level,id});return it==rows.end()?Result{}:it->second;};
 auto points=[&](int v,int m,int level,int start,int end){std::array<std::vector<Point>,2>p;for(int s=0;s<2;++s)for(int id=start;id<end;++id){auto a=get(s,v,m,level,id);p[s].push_back({double(id-start),a.eff,a.lo,a.hi,double(a.predicted)});}return p;};
 const std::string energy=cfg.text("energy","100")+"GeV";auto vs=variants(cfg);Geometry geometry;
 for(int matched=0;matched<2;++matched){std::string base="comparison/"+condition(0,matched,cfg)+"/"+energy;
 comparison(plot,base+"/layer_efficiency.png","Target layer excluded","Layer",points(0,matched,2,0,30),-.5,29.5);
 comparison(plot,base+"/chip_efficiency.png","Conditional chip aggregation","6*layer + chip",points(0,matched,1,0,180),-.5,179.5);
 std::array<std::vector<Point>,2>rp;for(int s=0;s<2;++s)for(int run:runs){auto a=get(s,0,matched,3,run);rp[s].push_back({double(run),a.eff,a.lo,a.hi,double(a.predicted)});}if(!runs.empty())comparison(plot,base+"/run_efficiency.png","Event-cluster intervals; MC source-run labels","Run",rp,*runs.begin()-1,*runs.rbegin()+1);
 std::array<std::vector<Point>,2>vp;for(int s=0;s<2;++s)for(size_t v=0;v<vs.size();++v){auto a=get(s,v,matched,4,0);vp[s].push_back({double(v),a.eff,a.lo,a.hi,double(a.predicted)});}comparison(plot,base+"/stability.png","0 nominal; 1 loose margin; 2 tight margin; 3 tight track","Variant",vp,-.5,vs.size()-.5);
 for(int l=0;l<30;++l){for(int chip=0;chip<6;++chip){int first=l*100000+chip*10000;auto p=points(0,matched,0,first,first+(chip==5?30:36));bool any=false;for(auto&v:p)for(auto&a:v)any|=a.n>0;if(any)comparison(plot,base+"/layer"+std::to_string(l)+"/chip"+std::to_string(chip)+"/channel_efficiency.png","Physical strip acquisition","Channel",p,-.5,chip==5?29.5:35.5);}
 for(int counts=0;counts<2;++counts){if(!plot.take())continue;double maximum=counts?1:1;for(auto&[id,strip]:geometry.strips)if(strip.layer==l)for(int s=0;s<2;++s)maximum=counts?std::max(maximum,double(get(s,0,matched,0,id).predicted)):1;
 TCanvas c("maps","",1800,600);c.SetCanvasSize(1800,600);c.Divide(3,1,.001,.001);std::vector<std::unique_ptr<TH2D>> frames;std::vector<std::unique_ptr<TBox>> boxes;std::vector<std::unique_ptr<TGaxis>>palettes;
 for(int panel=0;panel<3;++panel){c.cd(panel+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.16);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.15);std::string title="Layer "+std::to_string(l)+" "+(panel==0?"data":panel==1?"MC":"data - MC");auto h=std::make_unique<TH2D>(("frame"+std::to_string(panel)).c_str(),";x [mm];y [mm]",1,-120,120,1,-120,120);h->SetDirectory(nullptr);double low=panel==2?-maximum:0,high=maximum;h->SetMinimum(low);h->SetMaximum(high);h->Draw("AXIS");gPad->Update();
 for(auto&[id,strip]:geometry.strips)if(strip.layer==l){auto a=get(0,0,matched,0,id),b=get(1,0,matched,0,id);double value=NaN;if(counts){value=panel==0?a.predicted:panel==1?b.predicted:a.predicted-b.predicted;}else{value=panel==0?a.eff:panel==1?b.eff:std::isfinite(a.eff)&&std::isfinite(b.eff)?a.eff-b.eff:NaN;}
 double dx=l%2?2.5:22.5,dy=l%2?22.5:2.5;auto box=std::make_unique<TBox>(strip.p[0]-dx,strip.p[1]-dy,strip.p[0]+dx,strip.p[1]+dy);int color=kGray;if(std::isfinite(value)){int index=std::clamp(int((value-low)/(high-low)*(gStyle->GetNumberOfColors()-1)),0,gStyle->GetNumberOfColors()-1);color=gStyle->GetColorPalette(index);}box->SetFillColor(color);box->SetLineColor(kGray+2);box->SetLineWidth(1);box->Draw();boxes.push_back(std::move(box));}
 for(int i=0;i<100;++i){auto bar=std::make_unique<TBox>(124,-120+2.4*i,132,-120+2.4*(i+1));bar->SetFillColor(gStyle->GetColorPalette(i*(gStyle->GetNumberOfColors()-1)/99));bar->SetLineColor(bar->GetFillColor());bar->Draw();boxes.push_back(std::move(bar));} auto palette=std::make_unique<TGaxis>(132,-120,132,120,low,high,505,"+L");palette->SetLabelSize(.028);palette->SetLabelFont(42);palette->Draw();palettes.push_back(std::move(palette));frames.push_back(std::move(h));TLatex heading;heading.SetNDC();heading.SetTextSize(.045);heading.DrawLatex(.14,.94,title.c_str());heading.SetTextSize(.023);heading.DrawLatex(.14,.89,counts?"Absolute denominator; common data/MC scale":"Acquisition fraction; gray = unmeasured/masked");}
 plot.save(c,base+"/layer"+std::to_string(l)+(counts?"/denominator_map.png":"/efficiency_map.png"));
 }}
 }
 auto*dt=tree(*f,"Dependencies");int ds,dv,dm,axis,db;Long64_t dn,dk;double dp,dl,dh;
#define B(n) bind(dt,#n,&n);
 B(ds) B(dv) B(dm) B(axis) B(db) B(dn) B(dk) B(dp) B(dl) B(dh)
#undef B
 std::map<std::pair<int,int>,std::array<std::vector<Point>,2>>deps;std::array<double,8>widths={5,.5,.01,.01,.5,10,10,1};std::array<const char*,8>names={"longitudinal_mm","transverse_mm","angle_x_rad","angle_y_rad","track_RMS_mm","predicted_x_mm","predicted_y_mm","min_layers"};
 for(Long64_t i=0;i<dt->GetEntries();++i){read(dt,i);if(dv!=0)continue;deps[{dm,axis}][ds].push_back({(db+.5)*widths[axis],dp,dl,dh,double(dn)});}dt->ResetBranchAddresses();for(auto&[key,p]:deps){auto[m,ax]=key;double lo=1e20,hi=-1e20;for(auto&v:p)for(auto&a:v){lo=std::min(lo,a.x-widths[ax]);hi=std::max(hi,a.x+widths[ax]);}comparison(plot,"comparison/"+condition(0,m,cfg)+"/"+energy+"/"+names[ax]+".png","Predicted geometry; no target ADC selection",names[ax],p,lo,hi);}
 for(int s:samples)for(size_t m=0;m<groups(cfg).size();++m){if(!plot.take())continue;std::string prefix="events_s"+std::to_string(s)+"_run-1_v0_m"+std::to_string(m);auto*obs=dynamic_cast<TH1D*>(f->Get((prefix+"_observed").c_str()));auto*sup=dynamic_cast<TH1D*>(f->Get((prefix+"_supported").c_str()));auto*exp=dynamic_cast<TH1D*>(f->Get((prefix+"_expected").c_str()));auto*pred=dynamic_cast<TH1D*>(f->Get((prefix+"_predicted").c_str()));require(obs&&sup&&exp&&pred,"Missing event histograms");TCanvas c("event","",1500,650);c.Divide(2,1);c.cd(1);gPad->SetLeftMargin(.14);pred->SetLineColor(kBlack);pred->Draw("HIST");label(m?"All target-group layers excluded from tracking":"Each target excluded separately: selection bias remains");c.cd(2);gPad->SetLeftMargin(.14);obs->SetLineColor(kBlack);obs->SetMaximum(std::max(1.,obs->GetMaximum()*1.22));obs->Draw("HIST");sup->SetLineColor(kBlue+1);sup->SetLineStyle(2);sup->Draw("HIST SAME");exp->SetLineColor(kRed+1);exp->Draw("HIST SAME");TLegend legend(.36,.62,.88,.85);legend.AddEntry(obs,"Observed, >=2 predictions","l");legend.AddEntry(sup,"Observed, supported marginals","l");legend.AddEntry(exp,"Independent expectation, same events","l");legend.Draw();label("Absolute event counts; expectation uncertainty excludes rate estimation");plot.save(c,std::string(s?"mc":"data")+"/unadjusted_nominal/"+energy+"/missing_mode"+std::to_string(m)+".png");}
 std::cout<<"plot worker="<<worker<<" tasks="<<plot.task<<std::endl;return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
