#include "CalibrationAudit.hh"
#include "ResidualStats.hh"
#include "../adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <climits>

struct Options {
    std::string manifest,pedestal,mip,hl,out,sample,channels="all",layers="0-29";
    int eventsPerFile=0,minHits=30;
    double tempWidth=.5,adcMin=-100,adcMax=4200,adcWidth=50;
    double cut=2600,window=200,switchWidth=25,overlapMin=800,overlapMax=2200;
    bool plots=true;
};
Options options(int argc,char**argv) {
    Options o;
    for(int i=1;i<argc;i+=2) {
        if(i+1==argc)throw std::runtime_error("Missing option value");
        std::string k=argv[i],v=argv[i+1];
        if(k=="--manifest")o.manifest=v;else if(k=="--pedestal")o.pedestal=v;
        else if(k=="--mip")o.mip=v;else if(k=="--hl")o.hl=v;else if(k=="--out")o.out=v;
        else if(k=="--sample")o.sample=v;else if(k=="--channels")o.channels=v;else if(k=="--layers")o.layers=v;
        else if(k=="--events-per-file")o.eventsPerFile=std::stoi(v);else if(k=="--min-hits")o.minHits=std::stoi(v);
        else if(k=="--temp-width")o.tempWidth=std::stod(v);else if(k=="--adc-min")o.adcMin=std::stod(v);
        else if(k=="--adc-max")o.adcMax=std::stod(v);else if(k=="--adc-width")o.adcWidth=std::stod(v);
        else if(k=="--switch-adc")o.cut=std::stod(v);else if(k=="--switch-window")o.window=std::stod(v);
        else if(k=="--switch-bin-width")o.switchWidth=std::stod(v);
        else if(k=="--overlap-min")o.overlapMin=std::stod(v);else if(k=="--overlap-max")o.overlapMax=std::stod(v);
        else if(k=="--plots") { if(v!="0"&&v!="1")throw std::runtime_error("plots must be 0 or 1");o.plots=v=="1"; }
        else throw std::runtime_error("Unknown option "+k);
    }
    for(double v:{o.tempWidth,o.adcMin,o.adcMax,o.adcWidth,o.cut,o.window,o.switchWidth,o.overlapMin,o.overlapMax})
        if(!std::isfinite(v))throw std::runtime_error("Nonfinite option");
    if(o.manifest.empty()||o.pedestal.empty()||o.mip.empty()||o.hl.empty()||o.out.empty()||
       (o.sample!="data"&&o.sample!="mc")||o.eventsPerFile<0||o.minHits<3||o.tempWidth<=0||o.adcWidth<=0||
       o.adcMax<=o.adcMin||o.switchWidth<=0||o.window<3*o.switchWidth||o.overlapMax<=o.overlapMin||o.overlapMax>o.cut)
        throw std::runtime_error("Invalid options; see run.sh --help");
    if(std::abs(o.window/o.switchWidth-std::round(o.window/o.switchWidth))>1e-8)
        throw std::runtime_error("switch-window must be a multiple of switch-bin-width");
    return o;
}
struct Channel {
    std::map<int,ResidualBin> bins;
    std::map<int,SwitchBin> near;
    Regression selectedLeft,selectedRight,pairedNear,overlap;
    Moments temp,closureHG,closureLG;
    long long hits=0,adcOutside=0,nonfinite=0,selectionMismatch=0,closureInvalid=0;
    double closureHGMax=0,closureLGMax=0,leftGap=INFINITY,rightGap=INFINITY;
};
using Group=std::tuple<int,int,int>; // Run_Num, temperature bin, physical CellID; one input file at a time.
struct PngStats {
    long long saved=0,skippedEmpty=0;
    double rasterSeconds=0,encodeSeconds=0;
};
struct Tables {
    std::ofstream bins,near,summary,files;
    explicit Tables(const std::string& out):bins(out+"/residual_bins.tsv"),near(out+"/switch_bins.tsv"),summary(out+"/channels.tsv"),files(out+"/files.tsv") {
        if(!bins||!near||!summary||!files)throw std::runtime_error("Cannot create output tables");
        for(auto*p:{&bins,&near,&summary,&files})*p<<std::setprecision(17);
        const std::string key="sample\tfile_index\trun\ttemp_bin_low_C\ttemp_bin_high_C\tcellid\t";
        bins<<key<<"adc_bin_low\tadc_bin_high\tn\tmean_adc\tmean_delta_MeV\tsem_delta_MeV\tq16_MeV\tmedian_MeV\tq84_MeV\tcentral68_width_MeV\tstatus\n";
        near<<key<<"raw_hg_bin_low\traw_hg_bin_high\tn\tmean_raw_hg\thg_mean\thg_sem\thg_q16\thg_median\thg_q84\tlg_mean\tlg_sem\tlg_q16\tlg_median\tlg_q84\tselected_mean\tselected_sem\tselected_q16\tselected_median\tselected_q84\tstatus\n";
        summary<<key<<"hits\tnonfinite\tadc_outside\ttemp_mean_C\tselection_mismatch\tclosure_invalid\thg_recalc_minus_saved_mean\thg_closure_abs_max\tlg_recalc_minus_saved_mean\tlg_closure_abs_max\toverlap_n\toverlap_intercept_MeV\toverlap_slope_MeV_per_ADC\toverlap_median_bin_slope\tswitch_left_n\tswitch_right_n\tswitch_left_gap_ADC\tswitch_right_gap_ADC\tswitch_status\tselected_jump_OLS_MeV\tselected_jump_OLS_SE_MeV\tselected_jump_median_bins_MeV\tpaired_delta_at_switch_OLS_MeV\n";
        files<<"sample\tfile_index\traw\tcalibration\tmatching\traw_entries\tprocessed_events\trun_min\trun_max\tevent_time_min\tevent_time_max\tmissing_temperature_hits\n";
    }
    void check() { for(auto*p:{&bins,&near,&summary,&files}){p->flush();if(!*p)throw std::runtime_error("Output table write failed");} }
};
std::string rowKey(const Options&o,int index,int run,int tb,int id) {
    std::ostringstream s;s<<std::setprecision(17)<<o.sample<<'\t'<<index<<'\t'<<run<<'\t'<<tb*o.tempWidth<<'\t'<<(tb+1)*o.tempWidth<<'\t'<<id<<'\t';return s.str();
}
void addGraph(TGraphErrors&g,double x,const Moments&y) {
    int i=g.GetN();g.SetPoint(i,x,y.mean);g.SetPointError(i,0,y.sem());
}
void addBand(TGraphAsymmErrors&g,double x,const Distribution&y) {
    int i=g.GetN();double m=y.quantile(.5);g.SetPoint(i,x,m);g.SetPointError(i,0,0,m-y.quantile(.16),y.quantile(.84)-m);
}
void style(TGraph&g,int color,int marker) {g.SetLineColor(color);g.SetMarkerColor(color);g.SetMarkerStyle(marker);g.SetMarkerSize(.6);g.SetLineWidth(2);}

void writeChannel(const Options&o,int fi,const Group&key,Channel&ch,const CalibrationParameter&p,TFile&root,Tables&tables,PngStats&png) {
    auto[run,tb,id]=key;
    std::string prefix=rowKey(o,fi,run,tb,id);
    std::string path="file_"+std::to_string(fi)+"/run_"+std::to_string(run)+"/temperature_bin_"+std::to_string(tb)+"/channel_"+std::to_string(id);
    auto*dir=root.GetDirectory(path.c_str());if(!dir)dir=root.mkdir(path.c_str(),"",true);if(!dir)throw std::runtime_error("Cannot create ROOT directory");dir->cd();
    TGraphErrors mean,hgm,lgm,selm;TGraphAsymmErrors median,hgmed,lgmed,selmed;
    mean.SetName("delta_mean_sem");median.SetName("delta_median_q16_q84");
    hgm.SetName("switch_hg_mean_sem");lgm.SetName("switch_lg_mean_sem");selm.SetName("switch_selected_mean_sem");
    hgmed.SetName("switch_hg_median_q16_q84");lgmed.SetName("switch_lg_median_q16_q84");selmed.SetName("switch_selected_median_q16_q84");
    Regression robustOverlap,robustLeft,robustRight;
    int leftBins=0,rightBins=0;double leftBinGap=INFINITY,rightBinGap=INFINITY;
    for(auto&[b,v]:ch.bins) {
        v.delta.sort();auto&m=v.delta.m;double lo=o.adcMin+b*o.adcWidth,hi=lo+o.adcWidth;
        tables.bins<<prefix<<lo<<'\t'<<hi<<'\t'<<m.n<<'\t'<<v.x.mean<<'\t'<<m.mean<<'\t'<<m.sem()<<'\t'
            <<v.delta.quantile(.16)<<'\t'<<v.delta.quantile(.5)<<'\t'<<v.delta.quantile(.84)<<'\t'
            <<v.delta.quantile(.84)-v.delta.quantile(.16)<<'\t'<<(m.n>=o.minHits?"ok":"low_statistics")<<'\n';
        if(m.n<o.minHits)continue;
        addGraph(mean,v.x.mean,m);addBand(median,v.x.mean,v.delta);
        if(lo+p.ph>=o.overlapMin&&hi+p.ph<=o.overlapMax)robustOverlap.add(v.x.mean,v.delta.quantile(.5));
    }
    for(auto&[b,v]:ch.near) {
        v.hg.sort();v.lg.sort();v.selected.sort();double lo=o.cut+b*o.switchWidth;
        tables.near<<prefix<<lo<<'\t'<<lo+o.switchWidth<<'\t'<<v.x.n<<'\t'<<v.x.mean;
        for(auto*d:{&v.hg,&v.lg,&v.selected})tables.near<<'\t'<<d->m.mean<<'\t'<<d->m.sem()<<'\t'<<d->quantile(.16)<<'\t'<<d->quantile(.5)<<'\t'<<d->quantile(.84);
        tables.near<<'\t'<<(v.x.n>=o.minHits?"ok":"low_statistics")<<'\n';
        if(v.x.n<o.minHits)continue;
        addGraph(hgm,v.x.mean,v.hg.m);addGraph(lgm,v.x.mean,v.lg.m);addGraph(selm,v.x.mean,v.selected.m);
        addBand(hgmed,v.x.mean,v.hg);addBand(lgmed,v.x.mean,v.lg);addBand(selmed,v.x.mean,v.selected);
        if(b<0){++leftBins;leftBinGap=std::min(leftBinGap,o.cut-(lo+o.switchWidth));robustLeft.add(v.x.mean-o.cut,v.selected.quantile(.5));}
        else{++rightBins;rightBinGap=std::min(rightBinGap,lo-o.cut);robustRight.add(v.x.mean-o.cut,v.selected.quantile(.5));}
    }
    bool supported=leftBins>=3&&rightBins>=3&&leftBinGap<=o.switchWidth&&rightBinGap<=o.switchWidth&&
        ch.selectedLeft.n>=o.minHits&&ch.selectedRight.n>=o.minHits&&ch.selectedLeft.valid()&&ch.selectedRight.valid();
    std::string status=supported?"ok":"insufficient_support";
    double jump=supported?ch.selectedRight.at(0)-ch.selectedLeft.at(0):missing();
    double error=supported?std::hypot(ch.selectedLeft.errorAt(0),ch.selectedRight.errorAt(0)):missing();
    double robustJump=supported&&robustLeft.valid()&&robustRight.valid()?robustRight.at(0)-robustLeft.at(0):missing();
    tables.summary<<prefix<<ch.hits<<'\t'<<ch.nonfinite<<'\t'<<ch.adcOutside<<'\t'<<ch.temp.mean<<'\t'<<ch.selectionMismatch<<'\t'
        <<ch.closureInvalid<<'\t'<<(ch.closureHG.n?ch.closureHG.mean:missing())<<'\t'<<(ch.closureHG.n?ch.closureHGMax:missing())<<'\t'
        <<(ch.closureLG.n?ch.closureLG.mean:missing())<<'\t'<<(ch.closureLG.n?ch.closureLGMax:missing())<<'\t'
        <<ch.overlap.n<<'\t'<<ch.overlap.at(0)<<'\t'<<ch.overlap.slope()<<'\t'<<robustOverlap.slope()<<'\t'
        <<ch.selectedLeft.n<<'\t'<<ch.selectedRight.n<<'\t'<<ch.leftGap<<'\t'<<ch.rightGap<<'\t'<<status<<'\t'
        <<jump<<'\t'<<error<<'\t'<<robustJump<<'\t'<<(supported?ch.pairedNear.at(0):missing())<<'\n';
    style(mean,kBlack,20);style(median,kBlue+1,24);
    style(hgm,kBlue+1,20);style(lgm,kRed+1,24);style(selm,kBlack,21);
    style(hgmed,kBlue+1,20);style(lgmed,kRed+1,24);style(selmed,kBlack,21);
    const std::vector<TGraph*> graphs={&mean,&median,&hgm,&lgm,&selm,&hgmed,&lgmed,&selmed};
    for(TGraph*g:graphs)g->Write();
    TParameter<Long64_t>("hits",ch.hits).Write();TParameter<double>("raw_switch_adc",o.cut).Write();
    TParameter<double>("pedestal_hg",p.ph).Write();TParameter<double>("selected_jump_OLS_MeV",jump).Write();
    TParameter<double>("selected_jump_median_bins_MeV",robustJump).Write();TNamed("switch_status",status.c_str()).Write();
    // Keep ROOT/TSV statistics, but do not build a canvas when no panel has points.
    // Use all graphs so a switch-only plot is kept even if the residual ADC range is empty.
    if(!o.plots)return;
    if(std::none_of(graphs.begin(),graphs.end(),
            [](const TGraph* graph){ return graph->GetN()>0; })) {
        ++png.skippedEmpty;
        return;
    }
    fs::path figure=fs::path(o.out)/"figures"/path;fs::create_directories(figure.parent_path());
    std::ostringstream title;title<<o.sample<<", file "<<fi<<", run "<<run<<", CellID "<<id<<", "<<tb*o.tempWidth<<" <= T < "<<(tb+1)*o.tempWidth<<" C";
    TCanvas canvas("residual_canvas",title.str().c_str(),2100,1100);canvas.Divide(3,2);
    auto residualPanel=[&](int panel,double xmin,double xmax,const std::string&label){
        canvas.cd(panel);TGraphErrors gm;TGraphAsymmErrors gq;double ymin=0,ymax=0;
        for(const auto&[b,v]:ch.bins){(void)b;
            if(v.delta.m.n<o.minHits||v.x.mean<xmin||v.x.mean>=xmax)continue;
            addGraph(gm,v.x.mean,v.delta.m);addBand(gq,v.x.mean,v.delta);
            ymin=std::min({ymin,v.delta.quantile(.16),v.delta.m.mean-v.delta.m.sem()});
            ymax=std::max({ymax,v.delta.quantile(.84),v.delta.m.mean+v.delta.m.sem()});
        }
        double margin=std::max(.05,(ymax-ymin)*.15);ymin-=margin;ymax+=margin;
        gPad->DrawFrame(xmin,ymin,xmax,ymax,(label+(gm.GetN()?"":" (insufficient bins)")+";HG - pedestal [ADC];E_{LG} - E_{HG} [MeV]").c_str());
        style(gm,kBlack,20);style(gq,kBlue+1,24);
        if(gq.GetN())gq.DrawClone("P SAME");
        if(gm.GetN())gm.DrawClone("P SAME");
        TLine zero(xmin,0,xmax,0);zero.SetLineStyle(2);zero.DrawClone();
        if(o.cut-p.ph>=xmin&&o.cut-p.ph<=xmax){TLine cut(o.cut-p.ph,ymin,o.cut-p.ph,ymax);cut.SetLineStyle(3);cut.DrawClone();}
        TLegend legend(.13,.77,.61,.89);legend.AddEntry(&mean,"Mean +/- SEM","lp");legend.AddEntry(&median,"Median, q16-q84 (spread)","lp");legend.DrawClone();
    };
    residualPanel(1,o.adcMin,o.adcMax,title.str());
    residualPanel(2,o.overlapMin-p.ph,o.overlapMax-p.ph,"Residual: overlap zoom");
    residualPanel(3,o.cut-o.window-p.ph,o.cut+o.window-p.ph,"Residual: switch zoom");
    auto energyPanel=[&](int panel,TGraph&a,TGraph&b,TGraph&c,const std::string&label){
        canvas.cd(panel);
        // All branch distributions share a frame, so a shifted LG curve is never clipped by HG.
        double lo=INFINITY,hi=-INFINITY;
        for(auto&[idx,v]:ch.near){(void)idx;if(v.x.n<o.minHits)continue;
            for(auto*d:{&v.hg,&v.lg,&v.selected}){lo=std::min({lo,d->quantile(.16),d->m.mean-d->m.sem()});hi=std::max({hi,d->quantile(.84),d->m.mean+d->m.sem()});}}
        if(!std::isfinite(lo)||!std::isfinite(hi)){lo=-1;hi=1;}
        double margin=std::max(.1,(hi-lo)*.12);lo-=margin;hi+=margin;
        gPad->DrawFrame(o.cut-o.window,lo,o.cut+o.window,hi,(label+";Raw HG [ADC];Energy [MeV]").c_str());
        if(a.GetN())a.Draw("P SAME");
        if(b.GetN())b.Draw("P SAME");
        if(c.GetN())c.Draw("P SAME");
        TLine cut(o.cut,lo,o.cut,hi);cut.SetLineStyle(2);cut.DrawClone();
        TLegend legend(.12,.74,.46,.9);legend.AddEntry(&a,"Saved E_{HG}","lp");legend.AddEntry(&b,"Saved E_{LG}","lp");legend.AddEntry(&c,"Saved selected E","lp");legend.DrawClone();
    };
    energyPanel(4,hgm,lgm,selm,"Switch: mean +/- SEM");energyPanel(5,hgmed,lgmed,selmed,"Switch: median, q16-q84 (spread)");
    canvas.cd(6);TGraphErrors selectedFit;
    selectedFit.SetTitle(("Selected-energy fits at raw HG switch: "+status+";Raw HG - switch [ADC];Selected energy [MeV]").c_str());
    for(const auto&[b,v]:ch.near)if(v.x.n>=o.minHits)addGraph(selectedFit,v.x.mean-o.cut,v.selected.m);
    style(selectedFit,kBlack,20);
    if(selectedFit.GetN()) {selectedFit.Draw("AP");if(supported){
        TLine left(-o.window,ch.selectedLeft.at(-o.window),0,ch.selectedLeft.at(0));left.SetLineColor(kBlue+1);left.DrawClone();
        TLine right(0,ch.selectedRight.at(0),o.window,ch.selectedRight.at(o.window));right.SetLineColor(kRed+1);right.DrawClone();}}
    else gPad->DrawFrame(-o.window,-1,o.window,1,"Insufficient switch statistics;Raw HG - switch [ADC];Selected energy [MeV]");
    const auto timing=fastPng(canvas,(figure.string()+".png").c_str());
    ++png.saved;png.rasterSeconds+=timing.raster;png.encodeSeconds+=timing.encode;
    if(png.saved%1000==0)std::cout<<"  PNG progress: saved="<<png.saved
        <<" skipped_empty="<<png.skippedEmpty<<std::endl;
}

void run(const Options&o) {
    gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);
    // Keep axis labels readable independently of site ROOT startup styles.
    gStyle->SetTextFont(42);gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelSize(.035,"XYZ");gStyle->SetTitleSize(.04,"XYZ");
    gStyle->SetTitleOffset(1.2,"X");gStyle->SetTitleOffset(1.6,"Y");
    gStyle->SetPadLeftMargin(.17);gStyle->SetPadRightMargin(.04);
    gStyle->SetPadBottomMargin(.13);gStyle->SetPadTopMargin(.10);
    fs::create_directories(o.out);
    if(fs::exists(fs::path(o.out)/"residuals.root"))throw std::runtime_error("Output exists; use another label");
    const auto params=readCalibration(o.pedestal,o.mip,o.hl);writeCalibration(params,o.out+"/calibration_channels.tsv");
    auto layers=integerList(o.layers,true);auto cells=o.channels=="all"?std::set<int>{}:integerList(o.channels,false);
    for(int l:layers)if(l>=30)throw std::runtime_error("Only layers 0-29 supported");
    for(int id:cells)if(!params.count(id)||!layers.count(id/100000))throw std::runtime_error("Invalid selected physical CellID");
    auto selected=[&](int id){return params.count(id)&&layers.count(id/100000)&&(cells.empty()||cells.count(id));};
    Tables tables(o.out);auto root=openFile(o.out+"/residuals.root.tmp","RECREATE");
    TNamed("sample",o.sample.c_str()).Write();
    TNamed("energy_source","Saved Hit_HG_Energy, Hit_LG_Energy, Hit_Energy; no subtraction of an MC intercept").Write();
    TNamed("quantiles","Exact type-7 quantiles of untrimmed hit values in each ADC/run/temperature/file bin. q16-q84 is spread, not uncertainty. Low statistics bins remain in TSV.").Write();
    TNamed("provenance_status","Candidate constants only; see inputs/calibration_provenance.tsv and closure columns. Run/time coverage is NOT a constants validity interval.").Write();
    TNamed("fit_policy","Switch: separate hit OLS on each side evaluated at the same raw HG threshold; median-bin fit is unweighted OLS of bin medians. Formal SE assumes independent hits, excludes systematics/event correlations. No automatic pass/fail.").Write();
    std::ifstream input(o.manifest);if(!input)throw std::runtime_error("Cannot open manifest");
    std::vector<HitIndex> hitIndex(32*6*100*36);uint64_t generation=0;int fileIndex=0;
    for(std::string line;std::getline(input,line);) {
        if(line.empty())continue;
        auto sep=line.find('\t');if(sep==std::string::npos)throw std::runtime_error("Manifest requires raw TAB calib");
        auto rawFile=openFile(line.substr(0,sep));auto calFile=openFile(line.substr(sep+1));
        auto*raw=tree(*rawFile,"Raw_Hit");auto*cal=tree(*calFile,"Calib_Hit");
        int rr=0,rt=0,ri=0,cr=0,ct=0,ci=0;
        bind(raw,"Run_Num",&rr);bind(raw,"Event_Time",&rt);bind(raw,"TriggerID",&ri);
        bind(cal,"Run_Num",&cr);bind(cal,"Event_Time",&ct);bind(cal,"Event_Num",&ci);
        // Verify the complete sequence even for a bounded hit sample: CycleID is absent from calibration.
        std::vector<EventKey> keys;keys.reserve(cal->GetEntries());
        for(Long64_t e=0;e<cal->GetEntries();++e){if(cal->GetEntry(e)<=0)throw std::runtime_error("Cannot read calibration key");keys.emplace_back(cr,ct,ci);}
        bool ordered=raw->GetEntries()==cal->GetEntries();
        if(ordered)for(Long64_t e=0;e<raw->GetEntries();++e){if(raw->GetEntry(e)<=0)throw std::runtime_error("Cannot read raw key");if(EventKey{rr,rt,ri}!=keys[e]){ordered=false;break;}}
        std::map<EventKey,Long64_t> lookup;
        if(!ordered) {
            for(size_t e=0;e<keys.size();++e)if(!lookup.emplace(keys[e],e).second)throw std::runtime_error("Ambiguous duplicate calibration event key");
            std::set<EventKey> seen;
            for(Long64_t e=0;e<raw->GetEntries();++e){if(raw->GetEntry(e)<=0)throw std::runtime_error("Cannot read raw key");EventKey k{rr,rt,ri};
                if(!seen.insert(k).second||!lookup.count(k))throw std::runtime_error("Duplicate or unmatched raw event key");}
            if(seen.size()!=lookup.size())throw std::runtime_error("Unmatched calibration events");
        }
        std::vector<int>*rawIDs=nullptr,*tags=nullptr,*calIDs=nullptr;
        std::vector<double>*hg=nullptr,*eh=nullptr,*el=nullptr,*energy=nullptr,*temp=nullptr,*lg=nullptr;
        bind(raw,"CellID",&rawIDs);bind(raw,"HitTag",&tags);bind(raw,"HG_Charge",&hg);bind(raw,"LG_Charge",&lg);
        bind(cal,"CellID",&calIDs);bind(cal,"Hit_HG_Energy",&eh);bind(cal,"Hit_LG_Energy",&el);bind(cal,"Hit_Energy",&energy);bind(cal,"NewTemperature",&temp);
        cacheActiveBranches(raw);cacheActiveBranches(cal);
        std::map<Group,Channel> groups;long long processed=0,missingTemp=0;
        int minRun=INT_MAX,maxRun=INT_MIN,minTime=INT_MAX,maxTime=INT_MIN;
        Long64_t n=raw->GetEntries(),take=o.eventsPerFile?std::min<Long64_t>(o.eventsPerFile,n):n;
        std::cout<<o.sample<<" file "<<fileIndex<<": "<<rawFile->GetName()<<"; "<<take<<" / "<<n<<" events; "<<(ordered?"verified entry order":"unique event keys")<<std::endl;
        for(Long64_t j=0;j<take;++j) {
            // Systematic midpoint sampling covers the full run, not just its beginning.
            Long64_t e=take==n?j:static_cast<Long64_t>((j+.5L)*n/take);
            if(raw->GetEntry(e)<=0)throw std::runtime_error("Cannot read raw entry");
            Long64_t ce=ordered?e:lookup.at(EventKey{rr,rt,ri});
            if(cal->GetEntry(ce)<=0)throw std::runtime_error("Cannot read calibration entry");
            if(EventKey{rr,rt,ri}!=EventKey{cr,ct,ci})throw std::runtime_error("Event changed after key validation");
            if(rawIDs->size()!=tags->size()||rawIDs->size()!=hg->size()||rawIDs->size()!=lg->size()||calIDs->size()!=eh->size()||calIDs->size()!=el->size()||calIDs->size()!=energy->size()||calIDs->size()!=temp->size())throw std::runtime_error("Hit vector length mismatch");
            ++generation;
            for(size_t k=0;k<rawIDs->size();++k) {
                auto&slot=hitIndex.at(hitSlot(rawIDs->at(k)));
                if(slot.raw==generation)throw std::runtime_error("Duplicate full CellID in raw event");
                slot.raw=generation;slot.index=k;if(tags->at(k)==1)slot.signal=generation;
            }
            for(size_t k=0;k<calIDs->size();++k) {
                int full=calIDs->at(k),id=physicalID(full);auto&slot=hitIndex.at(hitSlot(full));
                if(slot.calibrated==generation)throw std::runtime_error("Duplicate full CellID in calibration event");
                slot.calibrated=generation;
                if(!selected(id))continue;
                if(slot.signal!=generation)throw std::runtime_error("Calibration hit missing raw HitTag=1 counterpart");
                double t=temp->at(k);if(!std::isfinite(t)||std::abs(t/o.tempWidth)>INT_MAX-1){++missingTemp;continue;}
                auto&ch=groups[{rr,static_cast<int>(std::floor(t/o.tempWidth)),id}];++ch.hits;ch.temp.add(t);
                size_t r=slot.index;double h=hg->at(r),l=lg->at(r),a=eh->at(k),b=el->at(k),s=energy->at(k);
                if(!std::isfinite(h)||!std::isfinite(l)||!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(s)){++ch.nonfinite;continue;}
                const auto&p=params.at(id);if(!p.pedPresent)throw std::runtime_error("Missing pedestal for selected observed channel "+std::to_string(id));
                double x=h-p.ph,d=b-a,expect=h<o.cut?a:b;
                if(std::abs(s-expect)>1e-7+1e-6*std::abs(expect))++ch.selectionMismatch;
                double correction=1-(t-20)*(sipmGroup(id)?1.6/135:3.5/230);
                double recHG=.305*(h-p.ph)*correction/p.mip,recLG=.305*((l-p.pl)*correction*p.gain+p.offset)/p.mip;
                if(std::isfinite(recHG)&&std::isfinite(recLG)){ch.closureHG.add(recHG-a);ch.closureLG.add(recLG-b);ch.closureHGMax=std::max(ch.closureHGMax,std::abs(recHG-a));ch.closureLGMax=std::max(ch.closureLGMax,std::abs(recLG-b));}else ++ch.closureInvalid;
                if(x>=o.adcMin&&x<o.adcMax)ch.bins[static_cast<int>(std::floor((x-o.adcMin)/o.adcWidth))].add(x,d);else ++ch.adcOutside;
                if(h>=o.overlapMin&&h<o.overlapMax)ch.overlap.add(x,d);
                if(h>=o.cut-o.window&&h<o.cut+o.window) {
                    ch.near[static_cast<int>(std::floor((h-o.cut)/o.switchWidth))].add(h,a,b,s);ch.pairedNear.add(h-o.cut,d);
                    if(h<o.cut){ch.selectedLeft.add(h-o.cut,s);ch.leftGap=std::min(ch.leftGap,o.cut-h);}
                    else{ch.selectedRight.add(h-o.cut,s);ch.rightGap=std::min(ch.rightGap,h-o.cut);}
                }
            }
            ++processed;minRun=std::min(minRun,rr);maxRun=std::max(maxRun,rr);minTime=std::min(minTime,rt);maxTime=std::max(maxTime,rt);
            if(processed%10000==0)std::cout<<"  processed "<<processed<<" events"<<std::endl;
        }
        tables.files<<o.sample<<'\t'<<fileIndex<<'\t'<<line<<'\t'<<(ordered?"verified_entry_order":"unique_event_key")<<'\t'<<n<<'\t'<<processed<<'\t'
            <<(processed?std::to_string(minRun):"NA")<<'\t'<<(processed?std::to_string(maxRun):"NA")<<'\t'<<(processed?std::to_string(minTime):"NA")<<'\t'<<(processed?std::to_string(maxTime):"NA")<<'\t'<<missingTemp<<'\n';
        PngStats png;
        for(auto&[key,ch]:groups)writeChannel(o,fileIndex,key,ch,params.at(std::get<2>(key)),*root,tables,png);
        if(o.plots)std::cout<<o.sample<<" file "<<fileIndex<<" PNG: saved="<<png.saved
            <<" skipped_empty="<<png.skippedEmpty<<" raster_seconds="<<png.rasterSeconds
            <<" encode_seconds="<<png.encodeSeconds<<std::endl;
        tables.check();raw->ResetBranchAddresses();cal->ResetBranchAddresses();
        delete rawIDs;delete tags;delete hg;delete lg;delete calIDs;delete eh;delete el;delete energy;delete temp;
        ++fileIndex;
    }
    if(!fileIndex)throw std::runtime_error("Empty manifest");
    tables.check();root->Write();if(root->TestBit(TFile::kWriteError))throw std::runtime_error("ROOT write failed");root->Close();
    fs::rename(o.out+"/residuals.root.tmp",o.out+"/residuals.root");
    std::cout<<"Saved "<<o.out<<std::endl;
}
int main(int argc,char**argv) {try{run(options(argc,argv));return 0;}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}}
