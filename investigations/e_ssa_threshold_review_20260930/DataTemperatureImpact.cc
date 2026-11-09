// Fixed-event/hit calibration sensitivity; no production constants are changed.
#define DIGITIZATION_IMPACT_REVIEW_NO_MAIN
#include "DigitizationImpactReview.cc"
#undef DIGITIZATION_IMPACT_REVIEW_NO_MAIN

static double responseRatio(int id, double temperature, int pdePower) {
    const bool ten = group(id);
    const double alpha = ten ? 1.6/135. : 3.5/230.;
    const double u0 = ten ? 4.5 : 4.;
    const double p0 = ten ? .10 : .15;
    const double du = temperature-20.;
    const double uRatio = 1-.06*du/u0;
    const double gRatio = 1-alpha*du;
    const double p = p0*uRatio*uRatio;
    need(uRatio>0 && gRatio>0 && p<1, "Temperature outside diagnostic model domain");
    return gRatio*(1-p0)/(1-p)*(pdePower ? uRatio : 1.);
}

int main(int argc,char**argv) try {
    need(argc==3,"Usage: DataTemperatureImpact LABEL_RESULT_DIR OUTPUT_DIR");
    fs::path base=argv[1],out=argv[2];
    need(!fs::exists(out/"summary.tsv"),"Output already exists");
    fs::create_directories(out);
    gROOT->SetBatch(true); TH1::AddDirectory(false);
    fs::path project=fs::canonical(fs::path(__FILE__)).parent_path().parent_path().parent_path();
    auto bad=mask(project);
    fs::path p=base/"digitization/inputs/reviewer_parameters";
    auto pars=currentParameters((p/"all_hl_electron2023_hlratio_v3.root").string(),
                                (p/"all_auto_muon_v4_trackfit.root").string(),
                                (p/"pedestal2023_SPS.root").string());
    const std::array<std::string,5> names={"current_calibration","temperature_factor_one",
      "gain_inverse","gain_CT_inverse_PDEconstant","gain_CT_PDElinear_inverse"};
    std::ofstream summary(out/"summary.tsv"),impact(out/"impact.tsv"),audit(out/"validation.tsv");
    summary<<std::setprecision(14)<<"energy_GeV\tcondition\tevents\tmean_GeV\tRMS_over_mean\tsigma68_over_median\tGaussian_sigma_over_mean\tGaussian_fit_error\tGaussian_fit_status\tmean_ratio_to_current\tnegative_recomputed_hits\n";
    impact<<std::setprecision(14)<<"energy_GeV\tcondition\tpaired_events\tmean_ratio\tmean_q16\tmean_q84\tRMS_delta_pp\tRMS_delta_q16_pp\tRMS_delta_q84_pp\tGaussian_delta_pp\tGaussian_delta_q16_pp\tGaussian_delta_q84_pp\n";
    audit<<std::setprecision(14)<<"energy_GeV\tselected_events\tchecked_hits\tmax_HG_LG_saved_closure_MeV\tmax_event_cache_error_GeV\toriginal_hit_selection_frozen\n";
    for(int beam:{5,20,100}) {
        auto all=records(base/"integrated/comparison/20mm"/(std::to_string(beam)+"GeV")/"events.tsv");
        std::vector<Record> selectedRecords;
        for(auto&r:all[0]) if(selected(r,2)) selectedRecords.push_back(r);
        std::map<int,std::vector<Record>> byFile;
        for(auto&r:selectedRecords) byFile[r.file].push_back(r);
        std::ifstream manifest(base/(std::to_string(beam)+"GeV")/"baseline/inputs/data.tsv");
        need(bool(manifest),"Missing data manifest");
        std::array<Values,5> values; Values hits;
        std::array<long long,5> negative{};
        long long checkedHits=0; double closure=0,cache=0;
        std::ofstream eventTable(out/(std::to_string(beam)+"GeV_events.tsv"));
        eventTable<<std::setprecision(14)<<"file\tentry\trun\tevent_time\tevent_number\tcurrent_GeV\tno_temperature_GeV\tgain_inverse_GeV\tgain_CT_inverse_PDEconstant_GeV\tgain_CT_PDElinear_inverse_GeV\thits\n";
        std::string line; int file=-1;
        for(;std::getline(manifest,line);) {
            ++file; if(!byFile.count(file))continue;
            auto f=fields(line); need(f.size()==2,"Malformed manifest");
            auto rf=checked(f[0]),cf=checked(f[1]); RawReader raw(*rf); CalReader cal(*cf);
            for(auto&r:byFile.at(file)) {
                raw.read(r.entry);cal.read(r.entry);
                need(raw.run==r.run&&cal.run==r.run&&raw.time==r.time&&cal.time==r.time&&raw.event==r.event&&cal.event==r.event,"Event identity differs from frozen selection");
                std::map<int,size_t> index;
                for(size_t k=0;k<raw.id->size();++k)if(raw.tag->at(k)==1)need(index.emplace(raw.id->at(k),k).second,"Duplicate raw CellID");
                std::array<double,5> E{}; long long n=0;
                for(size_t j=0;j<cal.id->size();++j) {
                    int id=physical(cal.id->at(j)); if(!connected(id))continue;
                    need(index.count(cal.id->at(j)),"Calibrated hit missing from raw");
                    size_t k=index.at(cal.id->at(j)); const auto&par=pars.at(id);
                    double T=cal.temp->at(j),a=group(id)?1.6/135.:3.5/230.;
                    double c=1-a*(T-20.),qh=raw.hg->at(k)-par.ph,ql=raw.lg->at(k)-par.pl;
                    double hg=.305/par.mip*c*qh,lg=.305/par.mip*(c*ql*par.gain+par.intercept);
                    bool useLG=raw.hg->at(k)>=2600;
                    closure=std::max({closure,std::abs(hg-cal.hg->at(j)),std::abs(lg-cal.lg->at(j)),std::abs(cal.e->at(j)-(useLG?lg:hg))});++checkedHits;
                    if(bad.count(id)||cal.e->at(j)<0)continue;
                    ++n; const std::array<double,5> C={c,1.,1./c,1./responseRatio(id,T,0),1./responseRatio(id,T,1)};
                    for(int mode=0;mode<5;++mode) {
                        double signal = useLG ? (mode==0 ? C[mode]*ql*par.gain+par.intercept : C[mode]*(ql*par.gain+par.intercept)) : C[mode]*qh;
                        double e=.305/par.mip*signal/1000.;
                        if(e<0)++negative[mode];
                        E[mode]+=e;
                    }
                }
                cache=std::max(cache,std::abs(E[0]-r.energy));need(n==r.hits,"Frozen hit count differs");
                need(closure<1e-6&&cache<1e-8,"Calibration/cache closure failed");
                for(int mode=0;mode<5;++mode)values[mode].push_back(E[mode]);
                hits.push_back(n);
                eventTable<<r.file<<'\t'<<r.entry<<'\t'<<r.run<<'\t'<<r.time<<'\t'<<r.event;
                for(auto e:E)eventTable<<'\t'<<e;
                eventTable<<'\t'<<n<<'\n';
            }
        }
        // Values are traversed in file/entry order, the same as selectedRecords.
        need(values[0].size()==selectedRecords.size(),"Selected event coverage differs");
        auto blocks=recordBlocks(selectedRecords);std::mt19937_64 rng(2026101300+beam);
        std::array<Metric,5> points;
        for(int m=0;m<5;++m){points[m]=metric(values[m],hits);auto&q=points[m];summary<<beam<<'\t'<<names[m]<<'\t'<<q.n<<'\t'<<q.mean<<'\t'<<q.rms<<'\t'<<q.core<<'\t'<<q.gaus.value<<'\t'<<q.gaus.error<<'\t'<<q.gaus.status<<'\t'<<q.mean/points[0].mean<<'\t'<<negative[m]<<'\n';}
        std::array<std::array<Values,3>,5> trials;
        for(int rep=0;rep<400;++rep){auto idx=resample(blocks,rng);std::array<Values,5> e;for(auto i:idx)for(int m=0;m<5;++m)e[m].push_back(values[m][i]);double baselineR=rms(e[0]),baselineMean=mean(e[0]);auto baseFit=rep<40?gaussian(e[0]):Gaussian{};for(int m=1;m<5;++m){trials[m][0].push_back(mean(e[m])/baselineMean);trials[m][1].push_back(100*(rms(e[m])-baselineR));if(rep<40){auto g=gaussian(e[m]);if(g.status==0&&baseFit.status==0)trials[m][2].push_back(100*(g.value-baseFit.value));}}}
        for(int m=1;m<5;++m){auto q=[&](int k,double x){return trials[m][k].empty()?NAN:quantile(trials[m][k],x);};impact<<beam<<'\t'<<names[m]<<'\t'<<values[m].size()<<'\t'<<points[m].mean/points[0].mean<<'\t'<<q(0,.16)<<'\t'<<q(0,.84)<<'\t'<<100*(points[m].rms-points[0].rms)<<'\t'<<q(1,.16)<<'\t'<<q(1,.84)<<'\t'<<100*(points[m].gaus.value-points[0].gaus.value)<<'\t'<<q(2,.16)<<'\t'<<q(2,.84)<<'\n';}
        audit<<beam<<'\t'<<values[0].size()<<'\t'<<checkedHits<<'\t'<<closure<<'\t'<<cache<<"\t1\n";
    }
    std::ofstream method(out/"methodology.tsv");
    method<<"quantity\tdefinition\nselection\tOriginal good nonnegative calibrated hit set and original CoG20/energy-window/data-tail-veto event set frozen;gainselectionrawHG2600 heldfixed;inversecorrectionsapplytoHG-equivalentsignalincludingLGoffset, while currentbaselinekeepsexistingplacement;newnegativevaluesremain in same-hit sum andcounted.\nmodels\tGainratio1-alphaDeltaT,alpha10=1.6/135,alpha15=3.5/230;Uratio1-.06DeltaT/U20,U20=4.5/4V;CTp=p0Uratio^2,p0=.1/.15;PDEconstantorlinearUratio sensitivities,NOTmeasuredparameterfit.\ncorrection\tCurrentC=1-alphaDeltaT;gaininverseC=1/(1-alphaDeltaT);fullresponseinverseC=1/[gainratio*(1-p0)/(1-p(T))*PDEratio]. Lowoccupancymeanmodel;notvalidsaturationcorrectionforallhits.\nuncertainty\t400paired100-entry-blockdrawsstratifiedwithinexistingdatafiles;40Gaussianfitdraws;16-84percentileintervals,notmodelsystematics.\nMC\tThisstagechangesdataonly;comparewithpairedpilotMCthermal study separately. No originaldata/MC ROOTs overwritten.\n";
    need(bool(summary)&&bool(impact)&&bool(audit)&&bool(method),"Scientific output write failed");
    std::cout<<"Completed exact data temperature calibration sensitivities for5/20/100GeV\n";
    return 0;
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
