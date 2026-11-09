#include "FitModel.hh"
int main(int argc,char**argv)try{
 using namespace mip;
 require(argc==6,"Usage: Fit histograms.root configuration.tsv shard.root worker workers");
 mmr::Config cfg(argv[2]);int worker=std::stoi(argv[4]),workers=std::stoi(argv[5]);
 require(workers>0&&worker>=0&&worker<workers,"Invalid worker index");fresh(argv[3]);
 auto input=mmr::openFile(argv[1]);auto meta=readMeta(*input);auto conditions=readConditions(*input);
 auto selected=selection(cfg.text("channels","all"));int mode=cfg.text("temperature-mode","current")=="current"?1:0;
 TFile output(argv[3],"CREATE");require(!output.IsZombie(),"Cannot create shard");
 Result primary,variant;
 TTree main("Diagnostics","Candidate values retained regardless of compatibility decision"),
       variants("FitVariants","1 plateau2,2 plateau3,3 plateau4,4 rebin2,5 upper85,6 seed70,7 seed130,8 threshold+sigma,9 threshold-sigma,width150,10 quadrature160,11 temperature counterpart,12 full ADC diagnostic,13 lower floor,14 higher floor,15 Gaussian-free diagnostic");
 primary.branches(main);variant.branches(variants);
 std::ofstream table(std::string(argv[3])+".tsv"),scans(std::string(argv[3])+".variants.tsv");resultHeader(table);resultHeader(scans);
 int index=0;
 for(auto&[id,m]:meta){
  if(index++%workers!=worker)continue;
  primary=Result{};primary.CellID=id;
  if(!selected.count(id)||!m.selected){primary.reason="not_selected";main.Fill();resultRow(table,primary);continue;}
  if(!m.ped_present){primary.reason="missing_pedestal";main.Fill();resultRow(table,primary);continue;}
  if(!(std::isfinite(m.threshold)&&std::isfinite(m.threshold_sigma)&&m.threshold_sigma>0)){primary.reason="missing_threshold";main.Fill();resultRow(table,primary);continue;}
  std::vector<InputCondition>histograms,other;
  for(auto&c:conditions[id]){auto h=getHist(*input,channelName(id)+"/conditions/"+c.name());(c.mode==mode?histograms:other).push_back({c,std::move(h)});}
  FitSpec nominal;nominal.rawLow=cfg.number("fit-min-adc",40);
  auto fitted=fitBest(histograms,m,cfg,nominal);primary=fitted.r;
  auto*d=directory(output,channelName(id));saveFit(*d,fitted,"nominal");
  std::vector<FitSpec>specs={{1,1,80,2},{2,1,80,3},{3,1,80,4},{4,2,80,-1},{5,1,80,-1,.85},{6,1,80,-1,1,.7},{7,1,80,-1,1,1.3},{8,1,80,-1,1,1,1},{9,1,80,-1,1,1,-1,1.5},{10,1,160,-1}};
  for(auto&spec:specs)spec.rawLow=nominal.rawLow;
  FitSpec low=nominal,high=nominal;low.variant=13;low.rawLow=std::max(0.,nominal.rawLow-20);high.variant=14;high.rawLow+=20;
  specs.push_back(low);specs.push_back(high);
  bool stable=primary.reason=="candidate";double spread=0;
  if(primary.chnEntries>=cfg.integer("min-fit-hits",200)){
   for(auto&spec:specs){
    if(spec.variant!=6&&spec.variant!=7)spec.initial={primary.LandauWidth,primary.LandauMPV,primary.GauSigma};
    auto alt=fit(histograms,m,cfg,spec);
    // A converged nominal fit supplies an efficient seed. Retry alternatives
    // with independent starts when they could change a candidate's acceptance.
    if(primary.reason=="candidate"&&(alt.r.reason!="candidate"||!std::isfinite(alt.r.LandauMPV)||std::abs(alt.r.LandauMPV/primary.LandauMPV-1)>cfg.number("max-stability",.05)))alt=fitBest(histograms,m,cfg,spec);
    variant=alt.r;variants.Fill();resultRow(scans,variant);saveFit(*d,alt,"variant"+std::to_string(spec.variant));
    if(variant.reason!="candidate")stable=false;
    if(std::isfinite(variant.LandauMPV)&&primary.LandauMPV>0)spread=std::max(spread,std::abs(variant.LandauMPV/primary.LandauMPV-1));else stable=false;
   }
   // These comparisons answer different model/data questions and are not
   // substituted for the three-free-parameter nominal measurement.
   for(int diagnostic:{11,12,15}){
    FitSpec spec=nominal;spec.variant=diagnostic;
    if(diagnostic==12)spec.rawLow=0;
    if(diagnostic==15)spec.noGaussian=true;
    auto alt=fit(diagnostic==11?other:histograms,m,cfg,spec);variant=alt.r;variants.Fill();resultRow(scans,variant);
    saveFit(*d,alt,diagnostic==11?"temperature_counterpart":"variant"+std::to_string(diagnostic));
   }
  }
  primary.stability=spread;
  if(primary.reason=="parameter_boundary"){
   if(primary.boundary_mask&2)primary.reason="MPV_at_lower_bound";
   else if(primary.boundary_mask&1)primary.reason="landau_width_unresolved";
   else if(primary.boundary_mask&4)primary.reason="gaussian_width_unresolved";
  }
  if(primary.reason=="candidate"){
   if(primary.LandauMPV-2*primary.LandauMPV_error<=std::max(primary.plateau,primary.fit_low)+primary.bin_width)primary.reason="MPV_not_observed_above_turnon";
   else if(m.threshold_fallback)primary.reason="threshold_chip_fallback";
   else if(!stable||spread>cfg.number("max-stability",.05))primary.reason="fit_condition_dependence";
   else {
    // Profile well beyond the +/-20% acceptance interval, including the
    // positive lower limit, to expose secondary/flat truncated-spectrum fits.
    std::vector<double>xs={.1,primary.LandauMPV*.25,primary.LandauMPV*.5,primary.LandauMPV*.8,primary.LandauMPV,primary.LandauMPV*1.2,primary.LandauMPV*1.5,primary.LandauMPV*2},ys;
    bool identified=true,lowerMinimum=false;
    for(double x:xs){
     if(x==primary.LandauMPV){ys.push_back(0);continue;}
     FitSpec spec=nominal;spec.hesse=false;spec.initial={primary.LandauWidth,primary.LandauMPV,primary.GauSigma};
     auto p=fitBest(histograms,m,cfg,spec,x);double delta=p.r.objective-primary.objective;ys.push_back(delta);
     if(p.r.status!=0||!std::isfinite(delta)||delta<3.84)identified=false;
     if(p.r.status==0&&std::isfinite(delta)&&delta<-.05)lowerMinimum=true;
    }
    d->cd();TGraph profile(xs.size(),xs.data(),ys.data());profile.SetName("MPV_profile");profile.Write();
    if(lowerMinimum)primary.reason="profile_lower_minimum";
    else if(!identified)primary.reason="profile_not_identifiable";
    else {primary.reason="accepted_direct";primary.decision=1;}
   }
  }
  main.Fill();resultRow(table,primary);
  std::cout<<"cell="<<id<<" n="<<primary.chnEntries<<" MPV="<<primary.LandauMPV<<" decision="<<primary.reason<<std::endl;
 }
 output.cd();main.Write();variants.Write();
 TNamed("fit_definition","Simultaneous binned Poisson shape likelihood; independent normalization per run/factor; data-only multi-start; low ADC component excluded with floor stability checks; nominal MPV and both widths free; Gaussian-free/full-ADC comparisons diagnostic only. Threshold transformed by the same factor as ADC. ChiSquare=pooled Pearson chi2; NDF=groups-normalizations-free parameters. Boundary mask bits 0,1,2=lower width,MPV,sigma; bits 3,4,5=upper. Model 0=Landau-Gaussian,1=Landau only.").Write();
 output.Close();require(bool(table)&&bool(scans),"Cannot write fit tables");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
