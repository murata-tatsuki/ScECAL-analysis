#include "Common.hh"
#include "DiskReduce.hh"
#include <TDecompSVD.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TFileCacheWrite.h>
#include <chrono>
using namespace ts;
using Key=std::array<int,10>;
using Block=std::array<Long64_t,3>;
using Group=std::array<int,4>;
using Channel=std::pair<int,int>;
Block block(const Base&b,const Config&c){return {b.run,b.file,b.entry/c.integer("block-events",250)};}
Key key(const Base&b,int g,int id,const Axis&a,int stage){return {b.sample,b.particle,g,id,a.axis,a.run,a.file,a.segment,a.bin,stage};}
struct Summary {
 Moment x,temp,value;long long denominator=0,observed=0,unique=0;
 double low=-1e100,high=1e100;bool threshold=true;
 Long64_t first=LLONG_MAX,last=0;int tmin=INT_MAX,tmax=INT_MIN;
 void merge(const Summary&b){mergeMoment(x,b.x);mergeMoment(temp,b.temp);mergeMoment(value,b.value);denominator+=b.denominator;observed+=b.observed;unique+=b.unique;low=std::max(low,b.low);high=std::min(high,b.high);threshold=threshold&&b.threshold;first=std::min(first,b.first);last=std::max(last,b.last);tmin=std::min(tmin,b.tmin);tmax=std::max(tmax,b.tmax);}
};
void count(Summary&v,const Base&b,const Axis&a,double temp){v.x.add(a.x);v.temp.add(temp);v.first=std::min(v.first,b.entry);v.last=std::max(v.last,b.entry);v.tmin=std::min(v.tmin,b.time);v.tmax=std::max(v.tmax,b.time);}
template<int N>struct Histogram:Summary {
 std::array<double,N+2> bins{};double stats[4]={};
 void add(double x,double lo,double hi){value.add(x);int i=x<lo?0:!(x<hi)?N+1:1+int(N*(x-lo)/(hi-lo));++bins[i];if(i>0&&i<=N){++stats[0];++stats[1];stats[2]+=x;stats[3]+=x*x;}}
 void merge(const Histogram&b){Summary::merge(b);for(int i=0;i<N+2;++i)bins[i]+=b.bins[i];for(int i=0;i<4;++i)stats[i]+=b.stats[i];}
 void restore(TH1D&h)const{for(int i=0;i<N+2;++i)h.SetBinContent(i,bins[i]);h.PutStats(const_cast<double*>(stats));h.SetEntries(value.n);}
};
#include "HistogramStore.hh"
struct Fitted {std::string status="insufficient_statistics";double mpv=NaN,me=NaN,w=NaN,we=NaN,sigma=NaN,se=NaN,chi=NaN;int code=-999,cov=-1;long long nfit=0;};
double langau(double*x,double*p){double center=p[1]+.22278298*p[0],sum=0,step=10*p[3]/60;for(int j=0;j<60;++j){double v=x[0]-5*p[3]+(j+.5)*step;sum+=TMath::Landau(v,center,p[0],true)*TMath::Gaus(x[0],v,p[3],true);}return p[2]*step*sum;}
Fitted fit(TH1D&h,const Summary&b,const Config&c){Fitted r;double low=std::max(5.,b.low),high=std::min(2495.,b.high);if(!b.threshold){r.status="missing_threshold";return r;}if(high-low<100){r.status="no_common_plateau";return r;}
 int first=h.FindFixBin(low),last=h.FindFixBin(high),mode=first;for(int j=first;j<=last;++j){r.nfit+=std::llround(h.GetBinContent(j));if(h.GetBinContent(j)>h.GetBinContent(mode))mode=j;}if(r.nfit<c.integer("min-fit-hits",100))return r;
 double seed=h.GetBinCenter(mode),scale=std::max(40.,seed),area=r.nfit*h.GetBinWidth(1);TF1 f("langau",langau,low,high,4);f.SetParameters(.12*scale,seed,area,.08*scale);f.SetParLimits(0,.5,500);f.SetParLimits(1,low,high);f.SetParLimits(2,.01*area,100*area);f.SetParLimits(3,.5,500);
 auto result=h.Fit(&f,"QLRS0N");r.code=int(result);r.cov=result.Get()?result->CovMatrixStatus():-1;r.chi=f.GetNDF()>0?f.GetChisquare()/f.GetNDF():NaN;
 if(r.code||r.cov<2){r.status="fit_failed_or_covariance";return r;}if(!std::isfinite(r.chi)||r.chi>c.number("fit-chi2max",5)){r.status="poor_fit_quality";return r;}
 for(int j=0;j<4;++j){double lo,hi;f.GetParLimits(j,lo,hi);if(f.GetParameter(j)-lo<.001*(hi-lo)||hi-f.GetParameter(j)<.001*(hi-lo)){r.status="parameter_boundary";return r;}}
 r.mpv=f.GetParameter(1);r.me=f.GetParError(1);r.w=f.GetParameter(0);r.we=f.GetParError(0);r.sigma=f.GetParameter(3);r.se=f.GetParError(3);r.status="ok";return r;
}
bool includeTickTrends=false;
std::vector<Axis> selectedAxes(const Base&b,double temperature,const Config&c){auto a=axes(b,temperature,c);if(!includeTickTrends)a.erase(std::remove_if(a.begin(),a.end(),[](const Axis&v){return v.axis==4;}),a.end());return a;}
#include "FitReuse.hh"
#include "FitPool.hh"
struct Tables {
 BufferedOutput trend,interval,fitout;int serial=0;
 struct Comparison {long long n=0,den=0;double v=NaN,e=NaN;std::string status="unmeasured";};
 std::map<std::tuple<int,int,int,int,std::string>,std::array<Comparison,2>>comp;
 Tables(const fs::path&o):trend(o/"trends.tsv"),interval(o/"intervals.tsv"),fitout(o/"fits.tsv"){
 std::string head="sample\tparticle\tgroup\tid\taxis\trun\tfile\tsegment\tbin\tstage";
 trend<<head<<"\tmetric\tx\tn\tdenominator\tvalue\terror\tstatus\n";
 interval<<head<<"\tn\tfirst_entry\tlast_entry\tmin_tick\tmax_tick\txmean\ttemperature_mean\n";
 fitout<<head<<"\thistogram\tn_all\tn_fit\tunderflow\toverflow\tlow\thigh\tMPV\tMPV_error\tLandau_width\twidth_error\tGaussian_sigma\tsigma_error\tchi2ndf\troot_status\tcovariance\tstatus\n";
 trend<<std::setprecision(12);fitout<<std::setprecision(12);interval<<std::setprecision(12);
 }
 void prefix(std::ostream&o,const Key&k){auto[s,p,g,id,axis,run,file,seg,bin,stage]=k;o<<s<<'\t'<<particleName(p)<<'\t'<<groupName(g)<<'\t'<<id<<'\t'<<axisName(axis)<<'\t'<<run<<'\t'<<file<<'\t'<<seg<<'\t'<<bin<<'\t'<<stage;}
 void row(const Key&k,const std::string&metric,const Summary&b,long long n,long long den,double val,double err,const std::string&status){prefix(trend,k);trend<<'\t'<<metric<<'\t'<<b.x.mean()<<'\t'<<n<<'\t'<<den<<'\t'<<val<<'\t'<<err<<'\t'<<status<<'\n';if(std::get<4>(k)==5)comp[{std::get<1>(k),std::get<2>(k),std::get<3>(k),std::get<9>(k),metric}][std::get<0>(k)]={n,den,val,err,status};}
 void range(const Key&k,const Summary&b){prefix(interval,k);interval<<'\t'<<b.x.n<<'\t'<<b.first<<'\t'<<b.last<<'\t'<<b.tmin<<'\t'<<b.tmax<<'\t'<<b.x.mean()<<'\t'<<b.temp.mean()<<'\n';}
};
#include "StreamingStats.hh"
// Bind only the scalar columns used by a pass. ROOT caches are bounded as well.
void hitColumns(TTree*t,Hit&h){t->SetBranchStatus("*",0);h.io(*t,true);for(auto name:{"source_run","event","selected","sensor","factor","x","y","length"})t->SetBranchStatus(name,0);t->SetCacheSize(8*1024*1024);}
void eventColumns(TTree*t,Event&e){t->SetBranchStatus("*",0);e.io(*t,true);for(auto name:{"source_run","event","paired","tail","cog_ok","virtual_hits","duplicate","cogx","cogy"})t->SetBranchStatus(name,0);t->SetCacheSize(8*1024*1024);}
int main(int argc,char**argv)try{
 require(argc>=3&&argc<=8,"Analyze config output_dir [memory_MiB=512] [compact|legacy] [reuse_fits.tsv|-] [standard|all] [fit_workers=1]");Config c(argv[1]);fs::path out=argv[2];size_t budget=memoryBudget(argc>=4?argv[3]:"512");require(!fs::exists(out/"statistics.root"),"Existing statistics");
 std::string trendMode=argc>6?argv[6]:"standard";require(trendMode=="standard"||trendMode=="all","Trend axes must be standard or all");includeTickTrends=trendMode=="all";
 WorkDirectory work(out);auto input=openFile((out/"flat.root").string());auto*version=dynamic_cast<TNamed*>(input->Get("schema"));require(version&&std::string(version->GetTitle())==schema,"Incompatible cached schema");
 std::cout<<"Response axes: "<<trendMode<<"; histogram storage: "<<(argc>4?argv[4]:"compact")<<std::endl;
 std::cout<<"Bounded analysis: aggregate budget "<<budget/1024/1024<<" MiB; scratch "<<work.path<<std::endl;
 TFile result((out/"statistics.root").c_str(),"CREATE","",1);require(!result.IsZombie(),"Cannot create statistics.root");new TFileCacheWrite(&result,8*1024*1024); // Registered with, and owned by, TFile.
 Tables tables(out);Histograms histograms(result,argc>4?argv[4]:"compact");FitReuse reused(argc>5?argv[5]:"-");AnalysisProgress progress(out/"inputs"/"analysis_progress.tsv");
 std::map<Channel,std::array<Moment,2>>centers;
 DiskReduce<Key,Histogram<540>> mu(work.path,"muon",budget/3);
 DiskReduce<RegBlockKey,Sums> regBlocks(work.path,"regression_blocks",budget/3);
 Hit h;auto*ht=tree(*input,"Crossings");hitColumns(ht,h);Long64_t nHits=ht->GetEntries(),nEvents=0,nBins=0;
 for(Long64_t i=0;i<nHits;++i){mmr::read(ht,i);if(h.valid){centers[{h.sample,h.cell}][0].add(h.before);centers[{h.sample,h.cell}][1].add(h.after);if(h.temp_measured&&std::isfinite(h.temperature))regBlocks.get({h.sample,h.cell,h.run,h.file,h.entry/c.integer("block-events",250)}).add(h);}
  auto aa=selectedAxes(h,h.temperature,c);for(auto[g,id]:groups(h.cell))for(auto&a:aa)for(int stage=0;stage<2;++stage){auto&b=mu.get(key(h,g,id,a,stage));count(b,h,a,h.temperature);++b.denominator;b.observed+=h.observed;b.unique+=h.unique;if(h.valid){b.add(stage?h.after:h.before,-100,2600);if(!std::isfinite(h.low)||!std::isfinite(h.high))b.threshold=false;else{b.low=std::max(b.low,h.low);b.high=std::min(b.high,h.high);}}}
  if(i&&i%1000000==0)std::cout<<"Muon aggregate "<<i<<" / "<<nHits<<std::endl;
 }ht->ResetBranchAddresses();
 BufferedOutput mpvslopes(out/"mpv_slopes.tsv");mpvslopes<<"sample\tgroup\tid\tstage\tmodel\tpoints\truns\tspan_C\tslope_ADC_per_C\terror_ADC_per_C\tstatus\n";
 using FitGroup=std::tuple<int,int,int,int>;std::map<FitGroup,std::vector<std::tuple<int,double,double,double>>>points;
 auto savePoints=[&](){
 for(auto&[g,v]:points){auto[s,group,id,stage]=g;for(int model=0;model<2;++model){std::map<int,std::array<double,3>>means;std::map<int,std::pair<double,double>>ranges;std::set<int>runs;for(auto[run,x,y,e]:v){int key=model?run:0;double w=1/(e*e);means[key][0]+=w;means[key][1]+=w*x;means[key][2]+=w*y;runs.insert(run);auto r=ranges.emplace(key,std::make_pair(x,x));r.first->second.first=std::min(r.first->second.first,x);r.first->second.second=std::max(r.first->second.second,x);}double xx=0,xy=0,span=0;for(auto&[run,r]:ranges)span=std::max(span,r.second-r.first);for(auto[run,x,y,e]:v){auto&m=means[model?run:0];x-=m[1]/m[0];y-=m[2]/m[0];xx+=x*x/(e*e);xy+=x*y/(e*e);}double beta=NaN,error=NaN;std::string status="insufficient_fit_points";long long dof=static_cast<long long>(v.size())-means.size()-1;if(span<c.number("min-temperature-span",.5))status=model?"run_temperature_confounded":"insufficient_temperature_span";else if(dof>0&&v.size()>=3&&xx>0){beta=xy/xx;double chi=0;for(auto[run,x,y,e]:v){auto&m=means[model?run:0];double residual=y-m[2]/m[0]-beta*(x-m[1]/m[0]);chi+=residual*residual/(e*e);}error=std::sqrt(std::max(1.,chi/dof)/xx);status="bin_fit_errors_position_angle_not_adjusted";}mpvslopes<<s<<'\t'<<groupName(group)<<'\t'<<id<<'\t'<<stage<<'\t'<<(model?"within_run":"pooled")<<'\t'<<v.size()<<'\t'<<runs.size()<<'\t'<<span<<'\t'<<beta<<'\t'<<error<<'\t'<<status<<'\n';}}
 points.clear();};
 auto fitTotal=mu.size();progress.write("muon_aggregate_complete",nHits,nHits);progress.write("muon_fit",0,fitTotal);long long newFits=0;auto checkpoint=std::chrono::steady_clock::now();
 int workerCount=argc>7?std::stoi(argv[7]):1;auto pool=std::make_unique<FitPool>(c,workerCount);std::cout<<"Fit workers: "<<workerCount<<"; parent is sole ROOT/TSV writer"<<std::endl;
 struct Pending{Key key;Histogram<540>value;};std::vector<Pending>batch;batch.reserve(256);
 std::array<int,3>previous={-1,-1,-1};
 auto saveFit=[&](const Key&k,const Histogram<540>&b,Fitted f){
  std::array<int,3>g={k[0],k[2],k[3]};if(g!=previous){savePoints();previous=g;}++nBins;tables.range(k,b);
  std::string name="mip_"+std::to_string(tables.serial++);TH1D hist(name.c_str(),"Same-hit ADC;path corrected ADC;hits / 5 ADC",540,-100,2600);hist.SetDirectory(nullptr);b.restore(hist);name=histograms.write(hist);
  if(k[4]==3&&f.status=="ok"&&f.me>0)points[{k[0],k[2],k[3],k[9]}].emplace_back(k[5],b.temp.mean(),f.mpv,f.me);
  tables.prefix(tables.fitout,k);tables.fitout<<'\t'<<name<<'\t'<<b.value.n<<'\t'<<f.nfit<<'\t'<<hist.GetBinContent(0)<<'\t'<<hist.GetBinContent(hist.GetNbinsX()+1)<<'\t'<<b.low<<'\t'<<b.high<<'\t'<<f.mpv<<'\t'<<f.me<<'\t'<<f.w<<'\t'<<f.we<<'\t'<<f.sigma<<'\t'<<f.se<<'\t'<<f.chi<<'\t'<<f.code<<'\t'<<f.cov<<'\t'<<f.status<<'\n';
  tables.row(k,"MPV_ADC",b,b.value.n,b.denominator,f.mpv,f.me,f.status);tables.row(k,"Landau_width_ADC",b,b.value.n,b.denominator,f.w,f.we,f.status);tables.row(k,"Gaussian_sigma_ADC",b,b.value.n,b.denominator,f.sigma,f.se,f.status);auto m=b.value;tables.row(k,"mean_ADC",b,m.n,b.denominator,m.mean(),m.se(),m.n>1?"ok_iid_SEM":"insufficient_statistics");
  if(k[9]==0){double p=double(b.observed)/b.denominator,pc=(b.observed+.5)/(b.denominator+1.);tables.row(k,"hit_fraction",b,b.observed,b.denominator,p,std::sqrt(pc*(1-pc)/(b.denominator+1.)),"conditional_efficiency_Jeffreys_SE");tables.row(k,"unique_fraction",b,b.unique,b.denominator,double(b.unique)/b.denominator,NaN,"conditional_fraction");}
  auto now=std::chrono::steady_clock::now();if(nBins%100000==0||std::chrono::duration<double>(now-checkpoint).count()>=60){tables.fitout.flush();tables.trend.flush();tables.interval.flush();require(bool(tables.fitout)&&bool(tables.trend)&&bool(tables.interval),"Checkpoint TSV write failed");progress.write("muon_fit",nBins,fitTotal,newFits,reused.reused);checkpoint=now;}
 };
 auto processBatch=[&](){if(batch.empty())return;std::vector<Fitted>results(batch.size());std::vector<FitPool::Task>tasks;for(size_t i=0;i<batch.size();++i){auto&v=batch[i];TH1D h("eligibility_adc","",540,-100,2600);h.SetDirectory(nullptr);v.value.restore(h);if(reused.get(v.key,h,v.value,c,results[i]))continue;if(requiresFit(h,v.value,c))tasks.push_back({i,&v.value});else results[i]=fit(h,v.value,c);}pool->run(tasks,results);newFits+=tasks.size();for(size_t i=0;i<batch.size();++i)saveFit(batch[i].key,batch[i].value,results[i]);batch.clear();};
 mu.each([&](const Key&k,const Histogram<540>&b){batch.push_back({k,b});if(batch.size()==256)processBatch();});processBatch();pool.reset();savePoints();reused.finish();mu.clear();tables.fitout.flush();progress.write("muon_fit_complete",nBins,fitTotal,newFits,reused.reused);std::cout<<"Muon fits complete: "<<nBins<<std::endl;
 // Fixed-effect regression: block sufficient statistics -> channel/run means -> OLS -> block scores.
 DiskReduce<MeanKey,Sums> means(work.path,"regression_means",budget/3);
 regBlocks.each([&](const RegBlockKey&k,const Sums&v){means.get({int(k[0]),int(k[1]),-1}).merge(v);means.get({int(k[0]),int(k[1]),int(k[2])}).merge(v);});
 std::map<RegKey,Regression> regressions;
 using RunKey=std::array<int,5>;DiskReduce<RunKey,RunSummary>runMeans(work.path,"run_means",budget/3);
 means.each([&](const MeanKey&k,const Sums&v){int model=k[2]>=0;for(auto[g,id]:groups(k[1]))for(int stage=0;stage<2;++stage){auto&r=regressions[{k[0],g,id,stage,model}];r.add(v,stage);if(model){for(int m=0;m<2;++m)regressions[{k[0],g,id,stage,m}].runs.insert(k[2]);auto&b=runMeans.get({k[0],g,id,stage,k[2]});double center=centers.at({k[0],k[1]})[stage].mean();b.temp.n+=v.n;b.temp.s+=v.s[0];b.temp.q+=v.q[0][0];b.value.n+=v.n;b.value.s+=v.s[5+stage]-center*v.n;long double offset=v.s[5+stage]/v.n-center;b.value.q+=v.centered[5+stage][5+stage]+v.n*offset*offset;}}});
 for(auto&[k,r]:regressions)r.solve(k[4],c);
 DiskReduce<ScoreKey,Score>scores(work.path,"regression_scores",budget/3);
 auto cursor=means.cursor();MeanKey mk;Sums mv;bool hasMean=cursor.next(mk,mv);Channel activeChannel={-1,-1};std::map<int,Sums>channelMeans;
 regBlocks.each([&](const RegBlockKey&k,const Sums&b){Channel ch{int(k[0]),int(k[1])};if(ch!=activeChannel){channelMeans.clear();activeChannel=ch;while(hasMean&&Channel{mk[0],mk[1]}==ch){channelMeans[mk[2]]=mv;hasMean=cursor.next(mk,mv);}}
  for(int model=0;model<2;++model){const auto&m=channelMeans.at(model?int(k[2]):-1);long double centered[7][7];for(int j=0;j<7;++j)for(int l=0;l<7;++l)centered[j][l]=b.centered[j][l]+b.n*(b.s[j]/b.n-m.s[j]/m.n)*(b.s[l]/b.n-m.s[l]/m.n);
   for(auto[g,id]:groups(int(k[1])))for(int stage=0;stage<2;++stage){auto&r=regressions.at({int(k[0]),g,id,stage,model});if(!r.fitted)continue;long double u=0;for(int j=0;j<5;++j){long double xy=centered[j][5+stage];for(int l=0;l<5;++l)xy-=r.coef[l]*centered[j][l];u+=r.influence[j]*xy;}scores.get({k[0],g,id,stage,model,k[2],k[3],k[4]}).u+=u;}
  }
 });scores.each([&](const ScoreKey&k,const Score&v){auto&r=regressions.at({int(k[0]),int(k[1]),int(k[2]),int(k[3]),int(k[4])});++r.blocks;r.clusterQ+=v.u*v.u;});scores.clear();means.clear();regBlocks.clear();
 BufferedOutput slopes(out/"slopes.tsv"),between(out/"run_means.tsv");slopes<<"sample\tgroup\tid\tstage\tmodel\tn\truns\tblocks\ttemperature_span_C\tslope_ADC_per_C\tiid_error\tblock_cluster_error\trank\tcorr_T_long\tcorr_T_trans\tcorr_T_angle\tstatus\n";between<<"sample\tgroup\tid\tstage\trun\tn\tmean_temperature_C\tchannel_centered_mean_ADC\tSEM_ADC\tstatus\n";
 for(auto&[k,r]:regressions){r.finish();auto[s,g,id,stage,model]=k;slopes<<s<<'\t'<<groupName(g)<<'\t'<<id<<'\t'<<stage<<'\t'<<(model?"within_channel_run_position_angle_adjusted":"pooled_channel_position_angle_adjusted")<<'\t'<<r.n<<'\t'<<r.runs.size()<<'\t'<<r.blocks<<'\t'<<r.range<<'\t'<<r.beta<<'\t'<<r.se<<'\t'<<r.cluster<<'\t'<<r.rank<<'\t'<<r.tx<<'\t'<<r.ty<<'\t'<<r.ta<<'\t'<<r.status<<'\n';}regressions.clear();
 runMeans.each([&](const RunKey&k,const RunSummary&b){between<<k[0]<<'\t'<<groupName(k[1])<<'\t'<<k[2]<<'\t'<<k[3]<<'\t'<<k[4]<<'\t'<<b.temp.n<<'\t'<<b.temp.mean()<<'\t'<<b.value.mean()<<'\t'<<b.value.se()<<"\tbetween_run_population_position_confounding_possible\n";});runMeans.clear();
 progress.write("regression_complete",1,1);std::cout<<"Regression complete"<<std::endl;
 // Only one physical event is retained for residual correlations.
 DiskReduce<ResidualKey,Summary>residualBlocks(work.path,"residual_blocks",budget/3);
 DiskReduce<CorrBlockKey,Pair>corrBlocks(work.path,"correlation_blocks",budget/3);
 BufferedOutput residual(out/"residuals.tsv"),common(out/"common_components.tsv");residual<<"sample\tfile\trun\tentry\tcell\tstage\tn_channel\tchannel_mean_ADC\tresidual_ADC\tresidual_fraction\tblock\n";common<<"sample\tparticle\tstage\tcomponent\tn_events\tvariance\tstatus\n";
 std::map<std::pair<int,int>,Moment>eventCommon;std::array<std::map<int,double>,2>eventResiduals;Base eventBase;bool haveEvent=false;
 auto saveEvent=[&](){if(!haveEvent)return;auto bl=block(eventBase,c);for(int stage=0;stage<2;++stage){auto&v=eventResiduals[stage];Moment m;for(auto[id,x]:v)m.add(x);if(m.n>=2)eventCommon[{eventBase.sample,stage}].add(m.mean());for(int group:{1,2}){std::map<int,double>q;for(auto[id,x]:v)q[group==1?id/100000*6+id%100000/10000:id/100000]=x;for(auto a=q.begin();a!=q.end();++a)for(auto b=a;b!=q.end();++b)corrBlocks.get({eventBase.sample,0,stage,group,a->first,b->first,bl[0],bl[1],bl[2]}).add(a->second,b->second);}v.clear();}};
 hitColumns(ht,h);for(Long64_t i=0;i<nHits;++i){mmr::read(ht,i);if(haveEvent&&std::tie(h.file,h.entry)!=std::tie(eventBase.file,eventBase.entry)){require(std::tie(eventBase.file,eventBase.entry)<std::tie(h.file,h.entry),"Crossings must be grouped by saved file/entry");saveEvent();}eventBase=h;haveEvent=true;if(!h.valid)continue;
  for(int stage=0;stage<2;++stage){auto m=centers.at({h.sample,h.cell})[stage];double delta=(stage?h.after:h.before)-m.mean(),fraction=m.n>=c.integer("min-residual-hits",10)&&m.mean()>0?delta/m.mean():NaN;residual<<h.sample<<'\t'<<h.file<<'\t'<<h.run<<'\t'<<h.entry<<'\t'<<h.cell<<'\t'<<stage<<'\t'<<m.n<<'\t'<<m.mean()<<'\t'<<delta<<'\t'<<fraction<<'\t'<<h.entry/c.integer("block-events",250)<<'\n';if(!std::isfinite(fraction))continue;eventResiduals[stage][h.cell]=fraction;auto bl=block(h,c);auto aa=selectedAxes(h,h.temperature,c);for(auto[g,id]:groups(h.cell))for(auto&a:aa){auto k=key(h,g,id,a,stage);ResidualKey rk;std::copy(k.begin(),k.end(),rk.begin());std::copy(bl.begin(),bl.end(),rk.begin()+10);auto&b=residualBlocks.get(rk);count(b,h,a,h.temperature);b.value.add(fraction);}
  }
 }saveEvent();ht->ResetBranchAddresses();
 Key lastKey{};bool haveResidual=false;ResidualSummary sum;
 auto saveResidual=[&](){if(!haveResidual)return;auto&m=sum.sum.value;long double mean=m.s/m.n;long double variance=std::max(0.L,sum.ss-2*mean*sum.sn+mean*mean*sum.nn);double se=sum.blocks>=4?std::sqrt(variance*sum.blocks/(sum.blocks-1.L))/m.n:NaN;tables.row(lastKey,"channel_centered_fraction",sum.sum,m.n,m.n,m.mean(),se,sum.blocks>=4?"block_cluster_SE":"insufficient_independent_blocks");};
 residualBlocks.each([&](const ResidualKey&k,const Summary&b){Key key;std::copy_n(k.begin(),10,key.begin());if(haveResidual&&key!=lastKey){saveResidual();sum=ResidualSummary{};}lastKey=key;haveResidual=true;sum.add(b);});saveResidual();residualBlocks.clear();
 for(auto&[g,m]:eventCommon){common<<g.first<<"\tmuon\t"<<g.second<<"\tevent_mean_channel_centered_fraction\t"<<m.n<<'\t'<<m.var()<<"\tmixture_of_shared_response_and_track_or_dEdx\n";for(auto name:{"within_chip","within_layer"})common<<g.first<<"\tmuon\t"<<g.second<<'\t'<<name<<"\t0\tnan\tunidentifiable_one_predicted_strip_per_layer_event\n";}
 progress.write("muon_residuals_complete",nHits,nHits);std::cout<<"Muon residuals complete"<<std::endl;
 DiskReduce<Key,Histogram<400>> energyBins(work.path,"energy",budget/3);
 DiskReduce<Key,Summary> eventBins(work.path,"event_metrics",budget/3);
 std::array<Pair,2>training;std::array<long long,2>validationCount{};std::array<double,2>trainLo={1e100,1e100},trainHi={-1e100,-1e100};long long dedicated=0;
 auto isTraining=[&](const Event&e){return c.text("split","blocks")=="runs"?e.run%2==0:(e.entry/c.integer("block-events",250))%2==0;};
 Event ev;auto*et=tree(*input,"Events");eventColumns(et,ev);nEvents=et->GetEntries();
 for(Long64_t i=0;i<nEvents;++i){mmr::read(et,i);dedicated+=ev.particle==2;if(ev.particle!=1||!ev.selected)continue;
  if(std::isfinite(ev.seconds)){if(isTraining(ev)){training[ev.sample].add(ev.seconds,ev.energy);trainLo[ev.sample]=std::min(trainLo[ev.sample],ev.seconds);trainHi[ev.sample]=std::max(trainHi[ev.sample],ev.seconds);}else ++validationCount[ev.sample];}
  auto aa=selectedAxes(ev,ev.temperature,c);for(auto&a:aa){std::array<double,8>v={ev.energy,double(ev.recorded),double(ev.strips),double(ev.raw_hits),double(ev.raw_strips),ev.signed_energy,double(ev.negative),double(ev.zero)};
   for(int metric=0;metric<8;++metric){if(metric==0){auto&b=energyBins.get(key(ev,3,-1,a,0));count(b,ev,a,ev.temperature);b.add(v[0],0,c.number("energy",100)*100);}else{auto&b=eventBins.get(key(ev,3,-1,a,metric));count(b,ev,a,ev.temperature);b.value.add(v[metric]);}}
   for(int l=0;l<30;++l){auto&b=energyBins.get(key(ev,2,l,a,0));count(b,ev,a,ev.temperature);b.add(ev.layer[l],0,c.number("energy",100)*100);}
  }
  auto bl=block(ev,c);for(int a=0;a<30;++a)for(int b=a;b<30;++b)corrBlocks.get({ev.sample,1,0,2,a,b,bl[0],bl[1],bl[2]}).add(ev.layer[a],ev.layer[b]);
  if(i&&i%100000==0)std::cout<<"Electron aggregate "<<i<<" / "<<nEvents<<std::endl;
 }et->ResetBranchAddresses();
 const std::array<std::string,8>metrics={"mean_energy_MeV","recorded_hits","physical_strips","raw_recorded_hits","raw_physical_strips","signed_energy_MeV","negative_hits","zero_hits"};
 auto eventRow=[&](const Key&k,const Summary&b){tables.range(k,b);auto m=b.value;tables.row(k,metrics[k[9]],b,m.n,m.n,m.mean(),m.se(),m.n>1?"ok_iid_SEM":"insufficient_statistics");};
 energyBins.each([&](const Key&k,const Histogram<400>&b){eventRow(k,b);auto m=b.value;double resolution=m.n>=20&&m.mean()>0?m.sd()/m.mean():NaN;tables.row(k,"resolution_RMS_mean",b,m.n,m.n,resolution,m.n>=20?resolution/std::sqrt(2.*(m.n-1)):NaN,m.n>=20?"moment_resolution_approx_iid_error":"insufficient_statistics");std::string name="energy_"+std::to_string(tables.serial++);TH1D hist(name.c_str(),"Absolute visible energy;MeV;events",400,0,c.number("energy",100)*100);hist.SetDirectory(nullptr);b.restore(hist);histograms.write(hist);});energyBins.clear();histograms.close();eventBins.each(eventRow);eventBins.clear();
 // Conditional pedestal distributions need sums/counts, never the individual ADC values.
 DiskReduce<Key,Summary>pedestal(work.path,"pedestal",budget/2);auto*pt=tree(*input,"Pedestal");hitColumns(pt,h);
 for(Long64_t i=0;i<pt->GetEntries();++i){mmr::read(pt,i);auto aa=selectedAxes(h,h.temperature,c);for(auto[g,id]:groups(h.cell))for(auto&a:aa)for(int gain=0;gain<2;++gain){auto&b=pedestal.get(key(h,g,id,a,h.memory*2+gain));count(b,h,a,h.temperature);++b.denominator;if(h.valid)b.value.add(gain?h.after:h.before);}if(i&&i%10000000==0)std::cout<<"Pedestal aggregate "<<i<<" / "<<pt->GetEntries()<<std::endl;}pt->ResetBranchAddresses();
 pedestal.each([&](const Key&k,const Summary&b){tables.range(k,b);auto m=b.value;std::string gain=k[9]%2?"LG":"HG";tables.row(k,"pedestal_"+gain+"_mean",b,m.n,b.denominator,m.mean(),m.se(),m.n>1?"HitTag0_conditional_bias_possible":"unmeasured_or_insufficient");tables.row(k,"pedestal_"+gain+"_width",b,m.n,b.denominator,m.sd(),m.n>1?m.sd()/std::sqrt(2.*(m.n-1)):NaN,m.n>1?"HitTag0_conditional_bias_possible":"unmeasured_or_insufficient");});pedestal.clear();progress.write("pedestal_complete",pt->GetEntries(),pt->GetEntries());std::cout<<"Pedestal complete"<<std::endl;
 std::map<CorrKey,Correlation>correlations;
 corrBlocks.each([&](const CorrBlockKey&k,const Pair&p){CorrKey key;std::copy_n(k.begin(),6,key.begin());auto&co=correlations[key];co.total.merge(p);++co.blocks;});
 corrBlocks.each([&](const CorrBlockKey&k,const Pair&p){CorrKey key;std::copy_n(k.begin(),6,key.begin());auto&co=correlations.at(key);co.jack.add((co.total-p).corr());});corrBlocks.clear();
 BufferedOutput corr(out/"correlations.tsv");corr<<"sample\tparticle\tstage\tgroup\ti\tj\tn_pairs\tblocks\tcovariance\tcorrelation\tblock_jackknife_error\tstatus\n";
 for(auto&[k,co]:correlations){auto[s,p,stage,g,a,b]=k;double rho=co.total.n>=20?co.total.corr():NaN;auto jack=co.jack;double se=jack.n>=4?jack.sd()*(jack.n-1)/std::sqrt(jack.n):NaN;corr<<s<<'\t'<<particleName(p)<<'\t'<<stage<<'\t'<<groupName(g)<<'\t'<<a<<'\t'<<b<<'\t'<<co.total.n<<'\t'<<co.blocks<<'\t'<<co.total.cov()<<'\t'<<rho<<'\t'<<se<<'\t'<<(std::isfinite(rho)?(p?"shower_covariance_not_electronic_noise":std::isfinite(se)?"channel_centered_pairwise_diagnostic":"insufficient_blocks_for_correlated_error"):"insufficient_pairs_or_zero_variance")<<'\n';}correlations.clear();
 // Held-out validation is represented by block moments; bootstrap samples whole blocks.
 BufferedOutput impact(out/"impact.tsv"),modelout(out/"correction_model.tsv");impact<<"sample\tsplit\tn_train\tn_validation\tvalidation_blocks\tmetric\tbefore\tafter\tpaired_difference\tblock_bootstrap_error\tstatus\n";modelout<<"sample\tn_train\tintercept_MeV\tslope_MeV_per_second\ttrain_mean_time\ttrain_span_seconds\tstatus\n";
 for(int sample=0;sample<2;++sample){auto p=training[sample];std::string status="disabled";double slope=NaN,intercept=NaN,meanT=NaN,span=NaN;
  if(c.integer("diagnose-correction",0)){status=sample?"MC_has_no_measured_time":"insufficient_training_or_validation";if(p.n>=40&&validationCount[sample]>=40){span=trainHi[sample]-trainLo[sample];meanT=p.x/p.n;double xx=p.xx-p.x*p.x/p.n;if(xx>0&&span>=c.number("time-bin",60)){slope=(p.xy-p.x*p.y/p.n)/xx;intercept=p.y/p.n;status="held_out_empirical_time_model_not_independent_calibration";}}}
  modelout<<sample<<'\t'<<p.n<<'\t'<<intercept<<'\t'<<slope<<'\t'<<meanT<<'\t'<<span<<'\t'<<status<<'\n';DiskReduce<Block,Validation>blocks(work.path,"validation_"+std::to_string(sample),budget/3);bool bad=false;
  if(std::isfinite(slope)){eventColumns(et,ev);for(Long64_t i=0;i<nEvents;++i){mmr::read(et,i);if(ev.sample!=sample||ev.particle!=1||!ev.selected||!std::isfinite(ev.seconds)||isTraining(ev))continue;double predicted=intercept+slope*(ev.seconds-meanT);if(predicted<=0){bad=true;break;}auto&b=blocks.get(block(ev,c));b.before.add(ev.energy);b.after.add(ev.energy*intercept/predicted);b.hits.add(ev.recorded);}et->ResetBranchAddresses();}
  if(bad){blocks.clear();status="invalid_prediction";}Validation total;size_t nb=0;auto blockFile=work.path/("bootstrap_"+std::to_string(sample));BufferedOutput bf(blockFile,std::ios::binary);blocks.each([&](const Block&,const Validation&v){total.merge(v);++nb;bf.write(reinterpret_cast<const char*>(&v),sizeof v);});bf.close();require(bool(bf),"Bootstrap scratch write failed");blocks.clear();
  auto b=performance(total.before,total.hits),a=performance(total.after,total.hits);std::array<Moment,3>boot;
  if(nb>=4){TRandom3 rng(41013);std::ifstream f(blockFile,std::ios::binary);for(int iter=0;iter<c.integer("bootstrap",200);++iter){Validation sum;for(size_t j=0;j<nb;++j){Validation v;f.seekg(std::streamoff(rng.Integer(nb))*sizeof v);f.read(reinterpret_cast<char*>(&v),sizeof v);require(bool(f),"Bootstrap scratch read failed");sum.merge(v);}auto pp=performance(sum.before,sum.hits),q=performance(sum.after,sum.hits);for(int j=0;j<3;++j)boot[j].add(q[j]-pp[j]);}}
  for(int j=0;j<3;++j)impact<<sample<<'\t'<<c.text("split","blocks")<<'\t'<<p.n<<'\t'<<validationCount[sample]<<'\t'<<nb<<'\t'<<(j==0?"mean_energy_MeV":j==1?"resolution_RMS_mean":"recorded_hits")<<'\t'<<b[j]<<'\t'<<a[j]<<'\t'<<a[j]-b[j]<<'\t'<<boot[j].sd()<<'\t'<<status<<(std::isfinite(slope)&&nb<4?"_insufficient_blocks":"")<<'\n';
 }
 long long gainRows=0;if(c.text("gain-table","-")!="-"){std::ifstream gain(c.text("gain-table"));require(bool(gain),"Missing gain table");BufferedOutput saved(out/"gain.tsv");saved<<gain.rdbuf();saved.close();std::ifstream scan(out/"gain.tsv");std::string line;std::getline(scan,line);while(std::getline(scan,line))if(!line.empty())++gainRows;}
 BufferedOutput availability(out/"availability.tsv");availability<<"quantity\tstatus\tdenominator\tnote\nMC_physical_time\tunmeasured\t0\tEvent_Time is generated index\nMC_measured_temperature\tunmeasured\t0\tnominal20C_is_not_a_measurement\ngain_physical_time\tunmeasured\t0\tImported LED table time_bin/period are not confirmed wall time\nSSA_virtual_hits\tunmeasured\t0\tCalib_Hit contains physical strip records\ngain\t"<<(gainRows?"imported_fit_statuses_preserved":"unmeasured")<<'\t'<<gainRows<<"\tindependent_LED_system; no_gain_to_crosstalk_inference\ndedicated_pedestal\t";availability<<(dedicated?"provided":"unmeasured")<<'\t'<<dedicated<<"\tHitTag0 is conditional even in supplied pedestal runs\n";
 BufferedOutput comparison(out/"comparison.tsv");comparison<<"particle\tgroup\tid\tstage\tmetric\tdata_n\tmc_n\tdata_denominator\tmc_denominator\tdata\tmc\tdata_minus_mc\tindependent_error\tdata_status\tmc_status\n";
 for(auto&[k,v]:tables.comp){auto[p,g,id,stage,metric]=k;comparison<<particleName(p)<<'\t'<<groupName(g)<<'\t'<<id<<'\t'<<stage<<'\t'<<metric<<'\t'<<v[0].n<<'\t'<<v[1].n<<'\t'<<v[0].den<<'\t'<<v[1].den<<'\t'<<v[0].v<<'\t'<<v[1].v<<'\t'<<v[0].v-v[1].v<<'\t'<<std::hypot(v[0].e,v[1].e)<<'\t'<<v[0].status<<'\t'<<v[1].status<<'\n';}
 result.cd();TNamed("histogram_storage",(argc>4?std::string(argv[4]):std::string("compact"))=="compact"?"sparse_bins_v1; use HistogramStore.hh readHistogram":"legacy_TH1").Write();TNamed("response_trend_axes",trendMode.c_str()).Write();TNamed("schema",schema).Write();TNamed("uncertainty","Fit covariance and trend SEM assume independent hits; slopes offer block-cluster errors, pair correlation block jackknife, optional correction paired block bootstrap. No electronics attribution from shower covariance.").Write();for(auto*stream:{&tables.trend,&tables.interval,&tables.fitout,&mpvslopes,&slopes,&between,&residual,&common,&corr,&impact,&modelout,&availability,&comparison}){stream->flush();require(bool(*stream),"Numerical TSV write failed (check disk space)");}
 require(!result.TestBit(TFile::kWriteError),"ROOT write failed");result.Close();require(!result.TestBit(TFile::kWriteError),"ROOT close failed");progress.write("analyzed_complete",1,1);std::cout<<"Analyzed "<<nEvents<<" events, "<<nHits<<" predicted crossings, "<<nBins<<" response bins\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
