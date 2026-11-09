// Layer covariance after held-out position/run regression. Diagnostic only.
#define DIGITIZATION_IMPACT_REVIEW_NO_MAIN
#include "DigitizationImpactReview.cc"
#include <TH2D.h>
namespace layer_review {
using Layers=std::array<double,30>;
struct Row{Record key;Layers layer{};};
using Cohort=std::vector<Row>;
struct Table{
 std::ifstream f;std::map<std::string,size_t> cols;
 explicit Table(const fs::path&p):f(p){need(bool(f),"Missing "+p.string());std::string l;need(bool(std::getline(f,l)),"Empty table");auto h=fields(l);for(size_t i=0;i<h.size();++i)cols[h[i]]=i;}
 bool next(std::vector<std::string>&r){std::string l;if(!std::getline(f,l))return false;r=fields(l);need(r.size()==cols.size(),"Malformed table");return true;}
 const std::string&at(const std::vector<std::string>&r,const std::string&k)const{return r.at(cols.at(k));}
 double num(const std::vector<std::string>&r,const std::string&k)const{return std::stod(at(r,k));}
};
double total(const Layers&x){return std::accumulate(x.begin(),x.end(),0.);}
Values sums(const std::vector<Layers>&v){Values a;for(auto&r:v)a.push_back(total(r));return a;}
std::array<Cohort,2> readLayers(const fs::path&p){
 auto keys=records(p);std::array<std::map<std::pair<int,long long>,Record>,2>lookup;
 for(int s=0;s<2;++s)for(auto&r:keys[s])if(selected(r,2))lookup[s][{r.file,r.entry}]=r;
 Table t(p);std::vector<std::string>a;std::array<Cohort,2>out;
 while(t.next(a)){int s=t.at(a,"sample")=="mc";auto it=lookup[s].find({int(t.num(a,"file")),(long long)t.num(a,"entry")});if(it==lookup[s].end())continue;Row r;r.key=it->second;for(int l=0;l<30;++l)r.layer[l]=t.num(a,"layer"+std::to_string(l)+"_GeV");need(std::abs(total(r.layer)-r.key.energy)<1e-8,"Layer sum mismatch");out[s].push_back(r);}return out;
}
Cohort readMC(const Cohort&original,const fs::path&manifest,int beam,const std::set<int>&bad,const std::string&condition,std::ostream&audit){
 Table m(manifest);std::vector<std::string>a;std::map<int,std::string>paths;
 while(m.next(a))if(int(m.num(a,"energy_GeV"))==beam&&(condition.empty()||m.at(a,"case")==condition))paths[int(m.num(a,"file_index"))]=m.at(a,"new_calib");
 need(paths.size()==2,"Expected two MC files");Cohort out;
 for(auto&[fi,path]:paths){auto f=checked(path);auto*t=object<TTree>(*f,"Calib_Hit");int time=0,event=0;std::vector<int>*id=nullptr;std::vector<double>*e=nullptr;t->SetBranchStatus("*",false);bind(t,"Event_Time",&time);bind(t,"Event_Num",&event);bind(t,"CellID",&id);bind(t,"Hit_Energy",&e);int count=0;
  for(auto&old:original)if(old.key.file==fi&&old.key.entry<t->GetEntries()){need(t->GetEntry(old.key.entry)>0&&id&&e&&id->size()==e->size(),"Cannot read MC");need(time==old.key.time&&event==old.key.event,"MC identity mismatch");Row r;r.key=old.key;r.key.energy=0;r.key.hits=0;for(size_t j=0;j<id->size();++j){int p=physical(id->at(j));double v=e->at(j)/1000;need(std::isfinite(v),"Nonfinite hit");if(!connected(p)||bad.count(p)||v<0)continue;r.layer[p/100000]+=v;r.key.energy+=v;++r.key.hits;}need(std::abs(total(r.layer)-r.key.energy)<1e-10,"MC layer closure");out.push_back(r);++count;}
  audit<<beam<<'\t'<<(condition.empty()?"CT_input":condition)<<'\t'<<fi<<'\t'<<count<<'\t'<<path<<"\t1\n";t->ResetBranchAddresses();delete id;delete e;
 }return out;
}
Event regressionEvent(const Row&r){return Event{r.key.file,r.key.run,r.key.entry,r.key.x,r.key.y,r.key.energy,0,r.key.hits,0,0,0,0,0,0};}
std::vector<Layers> residuals(const Cohort&rows,int mode,bool data,std::ostream&models,int beam,const std::string&name){
 std::vector<Layers>out;for(auto&r:rows)out.push_back(r.layer);if(mode<0)return out;
 std::vector<int>groups;for(auto&r:rows)groups.push_back(data?r.key.run:r.key.file);std::sort(groups.begin(),groups.end());groups.erase(std::unique(groups.begin(),groups.end()),groups.end());
 auto feat=[&](const Row&r){auto e=regressionEvent(r);e.run=data?r.key.run:r.key.file;return features(e,mode,groups);};int p=feat(rows.front()).size();
 for(int held=0;held<5;++held){int n=0;Layers avg{};for(auto&r:rows)if(fold(regressionEvent(r))!=held){++n;for(int l=0;l<30;++l)avg[l]+=r.layer[l];}need(n>20*p,"Insufficient train events");for(double&v:avg)v/=n;
  TMatrixD design(n,p);int j=0;for(auto&r:rows)if(fold(regressionEvent(r))!=held){auto f=feat(r);for(int k=0;k<p;++k)design(j,k)=f[k];++j;}TDecompSVD solver(design);
  for(int l=0;l<30;++l){TVectorD y(n);j=0;for(auto&r:rows)if(fold(regressionEvent(r))!=held)y[j++]=r.layer[l];Bool_t ok;auto beta=solver.Solve(y,ok);need(ok,"Layer regression failed");for(int k=0;k<p;++k)models<<beam<<'\t'<<name<<'\t'<<mode<<'\t'<<held<<'\t'<<l<<'\t'<<k<<'\t'<<beta[k]<<'\n';for(size_t i=0;i<rows.size();++i)if(fold(regressionEvent(rows[i]))==held){auto f=feat(rows[i]);double pred=0;for(int k=0;k<p;++k)pred+=f[k]*beta[k];out[i][l]=rows[i].layer[l]-pred+avg[l];}}
 }
 for(int l=0;l<30;++l){double delta=0;for(size_t i=0;i<rows.size();++i)delta+=rows[i].layer[l]-out[i][l];delta/=rows.size();for(auto&r:out)r[l]+=delta;}return out;
}
struct Covariance{std::array<Layers,30>c{};double mean=0,diagonal=0,near=0,far=0,variance=0,frontBack=0;};
Covariance covariance(const std::vector<Layers>&rows){
 Covariance out;Layers avg{};for(auto&r:rows)for(int l=0;l<30;++l)avg[l]+=r[l]/rows.size();out.mean=total(avg);need(out.mean>0,"Invalid mean");
 for(auto&r:rows)for(int l=0;l<30;++l)for(int m=0;m<30;++m)out.c[l][m]+=(r[l]-avg[l])*(r[m]-avg[m])/(rows.size()*out.mean*out.mean);
 double a=0,b=0,ab=0;for(int l=0;l<30;++l)for(int m=0;m<30;++m){double v=out.c[l][m];out.variance+=v;if(l==m)out.diagonal+=v;else if(std::abs(l-m)<5)out.near+=v;else out.far+=v;if(l<10&&m<10)a+=v;if(l>=20&&m>=20)b+=v;if(l<10&&m>=20)ab+=v;}out.frontBack=ab/std::sqrt(a*b);need(std::abs(out.variance-std::pow(rms(sums(rows)),2))<1e-10,"Covariance sum does not close");return out;
}
}
int main(int argc,char**argv)try{
 using namespace layer_review;need(argc==3,"Usage: LayerCovarianceReview RESULT_BASE OUTPUT_DIR");fs::path base=fs::absolute(argv[1]),out=fs::absolute(argv[2]);need(!fs::exists(out/"COMPLETE"),"Completed output exists");gROOT->SetBatch();gStyle->SetOptStat(0);TH1::AddDirectory(false);
 auto bad=mask(base/"digitization/inputs/reviewer_source");std::ofstream summary(out/"summary.tsv"),impact(out/"impact.tsv"),matrix(out/"covariance.tsv"),models(out/"models.tsv"),audit(out/"validation.tsv");for(auto*f:{&summary,&impact,&matrix,&models,&audit})*f<<std::setprecision(14);
 summary<<"energy_GeV\tsample\tcondition\tevents\tmean_GeV\tRMS_over_mean\tsigma68_over_median\tdiagonal\toffdiagonal_near\toffdiagonal_far\tfront_back_correlation\n";
 impact<<"energy_GeV\tmc_sample\tcondition\tvariance_gap\tgap_reduction\tq16\tq84\tdiagonal_gap\tnear_gap\tfar_gap\toffdiagonal_share\trequired_independent_common_scale_sigma\n";
 matrix<<"energy_GeV\tsample\tcondition\tlayer_i\tlayer_j\tnormalized_covariance\n";models<<"energy_GeV\tsample\tmode\tfold\tlayer\tcoefficient_index\tcoefficient\n";audit<<"energy_GeV\tsample\tfile\tevents\tpath\tidentity_equal\n";
 TFile root((out/"layer_covariance.root").c_str(),"RECREATE");need(!root.IsZombie(),"Cannot create ROOT");const std::array<std::string,4>names={"baseline","run_file_CV","position_CV","position_run_file_CV"};
 for(int beam:{5,20,100}){
  auto integrated=readLayers(base/"integrated/comparison/20mm"/(std::to_string(beam)+"GeV/events.tsv"));std::map<std::string,Cohort>cohorts;cohorts["data"]=integrated[0];cohorts["mc_threshold_full"]=integrated[1];cohorts["mc_CT_input"]=readMC(integrated[1],base/"digitization/inputs/manifest.tsv",beam,bad,"",audit);
  if(beam==100)cohorts["mc_joint"]=readMC(integrated[1],base/"joint_response/digitization/inputs/manifest.tsv",beam,bad,"norm_pixels19600_threshold10",audit);
  std::map<std::string,std::array<std::vector<Layers>,4>>res;std::map<std::string,std::array<Covariance,4>>cov;
  for(auto&[sample,rows]:cohorts)for(int q=0;q<4;++q){auto&r=res[sample][q];r=residuals(rows,q-1,sample=="data",models,beam,sample);auto&c=cov[sample][q];c=covariance(r);auto e=sums(r);summary<<beam<<'\t'<<sample<<'\t'<<names[q]<<'\t'<<rows.size()<<'\t'<<c.mean<<'\t'<<std::sqrt(c.variance)<<'\t'<<(quantile(e,.84)-quantile(e,.16))/(2*quantile(e,.5))<<'\t'<<c.diagonal<<'\t'<<c.near<<'\t'<<c.far<<'\t'<<c.frontBack<<'\n';
   fs::path dest=out/(sample=="data"?"data":"mc")/sample/names[q]/"20mm"/(std::to_string(beam)+"GeV");fs::create_directories(dest);std::ofstream ev(dest/"events.tsv");ev<<std::setprecision(14)<<"file\tentry\trun\toriginal_GeV\tresidual_sum_GeV\n";for(size_t i=0;i<rows.size();++i)ev<<rows[i].key.file<<'\t'<<rows[i].key.entry<<'\t'<<rows[i].key.run<<'\t'<<rows[i].key.energy<<'\t'<<e[i]<<'\n';for(int l=0;l<30;++l)for(int m=0;m<30;++m)matrix<<beam<<'\t'<<sample<<'\t'<<names[q]<<'\t'<<l<<'\t'<<m<<'\t'<<c.c[l][m]<<'\n';
  }
  for(auto&[sample,rows]:cohorts)if(sample!="data"){
   double baseGap=cov["data"][0].variance-cov[sample][0].variance;std::vector<Record>dr,mr;for(auto&r:cohorts["data"])dr.push_back(r.key);for(auto&r:rows)mr.push_back(r.key);auto db=recordBlocks(dr),mb=recordBlocks(mr);std::mt19937_64 rng(2026101300+beam);std::array<Values,4>boots,ed,em;for(int q=0;q<4;++q){ed[q]=sums(res["data"][q]);em[q]=sums(res[sample][q]);}
   for(int rep=0;rep<300;++rep){auto di=resample(db,rng),mi=resample(mb,rng);std::array<double,4>gap{};for(int q=0;q<4;++q){Values a,b;for(auto i:di)a.push_back(ed[q][i]);for(auto i:mi)b.push_back(em[q][i]);gap[q]=std::pow(rms(a),2)-std::pow(rms(b),2);}if(gap[0]>0)for(int q=0;q<4;++q)boots[q].push_back(1-gap[q]/gap[0]);}
   for(int q=0;q<4;++q){auto&d=cov["data"][q];auto&m=cov[sample][q];double gap=d.variance-m.variance;impact<<beam<<'\t'<<sample<<'\t'<<names[q]<<'\t'<<gap<<'\t'<<1-gap/baseGap<<'\t'<<quantile(boots[q],.16)<<'\t'<<quantile(boots[q],.84)<<'\t'<<d.diagonal-m.diagonal<<'\t'<<d.near-m.near<<'\t'<<d.far-m.far<<'\t'<<(d.near+d.far-m.near-m.far)/gap<<'\t'<<std::sqrt(std::max(0.,gap)/(1+m.variance))<<'\n';}
  }
  for(const std::string mc:{"mc_CT_input","mc_joint"})if(cohorts.count(mc)){
   TCanvas canvas("layer_covariance","",2250,1500);canvas.SetCanvasSize(2250,1500);canvas.Divide(3,2);std::vector<std::unique_ptr<TH2D>>h;double range=0;for(int q:{0,3})for(int l=0;l<30;++l)for(int m=0;m<30;++m)range=std::max(range,std::abs(cov["data"][q].c[l][m]-cov[mc][q].c[l][m]));
   for(int row=0;row<2;++row)for(int col=0;col<3;++col){int q=row?3:0;canvas.cd(row*3+col+1);gPad->SetRightMargin(.17);auto name=std::to_string(beam)+"_"+mc+"_"+std::to_string(row)+"_"+std::to_string(col);auto plot=std::make_unique<TH2D>(name.c_str(),"",30,-.5,29.5,30,-.5,29.5);std::string title=(col==0?"Data":col==1?mc:"Data - MC")+std::string(row?", position + run/file residuals":", baseline");plot->SetTitle((title+";Layer i;Layer j;Cov(E_{i}, E_{j}) / <E>^{2}").c_str());for(int l=0;l<30;++l)for(int m=0;m<30;++m)plot->SetBinContent(l+1,m+1,col==0?cov["data"][q].c[l][m]:col==1?cov[mc][q].c[l][m]:cov["data"][q].c[l][m]-cov[mc][q].c[l][m]);if(col==2){plot->SetMinimum(-range);plot->SetMaximum(range);}plot->Draw("COLZ");root.cd();plot->Write();h.push_back(std::move(plot));}
   fs::path dest=out/"comparison"/mc/"20mm"/(std::to_string(beam)+"GeV");fs::create_directories(dest);canvas.SaveAs((dest/"layer_covariance.png").c_str());canvas.Clear();
  }
  std::cout<<beam<<" GeV layer covariance completed\n";
 }
 TNamed("method","Frozen original CoG20/window and data tail veto. CT/joint MC use original CoG and selected truth slots. Layer OLS on quadratic shower CoG and run(DATA)/file(MC); five 100-entry-block folds. Additive residuals recentered to original layer means, possibly negative: NOT calibrated energies. Covariance accounting is not causal attribution. 300 paired block bootstrap with fixed fitted coefficients, omitting refit and between-file/model uncertainty. Required common scale sigma is a reverse width benchmark, not a measured noise.").Write();root.Close();need(bool(summary)&&bool(impact)&&bool(matrix)&&bool(models)&&bool(audit),"Output write failed");std::ofstream(out/"COMPLETE")<<"complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
