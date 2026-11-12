#include "Common.hh"
#include "DiskReduce.hh"
#include <sys/wait.h>
#include <cerrno>
#include "../adc_energy/FastPng.hh"
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TBox.h>
using namespace ts;
using Row=std::map<std::string,std::string>;
std::vector<Row> table(const fs::path&p){std::ifstream f(p);require(bool(f),"Missing saved table "+p.string());std::string line;std::getline(f,line);if(!line.empty()&&line.back()=='\r')line.pop_back();auto keys=split(line);std::vector<Row>v;while(std::getline(f,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();auto a=split(line);require(a.size()==keys.size(),"Bad TSV "+p.string());Row r;for(size_t j=0;j<a.size();++j)r[keys[j]]=a[j];v.push_back(std::move(r));}return v;}
template<class F>void streamTable(const fs::path&p,F fn){std::ifstream f(p);require(bool(f),"Missing saved table "+p.string());std::string line;std::getline(f,line);auto keys=split(line);while(std::getline(f,line)){auto a=split(line);require(a.size()==keys.size(),"Bad TSV "+p.string());Row r;for(size_t j=0;j<a.size();++j)r[keys[j]]=a[j];fn(r,line);}require(f.eof(),"Table read failed");}
void sortTrends(const fs::path&input,const fs::path&output,const fs::path&temp,size_t budget){
 std::string memory="--buffer-size="+std::to_string(budget/1024/1024)+"M",tmp="--temporary-directory="+temp.string(),dst="--output="+output.string();
 pid_t pid=fork();require(pid>=0,"Cannot start bounded TSV sort");if(pid==0){setenv("LC_ALL","C",1);execlp("sort","sort","--stable","--parallel=1",memory.c_str(),tmp.c_str(),dst.c_str(),"-t","\t","-k1,1n","-k2,2","-k3,3","-k4,4n","-k11,11","-k5,5","-k6,6n","-k7,7n","-k8,8n","-k9,9n","-k10,10n",input.c_str(),static_cast<char*>(nullptr));_exit(127);}
 int status=0;while(waitpid(pid,&status,0)<0)require(errno==EINTR,"Cannot wait for TSV sort");require(WIFEXITED(status)&&WEXITSTATUS(status)==0,"Bounded TSV sort failed");
}
struct PlotPoint {int run,file,segment,stage;long long n;double x,value,error;};
double val(const Row&r,const std::string&k){return std::stod(r.at(k));}
int num(const Row&r,const std::string&k){return std::stoi(r.at(k));}
int particle(const Row&r){return r.at("particle")=="muon"?0:r.at("particle")=="electron"?1:2;}
fs::path dir(const Row&r,const Config&c){auto p=subdir(num(r,"sample"),particle(r),c);if(r.at("group")=="channel"){int id=num(r,"id");p/="layer"+std::to_string(id/100000);p/="chip"+std::to_string(id%100000/10000);p/="channel"+std::to_string(id%100);}else if(r.at("group")=="chip"){int id=num(r,"id");p/="layer"+std::to_string(id/6);p/="chip"+std::to_string(id%6);}else if(r.at("group")=="layer")p/="layer"+r.at("id");return p;}
void save(TCanvas&canvas,const fs::path&p){fs::create_directories(p.parent_path());canvas.Modified();canvas.Update();fastPng(canvas,p.c_str());require(fs::file_size(p)>0,"PNG write failed");}
void note(const std::string&s){TLatex t;t.SetNDC();t.SetTextSize(.035);t.DrawLatex(.12,.5,s.c_str());}
int main(int argc,char**argv)try{
 require(argc==5||argc==6,"Draw config result_dir plot_dir channels [memory_MiB=512]");size_t budget=memoryBudget(argc==6?argv[5]:"512");Config c(argv[1]);fs::path out=argv[2],plots=argv[3];require(!fs::exists(plots),"Plot directory exists; choose a new --plot-label");fs::create_directories(plots);gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetCanvasPreferGL(false);
 std::set<int>selected;if(std::string(argv[4])!="all")for(auto&s:split(argv[4],','))selected.insert(std::stoi(s));
 WorkDirectory work(plots);auto filtered=work.path/"filtered.tsv",sorted=work.path/"sorted.tsv";BufferedOutput filter(filtered);
 streamTable(out/"trends.tsv",[&](const Row&r,const std::string&line){if(r.at("axis")=="pooled")return;
 if(!selected.empty()&&r.at("group")=="chip"&&!std::any_of(selected.begin(),selected.end(),[&](int id){return id/100000*6+id%100000/10000==num(r,"id");}))return;
 if(!selected.empty()&&r.at("group")=="layer"&&!std::any_of(selected.begin(),selected.end(),[&](int id){return id/100000==num(r,"id");}))return;
 if(r.at("group")=="channel"&&!selected.empty()&&!selected.count(num(r,"id")))return;
 if(!std::isfinite(val(r,"value"))||r.at("metric")=="unique_fraction")return;
 filter<<line<<'\n';});filter.close();require(bool(filter),"Trend scratch write failed");sortTrends(filtered,sorted,work.path,budget/2);
 BufferedOutput index(plots/"plot_index.tsv");index<<"path\tmeaning\tpoints\n";TCanvas canvas("canvas","time stability",960,720);
 auto drawTrend=[&](const Row&r,const std::vector<PlotPoint>&v){bool channel=r.at("group")=="channel";if(channel){long long total=0;for(auto&a:v)total+=a.n;if(total<c.integer("min-plot-hits",40))return;}
 canvas.Clear();canvas.SetRightMargin(.05);canvas.SetLeftMargin(.13);canvas.SetBottomMargin(.13);
 TMultiGraph graph;std::map<std::tuple<int,int,int,int>,std::unique_ptr<TGraphErrors>>lines;TLegend legend(.60,.71,.95,.92);legend.SetTextSize(.022);
 for(auto&a:v){int run=a.run,file=a.file,seg=a.segment,stage=a.stage;auto key=std::make_tuple(stage,r.at("axis")=="run"?0:run,file,seg);auto&g=lines[key];if(!g){g=std::make_unique<TGraphErrors>();int color=stage%2?kRed+1:kBlue+1;g->SetLineColor(color);g->SetMarkerColor(color);g->SetMarkerStyle(20+(std::abs(run)%5));g->SetMarkerSize(.65);if(lines.size()<=8){std::string label="run "+std::to_string(run)+" stage "+std::to_string(stage);if(r.at("metric").find("pedestal")!=std::string::npos)label="run "+std::to_string(run)+" memory "+std::to_string(stage/2);legend.AddEntry(g.get(),label.c_str(),"lp");}}
 int j=g->GetN();g->SetPoint(j,a.x,a.value);double e=a.error;g->SetPointError(j,0,std::isfinite(e)?e:0);}
 std::string x=r.at("axis");if(x=="event_order")x=num(r,"sample")?"generated event order (not physical time)":"raw entry order";if(x=="timestamp_ticks")x="saved integer timestamp (unit unconfirmed)";
 std::string title=sampleName(num(r,"sample"))+" "+r.at("particle")+" "+r.at("group")+r.at("id")+";"+x+";"+r.at("metric");graph.SetTitle(title.c_str());for(auto&[key,g]:lines)graph.Add(g.release(),"P");graph.Draw("A");legend.Draw();
 fs::path path=plots/dir(r,c)/(r.at("metric")+"_vs_"+r.at("axis")+".png");save(canvas,path);index<<fs::relative(path,plots).string()<<'\t'<<r.at("metric")<<"; absolute response; separate run/segment series\t"<<v.size()<<'\n';};
 // The external sort places one figure's points together; no full-table copies remain.
 using GK=std::tuple<int,std::string,std::string,int,std::string,std::string>;GK previous;Row metadata;std::vector<PlotPoint>points;
 std::ifstream head(out/"trends.tsv"),sortedInput(sorted);std::string header,line;std::getline(head,header);auto keys=split(header);
 while(std::getline(sortedInput,line)){auto a=split(line);require(a.size()==keys.size(),"Bad sorted trend row");Row r;for(size_t j=0;j<a.size();++j)r[keys[j]]=a[j];GK k={num(r,"sample"),r.at("particle"),r.at("group"),num(r,"id"),r.at("metric"),r.at("axis")};
 if(!points.empty()&&k!=previous){drawTrend(metadata,points);points.clear();}previous=k;metadata=r;
 require(points.size()<budget/(4*sizeof(PlotPoint)),"One figure exceeds point memory budget; increase --memory-mb (all numerical bins remain saved)");
 points.push_back({num(r,"run"),num(r,"file"),num(r,"segment"),num(r,"stage"),std::stoll(r.at("n")),val(r,"x"),val(r,"value"),val(r,"error")});}
 if(!points.empty())drawTrend(metadata,points);
 points.clear();points.shrink_to_fit();
 auto cr=table(out/"correlations.tsv");using CK=std::tuple<int,std::string,int,std::string>;std::map<CK,std::vector<Row>>corr;for(auto&r:cr)corr[{num(r,"sample"),r.at("particle"),num(r,"stage"),r.at("group")}].push_back(r);
 for(auto&[k,v]:corr){auto[s,p,stage,group]=k;int n=group=="layer"?30:180;canvas.Clear();canvas.SetRightMargin(.16);TH2D h("correlation",";index;index",n,-.5,n-.5,n,-.5,n-.5);h.SetDirectory(nullptr);h.SetMinimum(-1);h.SetMaximum(1);for(int a=1;a<=n;++a)for(int b=1;b<=n;++b)h.SetBinContent(a,b,NaN);for(auto&r:v){double rho=val(r,"correlation");h.SetBinContent(num(r,"i")+1,num(r,"j")+1,rho);h.SetBinContent(num(r,"j")+1,num(r,"i")+1,rho);}h.SetTitle((p=="electron"?"Electron shower covariance (not electronics noise)":"Muon channel-centered response correlation"));h.Draw("COLZ");
 // ROOT paints NaN TH2 bins at the palette minimum. Explicit gray boxes are required.
 std::vector<std::unique_ptr<TBox>> missing;
 for(int iy=1;iy<=n;++iy){int ix=1;while(ix<=n){if(std::isfinite(h.GetBinContent(ix,iy))){++ix;continue;}int first=ix;while(ix<=n&&!std::isfinite(h.GetBinContent(ix,iy)))++ix;auto box=std::make_unique<TBox>(first-1.5,iy-1.5,ix-1.5,iy-.5);box->SetFillColor(kGray);box->SetLineColor(kGray);box->Draw();missing.push_back(std::move(box));}}
 TLatex missingLabel;missingLabel.SetNDC();missingLabel.SetTextSize(.023);missingLabel.DrawLatex(.15,.92,"Gray = unmeasured / insufficient pairs");
 auto path=plots/subdir(s,p=="muon"?0:1,c)/(group+"_correlation_stage"+std::to_string(stage)+".png");save(canvas,path);index<<fs::relative(path,plots).string()<<"\tPairwise Pearson correlation, common data/MC -1..1; gray is unmeasured\t"<<v.size()<<'\n';}
 auto sr=table(out/"slopes.tsv");std::map<std::tuple<int,std::string,int,std::string>,std::vector<Row>>sl;
 for(auto&r:sr)sl[{num(r,"sample"),r.at("group"),num(r,"stage"),r.at("model")}].push_back(r);
 for(auto&[key,v]:sl){auto[s,g,stage,model]=key;canvas.Clear();canvas.SetRightMargin(.05);TGraphErrors graph;for(auto&r:v)if(std::isfinite(val(r,"slope_ADC_per_C"))){int j=graph.GetN();graph.SetPoint(j,num(r,"id"),val(r,"slope_ADC_per_C"));double err=val(r,"block_cluster_error");if(!std::isfinite(err))err=val(r,"iid_error");graph.SetPointError(j,0,std::isfinite(err)?err:0);}graph.SetMarkerStyle(20);graph.SetTitle((model+";"+g+" index;ADC / degree C").c_str());if(graph.GetN())graph.Draw("AP");else note("Temperature range / statistics insufficient; see slopes.tsv");auto path=plots/subdir(s,0,c)/(g+"_slopes_stage"+std::to_string(stage)+"_"+model+".png");save(canvas,path);index<<fs::relative(path,plots).string()<<"\tChannel means removed; position/angle adjusted; fit statuses retained\t"<<graph.GetN()<<'\n';}
 // Physical strip placement, square plot area and shared limits across corresponding samples.
 for(int stage=0;stage<2;++stage){double extent=0;std::map<std::pair<int,int>,double>values;for(auto&r:sr)if(r.at("group")=="channel"&&r.at("model")=="within_channel_run_position_angle_adjusted"&&num(r,"stage")==stage&&std::isfinite(val(r,"slope_ADC_per_C"))){double v=val(r,"slope_ADC_per_C");values[{num(r,"sample"),num(r,"id")}]=v;extent=std::max(extent,std::abs(v));}extent=std::max(extent,.01);
 std::set<int>layers;for(auto&[k,v]:values)layers.insert(k.second/100000);
 for(int layer:layers)for(int sample=0;sample<2;++sample){TCanvas mapCanvas("physical_map","physical map",800,800);mapCanvas.SetLeftMargin(.13);mapCanvas.SetRightMargin(.13);mapCanvas.SetTopMargin(.13);mapCanvas.SetBottomMargin(.13);TH2D frame("frame","Within-run ADC temperature slope;X [mm];Y [mm]",1,-120,120,1,-120,120);frame.SetDirectory(nullptr);frame.SetMinimum(-extent);frame.SetMaximum(extent);frame.Draw("AXIS");std::vector<std::unique_ptr<TBox>>boxes;int measured=0;
 for(int id:channels())if(id/100000==layer){auto pos=geometry(id);double dx=layer%2?2.5:22.5,dy=layer%2?22.5:2.5;auto b=std::make_unique<TBox>(pos[0]-dx,pos[1]-dy,pos[0]+dx,pos[1]+dy);auto it=values.find({sample,id});b->SetLineColor(kGray+1);if(it==values.end())b->SetFillColor(kGray);else{++measured;int n=gStyle->GetNumberOfColors();int color=std::clamp(int((it->second+extent)/(2*extent)*(n-1)),0,n-1);b->SetFillColor(gStyle->GetColorPalette(color));}b->Draw();boxes.push_back(std::move(b));}
 TLatex label;label.SetNDC();label.SetTextSize(.026);label.DrawLatex(.13,.95,("shared color: ["+std::to_string(-extent)+", "+std::to_string(extent)+"] ADC/C; gray=unmeasured").c_str());auto path=plots/subdir(sample,0,c)/("layer"+std::to_string(layer))/("physical_slope_stage"+std::to_string(stage)+".png");save(mapCanvas,path);index<<fs::relative(path,plots).string()<<"\tActual 45x5mm strips, equal XY scale, common data/MC color limits\t"<<measured<<'\n';}}
 if(fs::exists(out/"gain.tsv")){auto gain=table(out/"gain.tsv");std::map<std::tuple<int,std::string,std::string,std::string>,std::vector<Row>>gains;for(auto&r:gain){int id=num(r,"cellid");if(!selected.empty()&&!selected.count(id))continue;if(r.at("gain").empty()||r.at("temperature_C").empty())continue;if(std::isfinite(val(r,"gain"))&&std::isfinite(val(r,"temperature_C")))gains[{id,r.at("method"),r.at("period"),r.at("dac")}].push_back(r);}
 for(auto&[k,v]:gains){auto[id,method,period,dac]=k;canvas.Clear();TGraphErrors g;for(auto&r:v){int j=g.GetN();g.SetPoint(j,val(r,"temperature_C"),val(r,"gain"));double err=r.at("stat_error").empty()?NaN:val(r,"stat_error");g.SetPointError(j,0,std::isfinite(err)?err:0);}g.SetMarkerStyle(24);g.SetTitle(("LED "+period+" DAC"+dac+" "+method+" (candidates: see gain.tsv);temperature [C];gain [ADC / p.e.]").c_str());g.Draw("AP");fs::path path=plots/"data"/"LED_gain"/period/("DAC"+dac)/("layer"+std::to_string(id/100000))/("chip"+std::to_string(id%100000/10000))/("channel"+std::to_string(id%100)+"_"+method+".png");save(canvas,path);index<<fs::relative(path,plots).string()<<"\tExisting independent LED fit; review/failed status is NOT promoted; no CT inference\t"<<v.size()<<'\n';}}
 // Impact figures are created even when physical time or held-out statistics are unavailable.
 auto im=table(out/"impact.tsv");for(int s=0;s<2;++s){canvas.Clear();canvas.Divide(3,1);int pad=0;for(auto&r:im)if(num(r,"sample")==s){canvas.cd(++pad);double before=val(r,"before"),after=val(r,"after");if(std::isfinite(before)&&std::isfinite(after)){TGraphErrors g;g.SetPoint(0,0,before);g.SetPoint(1,1,after);g.SetMarkerStyle(20);g.SetTitle((r.at("metric")+";0=before 1=held-out correction;value").c_str());g.DrawClone("AP");}else note(r.at("status"));}auto path=plots/subdir(s,1,c)/"held_out_impact.png";save(canvas,path);index<<fs::relative(path,plots).string()<<"\tUnused event blocks/runs, paired correction diagnostic; status shown if unavailable\t3\n";}
 std::cout<<"PNG render complete: "<<plots<<"; no input data or numerical ROOT modified\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
