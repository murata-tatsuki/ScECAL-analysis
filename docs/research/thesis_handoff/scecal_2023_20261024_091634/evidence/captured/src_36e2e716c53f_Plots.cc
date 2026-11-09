#include "Common.hh"
#include <TBox.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>

namespace selected_gain {
namespace {
void label(double x,double y,const std::string& text,double size=.026){TLatex l;l.SetNDC();l.SetTextFont(42);l.SetTextSize(size);l.DrawLatex(x,y,text.c_str());}
void overviewPad(){gPad->SetTopMargin(.23);gPad->SetBottomMargin(.18);gPad->SetRightMargin(.03);}
void overviewAxes(TH1* frame){for(auto* axis:{frame->GetXaxis(),frame->GetYaxis()}){axis->SetLabelSize(.055);axis->SetTitleSize(.06);}frame->GetYaxis()->SetTitleOffset(.55);frame->GetXaxis()->SetTitleOffset(1.0);}
void save(TCanvas& canvas,TFile& file,const std::string& directory,const fs::path& png) {
 fs::create_directories(png.parent_path());canvas.SaveAs(png.c_str());require(fs::is_regular_file(png)&&fs::file_size(png)>0,"Plot PNG not written");
 auto* dir=file.GetDirectory(directory.c_str());if(!dir)dir=file.mkdir(directory.c_str());require(dir,"Cannot create plot ROOT directory");dir->cd();require(canvas.Write()>0,"Cannot write plot canvas");file.cd();
}
std::string temperatureText(const Json& interval) {
 if(!interval.value("valid",false))return "Temperature unavailable";
 const double width=num(interval,"width_C");int precision=1;
 while(precision<8&&std::abs(width*std::pow(10.,precision)-std::round(width*std::pow(10.,precision)))>1e-8)++precision;
 std::ostringstream text;text<<std::fixed<<std::setprecision(precision)
  <<"T: "<<num(interval,"lower_C")<<"-"<<num(interval,"upper_C")<<" ^{#circ}C";
 return text.str();
}
std::string shortCondition(const Json& c,const Json& temperatures){return str(c,"beam")+" "+str(c,"period")+"; "+temperatureText(temperatures.at(str(c,"condition")))+"; settings="+str(c,"settings");}
// Read the actual collection width rather than assuming the default 0.5 C.
Json temperatureMetadata(const Paths& p,const std::map<std::string,std::map<int,Json>>& decisions,
 const std::map<std::pair<std::string,int>,std::vector<Json>>& evidence) {
 const auto comparison=readJson(p.input/"inputs/configuration.json");
 const auto& frozen=comparison.at("input_snapshot").at("collections");
 Json samples=Json::object(),intervals=Json::object();
 for(const auto& [key,rows]:evidence){const auto& d=decisions.at(key.first).at(key.second);
  for(const auto& e:rows){const auto sample=token(str(e,"sample"));
   if(!samples.contains(sample)){
    const auto directory=p.project/"analysis/result/sipm_gain"/p.label/"data/led"/sample/"collection";
    Lock lock(directory/".lock",false,true);
    require(frozen.at(sample)==sha256(directory/"COLLECT_COMPLETE"),"Collection differs from the compared input: "+sample);
    const auto path=directory/"inputs/collection.json";const auto text=readText(path);const auto record=Json::parse(text);
    const auto& identity=record.at("manifest_row");double width=num(record.at("collection"),"temperature_bin");
    require(std::isfinite(width)&&width>0,"Invalid temperature-bin width: "+sample);
    require(str(identity,"sample")==sample&&str(identity,"kind")=="led","Unexpected LED collection identity");
    samples[sample]={{"path",path.string()},{"sha256",digest(text)},{"width_C",width},{"beam",str(identity,"beam")},{"period",str(identity,"period")}};
   }
   const auto& input=samples.at(sample);require(str(input,"beam")==str(d,"beam")&&str(input,"period")==str(d,"period"),"Temperature source condition mismatch");
   const double bin=num(d,"temperature_bin"),width=num(input,"width_C");
   require(std::isfinite(bin)&&bin==std::floor(bin)&&bin>=-1,"Invalid temperature-bin index");
   Json interval={{"valid",bin>=0},{"width_C",width},{"lower_C",bin>=0?Json(bin*width):Json()},{"upper_C",bin>=0?Json((bin+1)*width):Json()}};
   if(intervals.contains(key.first))require(intervals.at(key.first)==interval,"Mixed temperature-bin widths within one condition");
   else intervals[key.first]=interval;
   if(bin>=0){double t=num(e,"temperature_C");require(std::isfinite(t)&&t>=bin*width-1e-7&&t<=(bin+1)*width+1e-7,"Recorded temperature is outside the labelled interval");}
  }
 }
 for(const auto& [key,rows]:decisions){(void)rows;require(intervals.contains(key),"Missing temperature provenance for condition "+key);}
 return {{"conditions",intervals},{"collection_inputs",samples},{"interval_convention","lower inclusive, upper exclusive; floor(temperature_C / width_C)"}};
}
}
void plots(const Paths& p,const Json& config,int layer,int limit,const std::string& scope,const std::vector<std::string>& command) {
 auto name=layerName(layer);auto input=p.out/"layers"/name;Lock inputLock(input/".lock",false,true);auto signature=configHash(p);
 verifyComplete(input,"SELECTION_COMPLETE",signature);auto inputSHA=sha256(input/"SELECTION_COMPLETE");
 auto out=p.out/"plots/layers"/name;fs::create_directories(out);Lock lock(out/".lock");
 Json previous;
 if(fs::exists(out/"PLOTS_COMPLETE")){
  previous=verifyComplete(out,"PLOTS_COMPLETE",signature);require(previous.at("input_marker_sha256")==inputSHA&&previous.at("plot_limit")==limit&&previous.at("plot_scope")==scope,"Plot input/limit/scope changed; use the original plot options for this output label");
 }
 auto source=openRoot(input/"selection.root");std::map<std::string,std::map<int,Json>> decisions;std::map<std::pair<std::string,int>,std::vector<Json>> evidence;
 {auto* t=tree(*source,"Decisions");Reader reader(*t,DecisionSchema);for(Long64_t i=0;i<t->GetEntries();++i){auto r=reader.at(i);require(physicalCell(r)/100000==layer,"Wrong plot layer");require(decisions[str(r,"condition")].emplace(physicalCell(r),r).second,"Duplicate plot decision");}}
 {auto* t=tree(*source,"Evidence");Reader reader(*t,EvidenceSchema);for(Long64_t i=0;i<t->GetEntries();++i){auto r=reader.at(i);evidence[{str(r,"condition"),physicalCell(r)}].push_back(std::move(r));}}
 const auto temperature=temperatureMetadata(p,decisions,evidence);const auto& intervals=temperature.at("conditions");
 Json rendering={{"format_version",1},{"temperature_display","celsius_interval"},{"temperature",temperature},{"sources",sourceChecks(p)}};
 const auto renderingSHA=digest(rendering.dump());
 if(!previous.is_null()&&str(previous,"rendering_sha256")==renderingSHA){std::cout<<"Already complete: plots "<<name<<std::endl;return;}
 fs::remove(out/"PLOTS_COMPLETE");std::vector<std::string> files={"plots.root","plot_index.tsv","plot_configuration.json"};
 for(auto it=rendering.at("sources").begin();it!=rendering.at("sources").end();++it){const auto relative=fs::path("inputs/source")/it.key(),dest=out/relative;fs::create_directories(dest.parent_path());fs::copy_file(p.project/it.key(),dest,fs::copy_options::overwrite_existing);require(sha256(dest)==it.value(),"Rendering source changed while snapshotting");files.push_back(relative.string());}
 atomicJson(out/"plot_configuration.json",rendering);
 auto target=openRoot(out/"plots.root","RECREATE");Json index=Json::array();Settings settings(config);gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
 const std::array<int,4> colors={kBlue+1,kGreen+2,kRed+1,kMagenta+1};int channelPlots=0;
 auto record=[&](const fs::path& path,const std::string& condition,int cell,const std::string& kind){files.push_back(fs::relative(path,out).string());index.push_back({{"figure",fs::relative(path,p.out).string()},{"condition",condition},{"cellid",cell<0?Json():Json(cell)},{"kind",kind}});};
 for(auto& [key,rows]:decisions){auto c=conditionFields(rows.begin()->second);auto pngDir=p.out/"plots"/conditionPath(c,intervals.at(key))/name;
  TCanvas canvas("overview","Selected gain and coverage",1400,900);canvas.Divide(1,3);
  TGraphErrors gains;TGraphErrors uncertainty;std::array<TGraph,4> coverage;int n=0;
  for(int chip=0;chip<6;++chip)for(int ch=0;ch<(chip==5?30:36);++ch){int cell=layer*100000+chip*10000+ch,x=chip*36+ch,status=0;auto it=rows.find(cell);
   if(settings.dead.count(cell))status=3;
   else if(it!=rows.end())status=num(it->second,"accepted")==1?2:1;
   auto& g=coverage[status];g.SetPoint(g.GetN(),x,status);
   if(status==2){auto& r=it->second;gains.SetPoint(n,x,num(r,"gain"));gains.SetPointError(n,0,num(r,"gain_error"));uncertainty.SetPoint(n,x,100*num(r,"gain_error")/num(r,"gain"));++n;}
  }
  canvas.cd(1);overviewPad();double ymax=1;for(int i=0;i<gains.GetN();++i)ymax=std::max(ymax,gains.GetPointY(i)+gains.GetErrorY(i));
  overviewAxes(gPad->DrawFrame(-1,0,210,ymax*1.2,";chip*36+channel;Gain [ADC / p.e.]"));gains.SetName("selected_gain");gains.SetMarkerStyle(20);gains.SetMarkerSize(.65);gains.SetMarkerColor(kBlue+1);if(gains.GetN())gains.Draw("P SAME");
  label(.10,.92,name+" | "+shortCondition(c,intervals),.055);label(.10,.82,"Bars: conservative quality scale (not statistical #sigma)",.055);
  canvas.cd(2);overviewPad();double emax=5;for(int i=0;i<uncertainty.GetN();++i)emax=std::max(emax,uncertainty.GetPointY(i));overviewAxes(gPad->DrawFrame(-1,0,210,emax*1.2,";chip*36+channel;gainError / gain [%]"));uncertainty.SetName("relative_quality_scale");uncertainty.SetMarkerStyle(20);uncertainty.SetMarkerSize(.65);if(uncertainty.GetN())uncertainty.Draw("P SAME");
  canvas.cd(3);overviewPad();auto* frame=gPad->DrawFrame(-1,-.5,210,3.5,";chip*36+channel;Decision code");overviewAxes(frame);frame->GetYaxis()->SetNdivisions(4);
  const std::array<int,4> statusColors={kGray+1,kOrange+7,kBlue+1,kBlack};for(int i=0;i<4;++i){coverage[i].SetName(("decision_code_"+std::to_string(i)).c_str());coverage[i].SetMarkerStyle(20);coverage[i].SetMarkerColor(statusColors[i]);if(coverage[i].GetN())coverage[i].Draw("P SAME");}
  label(.10,.86,"0: no LED data   1: unresolved   2: selected   3: explicit dead mask",.055);
  save(canvas,*target,key,pngDir/"overview.png");record(pngDir/"overview.png",key,-1,"overview");
 }
 // Plot adopted channels first in a limited validation run, then unresolved ones.
 for(int selectedPass=1;selectedPass>=0;--selectedPass)for(auto& [key,rows]:decisions)for(auto& [cell,d]:rows){
  if((num(d,"accepted")==1)!=bool(selectedPass)||(limit&&channelPlots>=limit))continue;
  if(!selectedPass&&scope=="selected")continue;
  auto& values=evidence[{key,cell}];std::array<TGraph,4> candidates,reviews;TGraph used;double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
  for(auto& e:values){double x=num(e,"dac");if(!std::isfinite(x))continue;
   for(int m=0;m<4;++m){double y=num(e,Methods[m]+"_gain");if(!std::isfinite(y)||y<=0)continue;auto& g=str(e,Methods[m]+"_status")=="candidate"?candidates[m]:reviews[m];g.SetPoint(g.GetN(),x,y);xmin=std::min(xmin,x);xmax=std::max(xmax,x);ymin=std::min(ymin,y);ymax=std::max(ymax,y);}
   if(num(e,"used")==1)used.SetPoint(used.GetN(),x,num(e,"primary_gain"));
  }
  if(!std::isfinite(xmin)){xmin=0;xmax=1;ymin=0;ymax=1;}
  double dx=std::max(1.,(xmax-xmin)*.07),dy=std::max(1.,(ymax-ymin)*.13);xmin-=dx;xmax+=dx;ymin=std::max(0.,ymin-dy);ymax+=dy;
  auto canvasName="channel_"+std::to_string(cell);TCanvas canvas(canvasName.c_str(),"Gain versus LED DAC",1400,900);canvas.SetTopMargin(.19);canvas.SetRightMargin(.23);
  gPad->DrawFrame(xmin,ymin,xmax,ymax,("Physical CellID "+std::to_string(cell)+";LED DAC (internal setting);Gain [ADC / p.e.]").c_str());
  TBox band;TLine line;if(selectedPass){double g=num(d,"gain"),error=num(d,"gain_error");band.SetFillColorAlpha(kAzure-9,.35);band.DrawBox(xmin,g-error,xmax,g+error);line.SetLineColor(kBlue+3);line.SetLineWidth(2);line.DrawLine(xmin,g,xmax,g);}
  TLegend legend(.78,.34,.995,.80);legend.SetTextSize(.022);legend.SetBorderSize(0);
  for(int m=0;m<4;++m){candidates[m].SetName((Methods[m]+"_candidate").c_str());reviews[m].SetName((Methods[m]+"_review").c_str());candidates[m].SetMarkerColor(colors[m]);reviews[m].SetMarkerColorAlpha(colors[m],.45);candidates[m].SetMarkerStyle(20+m);reviews[m].SetMarkerStyle(24+m);candidates[m].SetMarkerSize(1);reviews[m].SetMarkerSize(.8);if(reviews[m].GetN())reviews[m].Draw("P SAME");if(candidates[m].GetN())candidates[m].Draw("P SAME");legend.AddEntry(&candidates[m],Methods[m].c_str(),"p");}
  used.SetName("used_primary");used.SetMarkerStyle(24);used.SetMarkerSize(1.65);used.SetMarkerColor(kBlack);if(used.GetN())used.Draw("P SAME");legend.AddEntry(&used,"used primary","p");legend.Draw();
  label(.10,.97,"Physical CellID "+std::to_string(cell)+" | "+shortCondition(d,intervals),.024);label(.10,.92,"Status: "+str(d,"status")+"; qualified DACs="+integerText(num(d,"n_dacs"))+"; beam="+str(d,"beam_check_status"),.025);
  label(.10,.87,selectedPass?"Band: gainError quality scale (not statistical #sigma)":"Reason: "+str(d,"reasons"),.024);label(.78,.27,"Filled: candidate",.022);label(.78,.23,"Open/faint: review",.022);
  auto path=p.out/"plots"/conditionPath(d,intervals.at(key))/name/(canvasName+".png");save(canvas,*target,key,path);record(path,key,cell,selectedPass?"selected":"unresolved");++channelPlots;
 }
 auto meta=metadataInfo(config,command);meta["layer"]=layer;meta["plot_limit"]=limit;meta["plot_scope"]=scope;meta["channel_plots"]=channelPlots;meta["input_marker_sha256"]=inputSHA;meta["rendering"]=rendering;meta["rendering_sha256"]=renderingSHA;metadata(*target,meta);closeRoot(*target);
 writeTable(out/"plot_index.tsv",index,{"figure","condition","cellid","kind"});
 seal(out,"PLOTS_COMPLETE",signature,files,{{"layer",layer},{"input_marker_sha256",inputSHA},{"plot_limit",limit},{"plot_scope",scope},{"channel_plots",channelPlots},{"overview_plots",decisions.size()},{"rendering_sha256",renderingSHA},{"command",commandLine(command)}});
 std::cout<<"Completed plots: "<<name<<", "<<decisions.size()<<" overviews + "<<channelPlots<<" channel plots"<<std::endl;
}
}
