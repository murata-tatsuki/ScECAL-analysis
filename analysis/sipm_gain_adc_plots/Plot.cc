#include "FftReconstruction.hh"
#include <TCanvas.h>
#include <TH1.h>
#include <TLatex.h>
#include <TStyle.h>
#include <array>
#include <cmath>
#include <tuple>

namespace adc_plots {
using namespace gain;
const std::array<std::string,4> methods={"multi_gaussian","envelope_comb","fft","autocorrelation"};
std::string padded(const std::string& prefix,int n){std::ostringstream s;s<<prefix<<std::setfill('0')<<std::setw(2)<<n;return s.str();}
std::string fixed(double x,int n=2){if(!std::isfinite(x))return "--";std::ostringstream s;s<<std::fixed<<std::setprecision(n)<<x;return s.str();}
std::string temperatureName(int bin,double width){
 if(bin<0)return "temperature_unavailable";
 int n=1;while(n<8&&std::abs(width*std::pow(10.,n)-std::round(width*std::pow(10.,n)))>1e-8)++n;
 return "temperature_"+fixed(bin*width,n)+"-"+fixed((bin+1)*width,n)+"C";
}
std::string statusCode(const Json& row){auto s=str(row,"status");return s=="candidate"?"C":s=="review"?"R":s=="insufficient_data"?"I":"missing";}
void label(double x,double y,const std::string& s,double size=.048,int color=kBlack){TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(size);t.SetTextColor(color);t.DrawLatex(x,y,s.c_str());}
void jsonObject(TDirectory& dir,const char* name,const Json& value){dir.cd();TObjString s(value.dump().c_str());require(s.Write(name)>0,"Cannot write ROOT metadata");}
struct Record {int sid,cell,temp,time;double entries,temperature;};
using Key=std::tuple<int,int,int>; // temperature, time bin, chip; layer and tags fixed per invocation
struct Input {
 fs::path sample;
 std::string sampleName,kind;
 Json collection,provenance=Json::object();
 double width;
 bool fullHash;
 std::map<std::string,Json> stamps;
 std::vector<std::unique_ptr<Lock>> locks;
 std::unique_ptr<TFile> raw;
 std::map<int,Record> records;
 std::map<Key,std::map<int,int>> pages;
 std::array<std::map<int,Json>,4> estimates;
 std::array<std::map<int,fs::path>,4> diagnosticPaths;
 std::map<std::string,std::unique_ptr<TFile>> diagnostics;
 void remember(const fs::path& path,const std::string& expected=""){
  auto key=path.string();if(stamps.count(key))return;
  auto fp=fingerprint(path);stamps[key]=fp;
  Json entry={{"fingerprint",fp}};
  if(!expected.empty()){
   entry["recorded_sha256"]=expected;entry["sha256_verified"]=fullHash;
   if(fullHash)require(sha256(path)==expected,"Input checksum differs: "+key);
  }
  provenance["files"][key]=entry;
 }
 void unchanged()const{for(auto& [path,fp]:stamps)require(fingerprint(path)==fp,"Input changed during drawing: "+path);}
 Input(const fs::path& root,const std::string& name,int layer,int chip,int hit,int gainTag,
       double requestedTemperature,bool verify):sampleName(name),fullHash(verify){
  int matches=0;for(auto k:{"led","electron","muon"})if(fs::is_directory(root/"data"/k/name)){sample=root/"data"/k/name;kind=k;++matches;}
  require(matches==1,"Sample not found or ambiguous: "+name);
  auto coll=sample/"collection";locks.push_back(std::make_unique<Lock>(coll/".lock",false,true));
  remember(coll/"COLLECT_COMPLETE");remember(coll/"inputs/collection.json");
  collection=readJson(coll/"inputs/collection.json");auto hash=trim(readText(coll/"COLLECT_COMPLETE"));
  require(str(collection.at("manifest_row"),"sample")==name,"Collection sample mismatch");
  width=num(collection.at("collection"),"temperature_bin");require(std::isfinite(width)&&width>0,"Invalid temperature width");
  require(num(collection.at("collection"),"memory_cell_selection")==0,"Expected memory cell 0");
  remember(coll/"spectra.root",hash);raw=openRoot(coll/"spectra.root");auto* t=tree(*raw,"Spectra");
  for(Long64_t i=0;i<t->GetEntries();++i){require(t->GetEntry(i)>0,"Cannot read Spectra");int cell=leaf(*t,"cellid"),c=cell%100000/10000,ch=cell%100;
   if(cell/100000!=layer||(chip>=0&&c!=chip)||leaf(*t,"hit_tag")!=hit||leaf(*t,"gain_tag")!=gainTag||leaf(*t,"memory_cell")!=0)continue;
   require(c>=0&&c<6&&ch>=0&&ch<(c==5?30:36)&&cell%10000/100==0,"Invalid physical CellID");
   int temp=leaf(*t,"temperature_bin"),time=leaf(*t,"time_bin");
   if(std::isfinite(requestedTemperature)&&(temp<0||temp!=int(std::floor(requestedTemperature/width))))continue;
   Record r{int(leaf(*t,"sid")),cell,temp,time,leaf(*t,"entries"),leaf(*t,"temperature_C")};
   require(records.emplace(r.sid,r).second,"Duplicate spectrum ID");
   require(pages[{temp,time,c}].emplace(ch,r.sid).second,"Duplicate channel within temperature/time condition");
  }
  require(!records.empty(),"No ADC spectra match the requested sample/layer/chip/temperature/tags");
  for(int m=0;m<4;++m){
   std::map<int,fs::path> parts;int groups=-1;std::regex pattern("part_([0-9]+)_of_([0-9]+)");
   require(fs::is_directory(sample/methods[m]),"Missing method: "+methods[m]);
   for(auto& d:fs::directory_iterator(sample/methods[m])){std::smatch match;auto s=d.path().filename().string();if(!d.is_directory()||!std::regex_match(s,match,pattern))continue;
    int group=std::stoi(match[1]),count=std::stoi(match[2]);require(count>0&&group>=0&&group<count&&(groups<0||groups==count),"Mixed fit partition counts");groups=count;
    require(parts.emplace(group,d.path()).second,"Duplicate fit partition");
   }
   require(groups>0,"No fit partitions: "+methods[m]);std::set<int> needed;
   for(auto& [sid,r]:records){(void)sid;needed.insert((layer*210+(r.cell%100000/10000)*36+r.cell%100)%groups);}
   for(int group:needed){require(parts.count(group),"Missing fit partition");auto part=parts.at(group);locks.push_back(std::make_unique<Lock>(part/".lock",false,true));
    remember(part/"FIT_COMPLETE");auto done=readJson(part/"FIT_COMPLETE");auto config=done.at("config");
    require(str(config,"method")==methods[m]&&num(config,"group")==group&&num(config,"groups")==groups&&str(config,"spectra_sha256")==hash,"Fit/collection identity mismatch");
    remember(part/"gain_estimates.tsv",done.at("files").at("gain_estimates.tsv").get<std::string>());
    remember(part/"diagnostics.root",done.at("files").at("diagnostics.root").get<std::string>());
    for(auto& row:readTable(part/"gain_estimates.tsv")){int sid=num(row,"sid");if(!records.count(sid))continue;auto& r=records.at(sid);
     require(num(row,"cellid")==r.cell&&num(row,"temperature_bin")==r.temp&&num(row,"time_bin")==r.time&&num(row,"hit_tag")==hit&&num(row,"gain_tag")==gainTag&&num(row,"memory_cell")==0&&str(row,"sample")==name&&str(row,"method")==methods[m],"Method row does not match collection spectrum");
     require(estimates[m].emplace(sid,row).second,"Duplicate method estimate");diagnosticPaths[m][sid]=part/"diagnostics.root";
    }
   }
  }
  provenance["collection"]=collection;
  provenance["verification"]=fullHash?"sha256 plus histogram-bin comparisons":"completion markers, file fingerprints, and histogram-bin comparisons; full file hashes not recomputed";
 }
 TFile& diagnostic(int method,int sid){auto p=diagnosticPaths[method].at(sid).string();auto it=diagnostics.find(p);if(it==diagnostics.end())it=diagnostics.emplace(p,openRoot(p)).first;return *it->second;}
 std::unique_ptr<TGraph> curve(int m,int sid,const TH1& rawHistogram){
  if(!estimates[m].count(sid))return {};
  auto& f=diagnostic(m,sid);auto prefix="s"+std::to_string(sid)+"/";
  auto* roi=dynamic_cast<TH1*>(f.Get((prefix+"spectrum").c_str()));
  if(roi)for(int i=1;i<=roi->GetNbinsX();++i){int bin=rawHistogram.FindFixBin(roi->GetBinCenter(i));require(bin>=1&&bin<=rawHistogram.GetNbinsX()&&std::abs(rawHistogram.GetBinContent(bin)-roi->GetBinContent(i))<1e-9,"Diagnostic ADC differs from raw histogram");}
  auto* g=dynamic_cast<TGraph*>(f.Get((prefix+methods[m]).c_str()));
  if(!g||g->GetN()==0){require(str(estimates[m].at(sid),"status")!="candidate","Candidate method has no saved curve");return {};}
  require(roi,"Saved curve has no ADC spectrum");
  for(int i=0;i<g->GetN();++i)require(std::isfinite(g->GetPointX(i))&&std::isfinite(g->GetPointY(i)),"Nonfinite saved curve");
  if(m<2)require(g->GetN()==roi->GetNbinsX(),"Fit curve/ADC dimension mismatch");
  return std::unique_ptr<TGraph>(static_cast<TGraph*>(g->Clone(methods[m].c_str())));
 }
};
std::pair<double,double> bounds(const TH1& h,const std::array<std::unique_ptr<TGraph>,4>& curves,const TGraph* fft,const std::string& mode){
 if(mode=="full")return {h.GetXaxis()->GetXmin(),h.GetXaxis()->GetXmax()};
 double lo=INFINITY,hi=-INFINITY;
 if(mode=="data")for(int i=1;i<=h.GetNbinsX();++i)if(h.GetBinContent(i)>0){lo=std::min(lo,h.GetBinCenter(i));hi=std::max(hi,h.GetBinCenter(i));}
 for(int m=0;m<2;++m)if(curves[m])for(int i=0;i<curves[m]->GetN();++i){lo=std::min(lo,curves[m]->GetPointX(i));hi=std::max(hi,curves[m]->GetPointX(i));}
 if(fft)for(int i=0;i<fft->GetN();++i){lo=std::min(lo,fft->GetPointX(i));hi=std::max(hi,fft->GetPointX(i));}
 if(!std::isfinite(lo)){if(mode=="fit")return bounds(h,curves,fft,"data");return {0,4096};}
 double pad=std::max(5.,.06*(hi-lo));return {std::max(h.GetXaxis()->GetXmin(),lo-pad),std::min(h.GetXaxis()->GetXmax(),hi+pad)};
}
void draw(Input& input,const Key& key,const fs::path& output,const Json& settings,const Json& sources,const std::vector<std::string>& command){
 auto [temp,time,chip]=key;int layer=num(settings,"layer"),hit=num(settings,"hit_tag"),gainTag=num(settings,"gain_tag");
 auto dir=output/"data"/input.kind/input.sampleName/temperatureName(temp,input.width)/("hit_"+std::to_string(hit)+"_gain_"+std::to_string(gainTag)+"_memory_0")/padded("layer_",layer);
 if(time!=0)dir/="time_bin_"+std::to_string(time);
 fs::create_directories(dir);auto stem=padded("chip_",chip);Lock lock(dir/(stem+".lock"));
 Json info={{"format_version",2},{"settings",settings},{"sources",sources},{"input",input.provenance},{"temperature_bin",temp},{"temperature_lower_C",temp<0?Json():Json(temp*input.width)},{"temperature_upper_C",temp<0?Json():Json((temp+1)*input.width)},{"time_bin",time},{"chip",chip},{"command",commandLine(command)},
  {"note","Raw ADC bin contents/errors and saved fit curves preserved; no refit, rebinning, pedestal subtraction, or calibration change. FFT band-pass reconstruction uses the saved FFT gain and original ADC phase; it is a display, not a fit. ACF gain is a label."}};
 auto marker=dir/(stem+".json");
 if(fs::is_regular_file(marker)){auto old=readJson(marker);bool same=true;for(auto k:{"format_version","settings","sources","input","temperature_bin","time_bin","chip"})if(old.at(k)!=info.at(k))same=false;
  if(same){for(auto ext:{"root","png"})require(sha256(dir/(stem+"."+ext))==old.at("outputs").at(ext),"Output checksum mismatch; remove only the chip JSON marker to regenerate");std::cout<<"Already complete: "<<dir/stem<<std::endl;return;}
 }
 auto suffix=".pending_"+std::to_string(getpid());auto rootPath=dir/(stem+suffix+".root"),pngPath=dir/(stem+suffix+".png");
 try{
  auto file=openRoot(rootPath,"RECREATE");
  std::vector<std::unique_ptr<TH1>> histograms;std::vector<std::unique_ptr<TGraph>> graphs;
  TCanvas canvas("adc_fits","Raw ADC and saved gain fits",3000,2400);canvas.Divide(6,6,.001,.001);
  for(int ch=0;ch<36;++ch){auto* pad=static_cast<TPad*>(canvas.cd(ch+1));double x=pad->GetXlowNDC(),y=pad->GetYlowNDC();pad->SetPad(x,y*.95,x+pad->GetWNDC(),(y+pad->GetHNDC())*.95);
   pad->SetLeftMargin(.13);pad->SetRightMargin(.035);pad->SetBottomMargin(.14);pad->SetTopMargin(.26);pad->SetTicks(1,1);
   int cell=layer*100000+chip*10000+ch;Json channel={{"channel",ch},{"cellid",cell}};
   const auto& page=input.pages.at(key);auto match=page.find(ch);
   if(chip==5&&ch>=30){label(.12,.93,"ch "+std::to_string(ch));label(.18,.48,"Not connected",.065,kGray+2);channel["status"]="not_connected";info["channels"].push_back(channel);continue;}
   if(match==page.end()){label(.12,.93,"ch "+std::to_string(ch)+" | "+std::to_string(cell));label(.12,.48,"No ADC in this condition",.055,kGray+2);channel["status"]="no_adc";info["channels"].push_back(channel);continue;}
   int sid=match->second;auto& record=input.records.at(sid);auto* original=dynamic_cast<TH1*>(input.raw->Get(("s"+std::to_string(sid)).c_str()));require(original&&original->GetNbinsX()==4096,"Missing raw ADC histogram");
   auto h=std::unique_ptr<TH1>(static_cast<TH1*>(original->Clone(("raw_adc_"+std::to_string(cell)).c_str())));h->SetDirectory(nullptr);
   std::array<std::unique_ptr<TGraph>,4> curves;for(int m=0;m<4;++m){curves[m]=input.curve(m,sid,*h);if(input.estimates[m].count(sid))channel["methods"][methods[m]]=input.estimates[m].at(sid);}
   FftReconstruction fft;
   if(!settings.at("fft_reconstruction").get<bool>())fft.info={{"status","disabled"}};
   else if(!curves[2])fft.info["reason"]="no_saved_fft_spectrum";
   else{
    auto& f=input.diagnostic(2,sid);auto prefix="s"+std::to_string(sid)+"/";
    auto* roi=dynamic_cast<TH1*>(f.Get((prefix+"spectrum").c_str()));auto* env=dynamic_cast<TGraph*>(f.Get((prefix+"env").c_str()));
    require(roi&&env,"Missing saved FFT ADC/envelope");fft=reconstructFft(*roi,*env,*curves[2],input.estimates[2].at(sid));
   }
   channel["fft_reconstruction"]=fft.info;
   auto range=bounds(*h,curves,fft.curve.get(),str(settings,"adc_range"));double ymax=1,ymin=0;
   for(int i=1;i<=h->GetNbinsX();++i)if(h->GetBinCenter(i)>=range.first&&h->GetBinCenter(i)<=range.second)ymax=std::max(ymax,h->GetBinContent(i)+h->GetBinError(i));
   for(int m=0;m<2;++m)if(curves[m])for(int i=0;i<curves[m]->GetN();++i)ymax=std::max(ymax,curves[m]->GetPointY(i));
   if(fft.curve)for(int i=0;i<fft.curve->GetN();++i){ymax=std::max(ymax,fft.curve->GetPointY(i));ymin=std::min(ymin,fft.curve->GetPointY(i));}
   h->SetTitle(";HG ADC;Entries / ADC");h->SetStats(false);h->SetLineColor(kBlack);h->SetMarkerColor(kBlack);h->SetMarkerStyle(20);h->SetMarkerSize(.38);h->SetMinimum(ymin<0?ymin-.05*(ymax-ymin):0);h->SetMaximum(ymax*1.12);h->GetXaxis()->SetRangeUser(range.first,range.second);
   for(auto* axis:{h->GetXaxis(),h->GetYaxis()}){axis->SetLabelSize(.048);axis->SetTitleSize(.052);axis->SetNdivisions(505);}h->GetYaxis()->SetTitleOffset(1.05);h->Draw("E");
   for(int m=0;m<2;++m)if(curves[m]){curves[m]->SetLineColor(m==0?kBlue+1:kGreen+2);curves[m]->SetLineWidth(2);curves[m]->SetLineStyle(str(input.estimates[m].at(sid),"status")=="candidate"?1:2);curves[m]->Draw("L SAME");}
   if(fft.curve){fft.curve->SetLineColor(kRed+1);fft.curve->SetLineWidth(2);fft.curve->SetLineStyle(str(input.estimates[2].at(sid),"status")=="candidate"?1:2);fft.curve->Draw("L SAME");}
   label(.13,.94,"ch "+std::to_string(ch)+" | cell "+std::to_string(cell)+" | N="+fixed(record.entries,0),.046);
   for(int m=0;m<4;++m){if(m>=2&&!settings.at("spectral_labels").get<bool>())continue;Json r=input.estimates[m].count(sid)?input.estimates[m].at(sid):Json::object();const char* shortNames[]={"MG","EC","FFT","ACF"};int colors[]={kBlue+1,kGreen+2,kRed+1,kMagenta+1};
    label(m%2==0?.13:.56,m<2?.865:.795,std::string(shortNames[m])+": "+fixed(num(r,"gain"))+" ["+statusCode(r)+"]",.043,colors[m]);}
   auto* sub=file->mkdir(("channel_"+std::to_string(ch)).c_str());require(sub,"Cannot create channel ROOT directory");sub->cd();require(h->Write("raw_adc")>0,"Cannot write raw ADC");
   for(int m=0;m<4;++m)if(curves[m]){require(curves[m]->Write(methods[m].c_str())>0,"Cannot write saved curve");graphs.push_back(std::move(curves[m]));}
   for(auto* g:{&fft.curve,&fft.envelope,&fft.filter,&fft.windowedResidual})if(*g){require((*g)->Write()>0,"Cannot write FFT reconstruction");graphs.push_back(std::move(*g));}
   channel["sid"]=sid;channel["status"]="drawn";channel["raw_entries"]=record.entries;channel["temperature_C"]=record.temperature;channel["display_adc_min"]=range.first;channel["display_adc_max"]=range.second;
   jsonObject(*sub,"metadata_json",channel);info["channels"].push_back(channel);histograms.push_back(std::move(h));
  }
  canvas.cd();label(.012,.981,input.sampleName+" | "+padded("layer ",layer)+" chip "+std::to_string(chip)+" | "+(temp<0?"Temperature unavailable":"T: "+fixed(temp*input.width)+"-"+fixed((temp+1)*input.width)+" ^{#circ}C")+" | HitTag="+std::to_string(hit)+" GainTag="+std::to_string(gainTag)+" Memory=0",.014);
  label(.012,.957,"Black: raw ADC; blue: multi-Gaussian; green: envelope-comb; red: FFT reconstruction (not fit) | C: candidate (solid), R: review (dashed), I: insufficient | ADC/p.e.",.012);
  file->cd();jsonObject(*file,"metadata_json",info);require(canvas.Write("adc_fits")>0,"Cannot write canvas");canvas.SaveAs(pngPath.c_str());require(fs::is_regular_file(pngPath)&&fs::file_size(pngPath)>0,"PNG was not written");
  require(!file->TestBit(TFile::kWriteError),"ROOT write error");file->Close();canvas.Close();input.unchanged();
  info["outputs"]={{"root",sha256(rootPath)},{"png",sha256(pngPath)}};
  fs::rename(rootPath,dir/(stem+".root"));fs::rename(pngPath,dir/(stem+".png"));auto pending=dir/(stem+suffix+".json");writeJson(pending,info);fs::rename(pending,marker);
  std::cout<<"Completed: "<<dir/(stem+".png")<<std::endl;
 }catch(...){fs::remove(rootPath);fs::remove(pngPath);throw;}
}
}
int main(int argc,char** argv)try{
 using namespace adc_plots;
 Options o;std::vector<std::string> command;for(int i=0;i<argc;++i)command.emplace_back(argv[i]);
 std::set<std::string> allowed={"label","output-label","sample","layer","chip","temperature","hit-tag","gain-tag","adc-range","spectral-labels","fft-reconstruction","verify-inputs","project"};
 for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--help"){
  std::cout<<"Plot --label LABEL --sample SAMPLE --layer 0..29 [--chip 0..5]\n"
   "  --temperature C          select the saved interval containing C; default all\n"
   "  --adc-range data|fit|full default data: all populated ADC bins + saved fit range\n"
   "  --spectral-labels 0|1     default 1: annotate FFT/ACF gains\n"
   "  --fft-reconstruction 0|1 default 1: overlay gain-selected inverse FFT (red)\n"
   "  --hit-tag 0|1 --gain-tag 0|1   defaults 1; memory cell always 0; raw branch HG_Charge\n"
   "  --verify-inputs metadata|sha256  default metadata; sha256 checks all consumed files\n"
   "  --output-label NAME --project DIRECTORY\n"
   "No refit. One 6x6 PNG and ROOT canvas per chip and saved temperature/time condition.\n";return 0;}
  require(a.rfind("--",0)==0&&allowed.count(a.substr(2))&&i+1<argc,"Unknown/missing option: "+a);a=a.substr(2);require(!o.has(a),"Duplicate option");o.values[a]=argv[++i];
 }
 auto here=fs::canonical("/proc/self/exe").parent_path();auto project=o.has("project")?fs::canonical(o.required("project")):here.parent_path().parent_path();
 auto labelName=token(o.required("label")),outputLabel=token(o.get("output-label",labelName)),sample=token(o.required("sample"));
 require(o.has("layer"),"Specify --layer 0..29");int layer=o.integer("layer",-1),chip=o.integer("chip",-1),hit=o.integer("hit-tag",1),gainTag=o.integer("gain-tag",1),spectral=o.integer("spectral-labels",1),fftOverlay=o.integer("fft-reconstruction",1);
 require(layer>=0&&layer<30&&chip>=-1&&chip<6&&(hit==0||hit==1)&&(gainTag==0||gainTag==1)&&(spectral==0||spectral==1)&&(fftOverlay==0||fftOverlay==1),"Invalid layer/chip/tag option");
 auto range=o.get("adc-range","data"),verification=o.get("verify-inputs","metadata");require(range=="data"||range=="fit"||range=="full","Invalid --adc-range");require(verification=="metadata"||verification=="sha256","Invalid --verify-inputs");
 auto output=project/"analysis/result/sipm_gain_adc_plots"/outputLabel;fs::create_directories(output/"inputs");
 Json sources=Json::object();for(auto path:{here/"Plot.cc",here/"FftReconstruction.hh",here/"Makefile",here.parent_path()/"sipm_gain/common/GainIO.hh",here.parent_path()/"sipm_gain/common/GainMethods.hh"})sources[fs::relative(path,project).string()]=sha256(path);
 auto executableHash=sha256(here/"Plot");sources["executable_sha256"]=executableHash;
 {Lock lock(output/"inputs/.lock",true);for(auto it=sources.begin();it!=sources.end();++it){if(it.key()=="executable_sha256")continue;auto dest=output/"inputs/source"/it.key();fs::create_directories(dest.parent_path());fs::copy_file(project/it.key(),dest,fs::copy_options::overwrite_existing);}writeJson(output/"inputs/rendering_sources.json",sources);}
 Json settings={{"label",labelName},{"output_label",outputLabel},{"sample",sample},{"layer",layer},{"hit_tag",hit},{"gain_tag",gainTag},{"adc_range",range},{"spectral_labels",spectral==1},{"fft_reconstruction",fftOverlay==1},{"verify_inputs",verification}};
 gROOT->SetBatch(kTRUE);gROOT->SetStyle("Plain");gStyle->SetOptStat(0);gStyle->SetOptTitle(0);TH1::AddDirectory(kFALSE);
 Input input(project/"analysis/result/sipm_gain"/labelName,sample,layer,chip,hit,gainTag,o.real("temperature",NaN),verification=="sha256");
 for(auto& [key,page]:input.pages){(void)page;draw(input,key,output,settings,sources,command);}return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
