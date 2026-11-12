#pragma once
// Fit rows are ordered by the same ten-dimensional key as the merged histogram
// stream. Only attempted fits are reused; all cheap statistics/status checks are
// recomputed. Shell verifies the collected input, configuration and definitions.
class FitReuse {
 std::ifstream input;std::unique_ptr<char[]> inputBuffer;bool enabled=false,have=false;Key pending{};Fitted fitted;
 double nall=0,under=0,over=0,low=0,high=0;unsigned long long scanned=0;
 static int axis(const std::string&s){for(int i=0;i<6;++i)if(s==ts::axisName(i))return i;throw std::runtime_error("Invalid cached fit axis");}
 void next(){have=false;std::string line;while(std::getline(input,line)){
   if(input.eof()){std::cout<<"Ignoring incomplete final fit row"<<std::endl;return;}
   ++scanned;auto end=line.rfind('\t');ts::require(end!=std::string::npos,"Invalid fit cache row");auto cov=line.rfind('\t',end-1),code=line.rfind('\t',cov-1);ts::require(code!=std::string::npos,"Invalid fit cache row");if(line.compare(code+1,cov-code-1,"-999")==0)continue;
   auto a=ts::split(line);ts::require(a.size()==27,"Invalid cached fit columns");ts::require(a[1]=="muon","Non-muon cached fit");int group=a[2]=="channel"?0:a[2]=="chip"?1:a[2]=="layer"?2:-1;ts::require(group>=0,"Invalid cached group");Key k={std::stoi(a[0]),0,group,std::stoi(a[3]),axis(a[4]),std::stoi(a[5]),std::stoi(a[6]),std::stoi(a[7]),std::stoi(a[8]),std::stoi(a[9])};if(!includeTickTrends&&k[4]==4){++excluded;continue;}ts::require(!havePrevious||previous<k,"Cached fits are not uniquely ordered");previous=k;havePrevious=true;pending=k;
   nall=std::stod(a[11]);fitted.nfit=std::stoll(a[12]);under=std::stod(a[13]);over=std::stod(a[14]);low=std::stod(a[15]);high=std::stod(a[16]);fitted.mpv=std::stod(a[17]);fitted.me=std::stod(a[18]);fitted.w=std::stod(a[19]);fitted.we=std::stod(a[20]);fitted.sigma=std::stod(a[21]);fitted.se=std::stod(a[22]);fitted.chi=std::stod(a[23]);fitted.code=std::stoi(a[24]);fitted.cov=std::stoi(a[25]);fitted.status=a[26];have=true;return;
  }ts::require(!input.bad(),"Cannot read fit cache");
 }
 Key previous{};bool havePrevious=false;
 static bool same(double a,double b){return (std::isnan(a)&&std::isnan(b))||(a==b)||std::abs(a-b)<=1e-9+1e-10*std::max(std::abs(a),std::abs(b));}
 public:unsigned long long reused=0,excluded=0;
 explicit FitReuse(const std::string&path){if(path.empty()||path=="-")return;enabled=true;inputBuffer=std::make_unique<char[]>(1024*1024);input.rdbuf()->pubsetbuf(inputBuffer.get(),1024*1024);input.open(path);ts::require(bool(input),"Cannot open fit cache");std::string header;std::getline(input,header);ts::require(header.find("\troot_status\tcovariance\tstatus")!=std::string::npos,"Invalid fit cache header");next();}
 ~FitReuse(){if(input.is_open())input.close();}
 bool get(const Key&k,TH1D&h,const Summary&b,const ts::Config&c,Fitted&r){if(!have)return false;ts::require(!(pending<k),"Cached fit key is missing from current aggregate");if(k!=pending)return false;
  ts::require(nall==h.GetEntries()&&under==h.GetBinContent(0)&&over==h.GetBinContent(h.GetNbinsX()+1)&&same(low,b.low)&&same(high,b.high)&&b.threshold,"Cached fit histogram summary mismatch");double lo=std::max(5.,b.low),hi=std::min(2495.,b.high);long long nfit=0;for(int i=h.FindFixBin(lo);i<=h.FindFixBin(hi);++i)nfit+=std::llround(h.GetBinContent(i));ts::require(nfit==fitted.nfit&&nfit>=c.integer("min-fit-hits",100)&&hi-lo>=100,"Cached fit plateau mismatch");r=fitted;++reused;next();return true;
 }
 void finish(){ts::require(!have,"Unused cached fit key");if(enabled)std::cout<<"Reused "<<reused<<" attempted fits from "<<scanned<<" complete cached rows; excluded-axis attempted fits "<<excluded<<std::endl;}
};
class AnalysisProgress {
 ts::BufferedOutput output;std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
 public:explicit AnalysisProgress(const ts::fs::path&p):output(p){output<<"phase\tcompleted\ttotal\tnew_fit_attempts\treused_fit_attempts\telapsed_seconds\n";}
 void write(const std::string&phase,long long done,long long total,long long attempted=0,long long reused=0){output<<phase<<'\t'<<done<<'\t'<<total<<'\t'<<attempted<<'\t'<<reused<<'\t'<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';output.flush();ts::require(bool(output),"Cannot save analysis progress");std::cout<<phase<<" "<<done<<" / "<<total<<"; new fits "<<attempted<<", reused "<<reused<<std::endl;}
};
