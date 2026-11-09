#pragma once
#include "GainMethods.hh"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TROOT.h>
#include <openssl/evp.h>
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cstdlib>
#include <TObjString.h>

namespace gain {
namespace fs=std::filesystem;
inline const std::vector<std::string> ManifestFields={"sample","kind","beam","energy_GeV","period","led_group","dac","settings","raw"};
inline void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
inline std::string readText(const fs::path& path){std::ifstream f(path);require(bool(f),"Cannot read "+path.string());return std::string(std::istreambuf_iterator<char>(f),{});}
inline void writeText(const fs::path& path,const std::string& text){if(!path.parent_path().empty())fs::create_directories(path.parent_path());std::ofstream f(path);require(bool(f),"Cannot write "+path.string());f<<text;f.close();require(bool(f),"Write failed "+path.string());}
inline std::string trim(std::string v){while(!v.empty()&&std::isspace(static_cast<unsigned char>(v.back())))v.pop_back();size_t first=v.find_first_not_of(" \t\r\n");return first==std::string::npos?"":v.substr(first);}
inline std::vector<std::string> split(const std::string& s,char delimiter=','){std::vector<std::string> v;std::string item;std::istringstream f(s);while(std::getline(f,item,delimiter))v.push_back(item);return v;}
inline std::string str(const Json& j,const std::string& key,const std::string& fallback=""){if(!j.contains(key)||j[key].is_null())return fallback;return j[key].is_string()?j[key].get<std::string>():j[key].dump();}
inline double num(const Json& j,const std::string& key,double fallback=NaN){if(!j.contains(key)||j[key].is_null())return fallback;if(j[key].is_number())return j[key].get<double>();try{size_t n=0;auto s=j[key].get<std::string>();double v=std::stod(s,&n);return n==s.size()?v:fallback;}catch(...){return fallback;}}
inline std::string token(const std::string& s){require(std::regex_match(s,std::regex("[A-Za-z0-9_.-]+"))&&s!="."&&s!="..","Unsafe name: "+s);return s;}
inline void writeJson(const fs::path& p,const Json& j){writeText(p,j.dump(2)+'\n');}
inline Json readJson(const fs::path& p){return Json::parse(readText(p));}
inline std::vector<std::string> tsvLine(const std::string& line){std::vector<std::string> out;std::string value;bool quoted=false;for(size_t i=0;i<line.size();++i){char c=line[i];if(c=='"'){if(quoted&&i+1<line.size()&&line[i+1]=='"'){value+='"';++i;}else quoted=!quoted;}else if(c=='\t'&&!quoted){out.push_back(value);value.clear();}else value+=c;}require(!quoted,"Multiline/invalid TSV field");out.push_back(value);return out;}
inline Json readTable(const fs::path& p){std::ifstream f(p);require(bool(f),"Cannot read TSV "+p.string());std::string line;require(bool(std::getline(f,line)),"Empty TSV "+p.string());if(!line.empty()&&line.back()=='\r')line.pop_back();auto names=tsvLine(line);require(std::set<std::string>(names.begin(),names.end()).size()==names.size(),"Duplicate TSV columns");Json rows=Json::array();while(std::getline(f,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();if(line.empty())continue;auto v=tsvLine(line);require(v.size()==names.size(),"TSV column count differs: "+p.string());Json r;for(size_t i=0;i<v.size();++i)r[names[i]]=v[i];rows.push_back(r);}return rows;}
inline std::string tsvValue(const Json& value){std::string s=(value.is_null()||(value.is_number_float()&&!std::isfinite(value.get<double>())))?"nan":value.is_string()?value.get<std::string>():value.dump();if(s.find_first_of("\t\n\r\"")!=std::string::npos){std::string escaped="\"";for(char c:s){if(c=='"')escaped+='"';escaped+=c;}s=escaped+'"';}return s;}
inline void writeTable(const fs::path& p,const Json& rows,std::vector<std::string> fields={}){
    if(fields.empty()){std::set<std::string> seen;for(const auto& r:rows)for(auto it=r.begin();it!=r.end();++it)if(seen.insert(it.key()).second)fields.push_back(it.key());}
    std::ostringstream s;for(size_t i=0;i<fields.size();++i)s<<(i?"\t":"")<<fields[i];s<<'\n';for(const auto& r:rows){for(size_t i=0;i<fields.size();++i)s<<(i?"\t":"")<<(r.contains(fields[i])?tsvValue(r[fields[i]]):"");s<<'\n';}writeText(p,s.str());
}
inline std::string sha256(const fs::path& path){
    std::ifstream f(path,std::ios::binary);require(bool(f),"Cannot checksum "+path.string());
    std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> ctx(EVP_MD_CTX_new(),EVP_MD_CTX_free);require(bool(ctx),"SHA256 allocation failed");
    require(EVP_DigestInit_ex(ctx.get(),EVP_sha256(),nullptr)==1,"SHA256 init failed");char buffer[65536];while(f){f.read(buffer,sizeof(buffer));if(f.gcount())require(EVP_DigestUpdate(ctx.get(),buffer,f.gcount())==1,"SHA256 update failed");}require(f.eof(),"Checksum read error");
    unsigned char bytes[EVP_MAX_MD_SIZE];unsigned n=0;require(EVP_DigestFinal_ex(ctx.get(),bytes,&n)==1,"SHA256 final failed");std::ostringstream out;out<<std::hex<<std::setfill('0');for(unsigned i=0;i<n;++i)out<<std::setw(2)<<unsigned(bytes[i]);return out.str();
}
inline Json fingerprint(const fs::path& p){auto path=fs::canonical(p);struct stat s{};require(::stat(path.c_str(),&s)==0,"Cannot stat input");return {{"path",path.string()},{"bytes",s.st_size},{"mtime_ns",int64_t(s.st_mtim.tv_sec)*1000000000LL+s.st_mtim.tv_nsec}};}
inline Json manifest(const fs::path& path){auto p=fs::canonical(path);auto rows=readTable(p);require(!rows.empty(),"Empty manifest");std::set<std::string> inputs,samples;for(auto& r:rows){
    for(auto& k:ManifestFields)require(r.contains(k),"Manifest column missing: "+k);
    auto sample=token(str(r,"sample")),kind=str(r,"kind"),beam=str(r,"beam");require(kind=="led"||kind=="electron"||kind=="muon","Unknown kind");require(beam=="ps"||beam=="sps","Unknown beam");
    double group=num(r,"led_group");require(std::isfinite(group)&&group==int(group)&&(kind=="led"?(group>=0&&group<14):group==-1),"Invalid LED group");
    fs::path raw=str(r,"raw");if(raw.is_relative())raw=p.parent_path()/raw;require(fs::is_regular_file(raw),"Missing raw "+raw.string());r["raw"]=fs::canonical(raw).string();
    require(inputs.insert(str(r,"raw")).second&&samples.insert(sample).second,"Duplicate raw or sample");
}return rows;}
struct Lock {
    int fd=-1;
    explicit Lock(const fs::path& p,bool blocking=false,bool shared=false){fd=::open(p.c_str(),O_CREAT|O_RDWR|O_CLOEXEC,0664);require(fd>=0,"Cannot open lock "+p.string());if(flock(fd,(shared?LOCK_SH:LOCK_EX)|(blocking?0:LOCK_NB))!=0){int error=errno;::close(fd);fd=-1;throw std::runtime_error("Cannot lock "+p.string()+": "+std::strerror(error));}}
    Lock(const Lock&)=delete;Lock& operator=(const Lock&)=delete;
    ~Lock(){if(fd>=0){flock(fd,LOCK_UN);::close(fd);}}
};
inline std::string quote(const std::string& s){std::string out="'";for(char c:s)out+=c=='\''?"'\"'\"'":std::string(1,c);return out+'\'';}
inline std::string commandLine(const std::vector<std::string>& argv){std::string s;for(auto& a:argv){if(!s.empty())s+=' ';s+=quote(a);}return s;}
inline std::string execute(const std::vector<std::string>& args,const fs::path& log={}){
    require(!args.empty(),"Empty command");int output=-1,pipefd[2]={-1,-1};
    if(!log.empty()){output=::open(log.c_str(),O_CREAT|O_WRONLY|O_TRUNC,0664);require(output>=0,"Cannot open command log");}else require(::pipe(pipefd)==0,"Cannot create output pipe");
    std::vector<char*> argv;for(auto& a:args)argv.push_back(const_cast<char*>(a.c_str()));argv.push_back(nullptr);
    pid_t pid=fork();require(pid>=0,"Cannot fork process");if(pid==0){if(output>=0){dup2(output,STDOUT_FILENO);dup2(output,STDERR_FILENO);close(output);}else{close(pipefd[0]);dup2(pipefd[1],STDOUT_FILENO);close(pipefd[1]);}execvp(argv[0],argv.data());perror(argv[0]);_exit(127);}
    std::string text;if(output>=0)close(output);else{close(pipefd[1]);char buffer[8192];ssize_t n;while((n=read(pipefd[0],buffer,sizeof(buffer)))>0)text.append(buffer,n);close(pipefd[0]);}
    int status=0;while(waitpid(pid,&status,0)<0){if(errno!=EINTR)throw std::runtime_error("Cannot wait for process");}
    require(WIFEXITED(status)&&WEXITSTATUS(status)==0,"Command failed: "+commandLine(args)+(log.empty()?"":"; see "+log.string()));return text;
}
struct Context {
    fs::path here,project,result,log;std::vector<std::string> command;
    explicit Context(int argc,char** argv){here=fs::canonical("/proc/self/exe").parent_path();project=std::getenv("SIPM_GAIN_PROJECT")?fs::canonical(std::getenv("SIPM_GAIN_PROJECT")):here.parent_path().parent_path();result=project/"analysis/result/sipm_gain";log=project/"analysis/log/sipm_gain";for(int i=0;i<argc;++i)command.emplace_back(argv[i]);}
    Json sources(bool collectorOnly=false)const{Json j=Json::object();for(auto& e:fs::recursive_directory_iterator(here)){auto p=e.path();auto ext=p.extension().string();if(!e.is_regular_file())continue;if(collectorOnly?(p.filename()=="Collect.cc"||p.filename()=="LEDGroups.hh"):(ext==".cc"||ext==".hh"||ext==".sh"||p.filename()=="Makefile"))j[(fs::path("analysis/sipm_gain")/fs::relative(p,here)).string()]=sha256(p);}auto geometry=project/"analysis/channel_response/temperature_calibration/TemperatureGeometry.hh";j[fs::relative(geometry,project).string()]=sha256(geometry);return j;}
};
struct Options {
    std::string action;std::map<std::string,std::string> values;std::set<std::string> flags;
    std::string get(const std::string& k,const std::string& fallback="")const{auto it=values.find(k);return it==values.end()?fallback:it->second;}
    bool has(const std::string& k)const{return flags.count(k)||values.count(k);}
    double real(const std::string& k,double fallback)const{if(!values.count(k))return fallback;size_t used=0;double x=std::stod(values.at(k),&used);require(used==values.at(k).size()&&std::isfinite(x),"Invalid numeric --"+k);return x;}
    long long integer(const std::string& k,long long fallback)const{double x=real(k,double(fallback));require(x>=-9e15&&x<=9e15&&std::floor(x)==x,"Invalid integer --"+k);return static_cast<long long>(x);}
    std::string required(const std::string& k)const{auto s=get(k);require(!s.empty(),"Required --"+k);return s;}
};
inline std::unique_ptr<TFile> openRoot(const fs::path& p,const char* mode="READ"){
    std::unique_ptr<TFile> f(TFile::Open(p.c_str(),mode));require(f&&!f->IsZombie()&&!f->TestBit(TFile::kRecovered),"Bad ROOT file "+p.string());return f;
}
inline TTree* tree(TFile& f,const char* name){auto* t=dynamic_cast<TTree*>(f.Get(name));require(t,"Missing tree "+std::string(name));return t;}
inline double leaf(TTree& t,const char* name){auto* p=t.GetLeaf(name);require(p,"Missing leaf "+std::string(name));return p->GetValue();}
inline TGraph graph(const Vec& x,const Vec& y,const std::string& name){require(x.size()==y.size(),"Graph dimensions differ");TGraph g(x.size(),x.data(),y.data());g.SetName(name.c_str());return g;}
void compare(const Options& options,const Context& context);
}
