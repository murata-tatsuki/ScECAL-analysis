#pragma once
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Link the established detector coordinate map; never infer geometry from hits.
double* EBUdecode(int,int,int);
namespace mmr {
namespace fs = std::filesystem;
inline void require(bool ok,const std::string& msg) {if(!ok) throw std::runtime_error(msg);}
using File = std::unique_ptr<TFile>;
inline File openFile(const std::string& path,bool allow=false) {
    File f(TFile::Open(path.c_str(),"READ"));
    require(f&&!f->IsZombie(),"Cannot open "+path);
    require(allow||!f->TestBit(TFile::kRecovered),"Recovered ROOT file: "+path+"; use --allow-recovered for diagnostics");
    return f;
}
inline TTree* tree(TFile& f,const char* name) {
    auto* t=dynamic_cast<TTree*>(f.Get(name));require(t,"Missing tree "+std::string(name)+" in "+f.GetName());
    t->SetBranchStatus("*",false);return t;
}
template<class T> inline void bind(TTree* t,const char* name,T* ptr) {
    require(t->GetBranch(name),"Missing branch "+std::string(name));
    t->SetBranchStatus(name,true);require(t->SetBranchAddress(name,ptr)>=0,"Type mismatch: "+std::string(name));
}
inline void read(TTree* t,Long64_t i) {require(t->GetEntry(i)>0,"Cannot read "+std::string(t->GetName())+" entry "+std::to_string(i));}
inline int physical(int id) {if(id<0)return -1;return id/100000*100000+id%100000/10000*10000+id%100;}
inline bool connected(int id) {return id>=0&&id/100000<30&&id%100000/10000<6&&id%100<(id%100000/10000==5?30:36);}
inline int group(int id) {return id/100000>=4&&id/100000<=27;}
inline std::array<double,3> geometry(int id) {
    require(connected(id),"Invalid geometry channel");auto* p=::EBUdecode(id/100000,id%100000/10000,id%100);return {p[0],p[1],p[2]};
}
inline std::vector<int> channels() {std::vector<int> out;for(int l=0;l<30;++l)for(int c=0;c<6;++c)for(int ch=0;ch<(c==5?30:36);++ch)out.push_back(l*100000+c*10000+ch);return out;}
inline double tempFactor(int id,double t) {return 1-(t-20)*(group(id)?1.6/135:3.5/230);}
class Config {
    std::map<std::string,std::string> values;
public:
    explicit Config(const std::string& path) {
        std::ifstream f(path);require(bool(f),"Cannot open configuration: "+path);
        for(std::string line;std::getline(f,line);){if(line.empty()||line[0]=='#')continue;auto at=line.find('\t');require(at!=std::string::npos,"Malformed configuration line");require(values.emplace(line.substr(0,at),line.substr(at+1)).second,"Duplicate configuration key");}
    }
    std::string text(const std::string& key,const std::string& fallback="")const {auto i=values.find(key);return i==values.end()?fallback:i->second;}
    double number(const std::string& key,double fallback=0)const {auto s=text(key);if(s.empty())return fallback;size_t n=0;double x=std::stod(s,&n);require(n==s.size()&&std::isfinite(x),"Invalid number: "+key);return x;}
    long long integer(const std::string& key,long long fallback=0)const {auto s=text(key);if(s.empty())return fallback;size_t n=0;auto x=std::stoll(s,&n);require(n==s.size(),"Invalid integer: "+key);return x;}
};
inline std::vector<std::string> split(const std::string& line,char delimiter='\t') {std::vector<std::string> v;std::istringstream in(line);for(std::string s;std::getline(in,s,delimiter);)v.push_back(s);return v;}
} // namespace mmr
