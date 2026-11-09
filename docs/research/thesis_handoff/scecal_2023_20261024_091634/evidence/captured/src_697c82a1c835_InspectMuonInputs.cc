#include "MuonMIPCommon.hh"
#include <TROOT.h>
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace mmr;
struct InputRow {std::vector<std::string> fields;bool valid=true;};
int main(int argc,char**argv)try{
    require(argc==5,"Usage: InspectMuonInputs requested.tsv config.tsv manifest.tsv file_metadata.tsv");
    Config cfg(argv[2]);const bool allow=cfg.integer("allow-recovered",0);gROOT->SetBatch(true);
    std::ifstream in(argv[1]);require(bool(in),"Cannot read requested manifest");
    std::ofstream metadata(argv[4]),out(argv[3]);require(bool(metadata)&&bool(out),"Cannot write input audit");
    metadata<<"sample\trun\tstage\tpath\tentries\trecovered\tbytes\tmtime_unix_seconds\tstatus\n";
    std::vector<InputRow> rows;std::map<int,bool> runOK;
    for(std::string line;std::getline(in,line);){if(line.empty()||line[0]=='#')continue;InputRow row{split(line),true};require(row.fields.size()==6,"Malformed requested manifest");
        int sample=std::stoi(row.fields[0]),run=std::stoi(row.fields[1]);require(sample==0||sample==1,"Invalid sample");if(!runOK.count(run))runOK[run]=true;
        const std::array<const char*,3> stages={"decode","calib","truth"},names={"Raw_Hit","Calib_Hit","MC_Truth"};
        for(int stage=0;stage<(sample?3:2);++stage){const auto&path=row.fields[stage==2?5:stage+2];long long entries=0,bytes=0,mtime=0;bool recovered=false;std::string status="ok";
            std::error_code ec;if(!fs::is_regular_file(path,ec)){status="missing";}else{
                bytes=fs::file_size(path,ec);auto ft=fs::last_write_time(path,ec);if(!ec){auto st=std::chrono::time_point_cast<std::chrono::seconds>(ft-fs::file_time_type::clock::now()+std::chrono::system_clock::now());mtime=st.time_since_epoch().count();}
                File f(TFile::Open(path.c_str(),"READ"));if(!f||f->IsZombie())status="unreadable";else{recovered=f->TestBit(TFile::kRecovered);auto*t=dynamic_cast<TTree*>(f->Get(names[stage]));if(!t)status="missing_tree";else{entries=t->GetEntries();if(entries<=0)status="empty_tree";}if(recovered&&!allow)status="recovered_excluded";}
            }
            metadata<<(sample?"mc":"data")<<'\t'<<run<<'\t'<<stages[stage]<<'\t'<<path<<'\t'<<entries<<'\t'<<recovered<<'\t'<<bytes<<'\t'<<mtime<<'\t'<<status<<'\n';
            if(status!="ok"){row.valid=false;runOK[run]=false;std::cerr<<"Exclude run "<<run<<": "<<status<<" "<<path<<'\n';}
        }rows.push_back(std::move(row));
    }
    long long n[2]={};for(const auto&r:rows)if(r.valid&&runOK.at(std::stoi(r.fields[1]))){for(size_t i=0;i<r.fields.size();++i)out<<(i?"\t":"")<<r.fields[i];out<<'\n';++n[std::stoi(r.fields[0])];}
    auto sample=cfg.text("sample","both");require((sample=="mc"||n[0]>0)&&(sample=="data"||n[1]>0),"No usable input run for a requested sample; inspect file_metadata.tsv");
    require(bool(out)&&bool(metadata),"Input audit write failed");std::cout<<"Usable inputs: data "<<n[0]<<", MC "<<n[1]<<"; recovered="<<allow<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
