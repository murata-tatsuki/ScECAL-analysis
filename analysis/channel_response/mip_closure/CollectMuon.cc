// Read-only joins. All selections are independent of the MIP policy under test.
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using Key=std::pair<int,int>; // Event_Time, trigger; scoped to ONE input file/run
using File=std::unique_ptr<TFile>;
void require(bool ok,const std::string& s) { if(!ok) throw std::runtime_error(s); }
File open(const std::string& p,bool recovered) {
    File f(TFile::Open(p.c_str(),"READ"));
    require(f&&!f->IsZombie(),"Cannot open "+p);
    require(recovered||!f->TestBit(TFile::kRecovered),"Recovered ROOT file: "+p);
    return f;
}
TTree* tree(TFile& f,const char* n) {
    auto*t=dynamic_cast<TTree*>(f.Get(n)); require(t,std::string("Missing ")+n+" in "+f.GetName());
    t->SetBranchStatus("*",false); return t;
}
template<class T> void bind(TTree*t,const char*n,T* x) {
    require(t->GetBranch(n),std::string("Missing branch ")+n);
    t->SetBranchStatus(n,true);require(t->SetBranchAddress(n,x)>=0,std::string("Type mismatch ")+n);
}
void read(TTree*t,Long64_t i) {require(t->GetEntry(i)>0,"Cannot read entry "+std::to_string(i)+" in "+t->GetName());}
int physical(int c) {
    require(c>=0&&c/100000<32&&c%100000/10000<6&&c%100<36,"Invalid CellID");
    return c/100000*100000+c%100000/10000*10000+c%100;
}
bool connected(int c) {return c/100000<30&&(c%100000/10000!=5||c%100<30);}
std::map<Key,Long64_t> index(TTree*t,const char* event) {
    int ti=0,id=0;bind(t,"Event_Time",&ti);bind(t,event,&id);
    std::map<Key,Long64_t> m;
    for(Long64_t i=0;i<t->GetEntries();++i) {read(t,i);require(m.emplace(Key{ti,id},i).second,"Duplicate event key in "+std::string(t->GetName()));}
    t->ResetBranchAddresses();t->SetBranchStatus("*",false);return m;
}
struct Raw {
    int id=0,time=0; std::vector<int>* cells=nullptr,*tags=nullptr;
    std::vector<double>* hg=nullptr;
    void attach(TTree*t){bind(t,"TriggerID",&id);bind(t,"Event_Time",&time);bind(t,"CellID",&cells);bind(t,"HitTag",&tags);bind(t,"HG_Charge",&hg);}
};
struct Cal {
    int id=0,time=0;std::vector<int>*cells=nullptr;
    std::vector<double>* energy=nullptr,*hg=nullptr,*temp=nullptr,*x=nullptr,*y=nullptr,*z=nullptr;
    void attach(TTree*t) {bind(t,"Event_Num",&id);bind(t,"Event_Time",&time);bind(t,"CellID",&cells);bind(t,"Hit_Energy",&energy);bind(t,"Hit_HG_Energy",&hg);bind(t,"NewTemperature",&temp);bind(t,"Hit_X",&x);bind(t,"Hit_Y",&y);bind(t,"Hit_Z",&z);}
};
struct Track {
    int id=0,time=0;double savedLength=0;std::vector<int>*cells=nullptr;
    std::vector<double>*pars=nullptr;
    void attach(TTree*t) {bind(t,"TriggerID",&id);bind(t,"Event_Time",&time);bind(t,"hitCellnew",&cells);bind(t,"trackFitPars",&pars);bind(t,"realLength",&savedLength);}
};
struct Truth {
    int pdg=0;double px=0,py=0,pz=0;
    std::vector<int>* cells=nullptr;std::vector<double>* energy=nullptr,*direction=nullptr;
    void attach(TTree*t) {bind(t,"CellID",&cells);bind(t,"Hit_Energy",&energy);bind(t,"PrimaryDirection",&direction);bind(t,"PrimaryPDG",&pdg);bind(t,"PrimaryPosX",&px);bind(t,"PrimaryPosY",&py);bind(t,"PrimaryPosZ",&pz);}
};

int main(int argc,char**argv) try {
    // manifest, output ROOT, validation TSV, max tracks/file, angle [rad], rms [mm], allow recovered
    require(argc==8,"Usage: CollectMuon manifest.tsv flat.root validation.tsv max_tracks max_angle max_rms allow_recovered");
    const long long limit=std::stoll(argv[4]);const double maxAngle=std::stod(argv[5]),maxRms=std::stod(argv[6]);
    const bool allow=std::stoi(argv[7]);gROOT->SetBatch(true);
    std::ifstream manifest(argv[1]);require(bool(manifest),"Cannot read manifest");
    File output(TFile::Open(argv[2],"CREATE"));require(output&&!output->IsZombie(),"Output already exists or is not writable");
    std::ofstream validation(argv[3]);require(bool(validation),"Cannot write validation");
    validation<<"sample\trun\traw_entries\tcalib_entries\ttrack_entries\ttruth_entries\ttracks_scanned\ttracks_selected\tmissing_raw_events\tmissing_calib_events\tmissing_raw_hits\tambiguous_hits\tselected_hits\trecovered\n";
    int sample=0,run=0,event=0,event_time=0,cellid=0,full_cellid=0,mip_selected=0,raw_ok=0;
    Long64_t event_index=0;
    double raw_hg=NAN,saved_hg=NAN,saved_energy=NAN,temperature=NAN,x=NAN,y=NAN,z=NAN;
    double sx=NAN,sy=NAN,x0=NAN,y0=NAN,length=NAN,saved_length=NAN,track_rms=NAN;
    double longitudinal=NAN,transverse=NAN,truth_edep=NAN;
    TTree hits("Hits","All calibrated hits in selected muon events; no MIP-dependent cuts");
    for(auto item:std::vector<std::pair<const char*,int*>>{{"sample",&sample},{"run",&run},{"event",&event},{"event_time",&event_time},{"cellid",&cellid},{"full_cellid",&full_cellid},{"mip_selected",&mip_selected},{"raw_ok",&raw_ok}}) hits.Branch(item.first,item.second);
    hits.Branch("event_index",&event_index);
    for(auto item:std::vector<std::pair<const char*,double*>>{{"raw_hg",&raw_hg},{"saved_hg",&saved_hg},{"saved_energy",&saved_energy},{"temperature",&temperature},{"x",&x},{"y",&y},{"z",&z},{"sx",&sx},{"sy",&sy},{"length",&length},{"longitudinal",&longitudinal},{"transverse",&transverse},{"truth_edep",&truth_edep}}) hits.Branch(item.first,item.second);
    TTree events("Events","Track parameters and join accounting");
    int n_calib=0,n_raw_missing=0,n_ambiguous=0,n_mip=0;
    events.Branch("event_index",&event_index);
    for(auto item:std::vector<std::pair<const char*,int*>>{{"sample",&sample},{"run",&run},{"event",&event},{"event_time",&event_time},{"n_calib",&n_calib},{"n_raw_missing",&n_raw_missing},{"n_ambiguous",&n_ambiguous},{"n_mip",&n_mip}}) events.Branch(item.first,item.second);
    for(auto item:std::vector<std::pair<const char*,double*>>{{"sx",&sx},{"sy",&sy},{"x0",&x0},{"y0",&y0},{"length",&length},{"saved_length",&saved_length},{"track_rms",&track_rms}}) events.Branch(item.first,item.second);
    // Audit primary direction for every truth entry, independently of reconstructed selection.
    TTree beam("Beam","All generated primaries"); double px=0,py=0,pz=0,dx=0,dy=0,dz=0;
    beam.Branch("run",&run);
    for(auto item:std::vector<std::pair<const char*,double*>>{{"px",&px},{"py",&py},{"pz",&pz},{"dx",&dx},{"dy",&dy},{"dz",&dz}}) beam.Branch(item.first,item.second);
    for(std::string line;std::getline(manifest,line);) {
        if(line.empty()||line[0]=='#')continue;
        std::vector<std::string> fields;std::istringstream is(line);
        for(std::string p;std::getline(is,p,'\t');)fields.push_back(p);
        require(fields.size()==6,"Expected sample,run,raw,calib,track,truth columns");
        sample=std::stoi(fields[0]);run=std::stoi(fields[1]);
        auto fr=open(fields[2],allow),fc=open(fields[3],allow),fk=open(fields[4],allow);
        auto*tr=tree(*fr,"Raw_Hit");auto*tc=tree(*fc,"Calib_Hit");auto*tk=tree(*fk,"T_Event");
        auto ri=index(tr,"TriggerID"),ci=index(tc,"Event_Num");
        Raw r;Cal c;Track k;Truth v;r.attach(tr);c.attach(tc);k.attach(tk);
        File ft;TTree*tt=nullptr;
        if(sample==1) {
            ft=open(fields[5],allow);tt=tree(*ft,"MC_Truth");v.attach(tt);
            // Only read primary branches for this full-file audit.
            tt->SetBranchStatus("CellID",false);tt->SetBranchStatus("Hit_Energy",false);
            for(Long64_t i=0;i<tt->GetEntries();++i) {
                read(tt,i);require(v.direction->size()==3&&std::abs(v.pdg)==13,"Not a muon truth sample");
                px=v.px;py=v.py;pz=v.pz;dx=v.direction->at(0);dy=v.direction->at(1);dz=v.direction->at(2);beam.Fill();
            }
            tt->SetBranchStatus("CellID",true);tt->SetBranchStatus("Hit_Energy",true);
        }
        Long64_t scanned=0,selected=0,missR=0,missC=0,missH=0,ambH=0,selH=0;
        const Long64_t total=tk->GetEntries(),n=limit>0?std::min<Long64_t>(limit,total):total;
        std::set<Key> trackKeys;
        for(Long64_t j=0;j<n;++j) {
            const Long64_t entry=limit>0?static_cast<Long64_t>((j+.5L)*total/n):j;
            read(tk,entry);++scanned;Key key{k.time,k.id};
            require(trackKeys.insert(key).second,"Duplicate sampled track event key");
            if(k.pars->size()<10)continue;
            const auto&p=*k.pars;
            if(!std::all_of(p.begin(),p.end(),[](double a){return std::isfinite(a);})||p[9]<=0||p[8]<0)continue;
            sx=p[0];sy=p[4];x0=p[2];y0=p[6];track_rms=std::sqrt(p[8]/p[9]);
            if(std::abs(std::atan(sx))>maxAngle||std::abs(std::atan(sy))>maxAngle||track_rms>maxRms)continue;
            std::set<int> planes[2],ontrack;
            for(int id:*k.cells)if(id>=0&&connected(physical(id))) {id=physical(id);ontrack.insert(id);planes[id/100000%2].insert(id/100000);}
            if(planes[0].size()<7||planes[1].size()<7)continue;
            if(!ri.count(key)){++missR;continue;}if(!ci.count(key)){++missC;continue;}
            read(tr,ri.at(key));read(tc,ci.at(key));
            require(r.id==k.id&&r.time==k.time&&c.id==k.id&&c.time==k.time,"Event join mismatch");
            require(r.cells->size()==r.hg->size()&&r.cells->size()==r.tags->size(),"Raw vector size mismatch");
            const size_t nc=c.cells->size();
            require(nc==c.hg->size()&&nc==c.energy->size()&&nc==c.temp->size()&&nc==c.x->size()&&nc==c.y->size()&&nc==c.z->size(),"Calib vector size mismatch");
            std::map<int,double> truth;
            if(tt) {
                // MCDigi.cxx explicitly sets TriggerID=i_entry. Never use track entry order.
                require(k.id>=0&&k.id<tt->GetEntries()&&r.time==r.id,"MC truth mapping does not match MCDigi convention");
                read(tt,k.id);require(v.cells->size()==v.energy->size(),"Truth vector size mismatch");
                for(size_t i=0;i<v.cells->size();++i)truth[physical(v.cells->at(i))]+=v.energy->at(i);
            }
            std::map<int,std::vector<size_t>> rawIndex,calCounts;
            for(size_t i=0;i<r.cells->size();++i)if(r.tags->at(i))rawIndex[r.cells->at(i)].push_back(i);
            for(size_t i=0;i<nc;++i)calCounts[physical(c.cells->at(i))].push_back(i);
            ++selected;event=k.id;event_time=k.time;length=2*std::sqrt(1+sx*sx+sy*sy);saved_length=k.savedLength;
            n_calib=n_raw_missing=n_ambiguous=n_mip=0;
            for(size_t i=0;i<nc;++i) {
                full_cellid=c.cells->at(i);cellid=physical(full_cellid);if(!connected(cellid))continue;
                ++n_calib;raw_ok=0;raw_hg=NAN;
                auto it=rawIndex.find(full_cellid);
                if(it==rawIndex.end()){++n_raw_missing;++missH;}
                else if(it->second.size()!=1){++n_ambiguous;++ambH;}
                else {raw_ok=1;raw_hg=r.hg->at(it->second[0]);}
                saved_hg=c.hg->at(i);saved_energy=c.energy->at(i);temperature=c.temp->at(i);
                x=c.x->at(i);y=c.y->at(i);z=c.z->at(i);
                const bool even=cellid/100000%2==0;
                longitudinal=even?x0+sx*z-x:y0+sy*z-y;
                transverse=even?y0+sy*z-y:x0+sx*z-x;
                // Physical active strip 45 x 5 x 2 mm. Require the entire segment inside.
                mip_selected=raw_ok&&raw_hg<2600&&ontrack.count(cellid)&&calCounts[cellid].size()==1
                    &&std::abs(longitudinal)+std::abs(even?sx:sy)<21.5
                    &&std::abs(transverse)+std::abs(even?sy:sx)<2.25;
                truth_edep=tt&&truth.count(cellid)?truth.at(cellid):NAN;
                if(mip_selected){++n_mip;++selH;}
                hits.Fill();
            }
            events.Fill();++event_index;
        }
        const bool recovered=fr->TestBit(TFile::kRecovered)||fc->TestBit(TFile::kRecovered)||fk->TestBit(TFile::kRecovered)||(ft&&ft->TestBit(TFile::kRecovered));
        validation<<sample<<'\t'<<run<<'\t'<<tr->GetEntries()<<'\t'<<tc->GetEntries()<<'\t'<<total<<'\t'<<(tt?tt->GetEntries():0)<<'\t'<<scanned<<'\t'<<selected<<'\t'<<missR<<'\t'<<missC<<'\t'<<missH<<'\t'<<ambH<<'\t'<<selH<<'\t'<<recovered<<'\n';
        std::cout<<(sample?"MC":"data")<<" run "<<run<<": "<<selected<<" selected tracks, "<<selH<<" MIP hits"<<std::endl;
        tr->ResetBranchAddresses();tc->ResetBranchAddresses();tk->ResetBranchAddresses();if(tt)tt->ResetBranchAddresses();
        delete r.cells;delete r.tags;delete r.hg;
        delete c.cells;delete c.energy;delete c.hg;delete c.temp;delete c.x;delete c.y;delete c.z;
        delete k.cells;delete k.pars;delete v.cells;delete v.energy;delete v.direction;
    }
    output->cd();hits.Write();events.Write();beam.Write();require(bool(validation),"Validation write failed");
    return 0;
} catch(const std::exception&e) {std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
