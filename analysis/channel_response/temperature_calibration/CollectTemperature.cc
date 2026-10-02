// Read-only raw/calibrated/track joins, based on CollectMuon and CalibrationResidual.
#include "../calibration_residual/RootInput.hh"
#include "../calibration_residual/ResidualStats.hh"
#include "TemperatureGeometry.hh"
#include <TROOT.h>
#include <climits>
#include <iomanip>

void require(bool b,const std::string&s){if(!b)throw std::runtime_error(s);}
void read(TTree*t,Long64_t e){require(t->GetEntry(e)>0,"Unreadable tree entry");}
struct Event {
    int run=0,time=0,id=0;
    void attach(TTree*t,bool cal){bind(t,"Run_Num",&run);bind(t,"Event_Time",&time);bind(t,cal?"Event_Num":"TriggerID",&id);}
    EventKey key()const{return {run,time,id};}
};
struct Ped {Moments t,h,l;long long missing=0;};
int main(int argc,char**argv)try{
    require(argc==12,"CollectTemperature kind raw calib track output.root pedestal.tsv counters.tsv channels events temp_width max_angle");
    const std::string kind=argv[1];require(kind=="mip"||kind=="response"||kind=="pedestal","Invalid kind");
    const auto cells=std::string(argv[8])=="all"?std::set<int>{}:integerList(argv[8],false);
    const auto limit=std::stoll(argv[9]);const double tw=std::stod(argv[10]),angle=std::stod(argv[11]);
    gROOT->SetBatch(true);
    auto rf=openFile(argv[2]);require(!rf->TestBit(TFile::kRecovered),"Recovered raw file");
    auto*r=tree(*rf,"Raw_Hit");Event re;re.attach(r,false);
    FilePtr cf,kf;TTree*c=nullptr,*k=nullptr;Event ce;
    std::map<EventKey,Long64_t> lookup;bool ordered=false;
    std::map<std::pair<int,int>,Long64_t> tracks;
    if(kind!="pedestal"){
        cf=openFile(argv[3]);require(!cf->TestBit(TFile::kRecovered),"Recovered calibrated file");c=tree(*cf,"Calib_Hit");ce.attach(c,true);
        std::vector<EventKey> keys;for(Long64_t e=0;e<c->GetEntries();++e){read(c,e);keys.push_back(ce.key());}
        ordered=r->GetEntries()==c->GetEntries();
        if(ordered)for(Long64_t e=0;e<r->GetEntries();++e){read(r,e);if(re.key()!=keys[e]){ordered=false;break;}}
        if(!ordered){
            for(size_t e=0;e<keys.size();++e)require(lookup.emplace(keys[e],e).second,"Duplicate calibration key");
            std::set<EventKey> seen;
            for(Long64_t e=0;e<r->GetEntries();++e){read(r,e);require(seen.insert(re.key()).second&&lookup.count(re.key()),"Duplicate/unmatched raw key");}
            require(seen.size()==lookup.size(),"Unmatched calibrated events");
        }
    }
    int kt=0,ki=0;std::vector<int>*kc=nullptr;std::vector<double>*pars=nullptr;
    if(kind=="mip"){
        kf=openFile(argv[4]);require(!kf->TestBit(TFile::kRecovered),"Recovered track file");k=tree(*kf,"T_Event");
        bind(k,"Event_Time",&kt);bind(k,"TriggerID",&ki);
        for(Long64_t e=0;e<k->GetEntries();++e){read(k,e);require(tracks.emplace(std::make_pair(kt,ki),e).second,"Duplicate track key");}
        // Track files omit Run_Num: require unique (time,trigger) in this input raw file.
        std::set<std::pair<int,int>> seen;
        for(Long64_t e=0;e<r->GetEntries();++e){read(r,e);require(seen.emplace(re.time,re.id).second,"Ambiguous raw key for track join");}
        bind(k,"hitCellnew",&kc);bind(k,"trackFitPars",&pars);
    }
    std::vector<int>*ids=nullptr,*tags=nullptr,*ci=nullptr;
    std::vector<double>*hg=nullptr,*lg=nullptr,*eh=nullptr,*el=nullptr,*temp=nullptr,*xx=nullptr,*yy=nullptr,*zz=nullptr;
    std::vector<std::vector<double>>*sensors=nullptr;
    bind(r,"CellID",&ids);bind(r,"HitTag",&tags);bind(r,"HG_Charge",&hg);bind(r,"LG_Charge",&lg);
    bind(r,"Temperature",&sensors);
    if(c){bind(c,"CellID",&ci);bind(c,"Hit_HG_Energy",&eh);bind(c,"Hit_LG_Energy",&el);bind(c,"NewTemperature",&temp);
        if(k){bind(c,"Hit_X",&xx);bind(c,"Hit_Y",&yy);bind(c,"Hit_Z",&zz);}}
    cacheActiveBranches(r);if(c)cacheActiveBranches(c);if(k)cacheActiveBranches(k);
    auto out=openFile(argv[5],"CREATE");TTree hits("Hits","Joined hits; no cut on corrected energy");
    int run=0,cellid=0,memory=0,temperature_valid=0;double raw_hg=0,raw_lg=0,saved_hg=0,saved_lg=0,temperature=0,sensor_temperature=0,path_factor=1;
    for(auto v:std::vector<std::pair<const char*,int*>>{{"run",&run},{"cellid",&cellid},{"memory",&memory},{"temperature_valid",&temperature_valid}})hits.Branch(v.first,v.second);
    for(auto v:std::vector<std::pair<const char*,double*>>{{"raw_hg",&raw_hg},{"raw_lg",&raw_lg},{"saved_hg",&saved_hg},{"saved_lg",&saved_lg},{"temperature",&temperature},{"sensor_temperature",&sensor_temperature},{"path_factor",&path_factor}})hits.Branch(v.first,v.second);
    using PK=std::tuple<int,int,int,int>; // run, channel, memory, temperature bin; INT_MIN = unknown
    std::map<PK,Ped> ped;
    long long scanned=0,selectedTracks=0,missingSensors=0,nonfinite=0,signal=0;
    Long64_t total=r->GetEntries(),take=limit?std::min<Long64_t>(limit,total):total;
    for(Long64_t j=0;j<take;++j){
        auto e=take==total?j:static_cast<Long64_t>((j+.5L)*total/take);read(r,e);++scanned;
        require(ids->size()==hg->size()&&ids->size()==lg->size()&&ids->size()==tags->size(),"Raw vector mismatch");
        std::map<int,size_t> raw;
        for(size_t i=0;i<ids->size();++i)require(raw.emplace(ids->at(i),i).second,"Duplicate full raw CellID");
        auto selected=[&](int id){return id/100000<30&&connected(id)&&(cells.empty()||cells.count(id));};
        auto sensorT=[&](int id)->double{
            int layer=id/100000;
            if(layer>=int(sensors->size())||sensors->at(layer).size()!=16)return NAN;
            const auto&ts=sensors->at(layer);
            if(!std::all_of(ts.begin(),ts.end(),[](double t){return std::isfinite(t)&&t>0&&t<60;}))return NAN;
            return tempReconstruction(layer,EBUdecode(layer,id%100000/10000,id%100),ts);
        };
        for(size_t i=0;i<ids->size();++i){
            int id=physicalID(ids->at(i));if(!selected(id)||tags->at(i)!=0)continue;
            double t=sensorT(id),h=hg->at(i),l=lg->at(i);
            if(!std::isfinite(h)||!std::isfinite(l)){++nonfinite;continue;}
            auto&p=ped[{re.run,id,ids->at(i)%10000/100,std::isfinite(t)?int(std::floor(t/tw)):INT_MIN}];
            p.h.add(h);p.l.add(l);if(std::isfinite(t))p.t.add(t);else ++p.missing;
        }
        if(!c)continue;
        std::set<int> ontrack;double sx=0,sy=0,x0=0,y0=0;
        if(k){
            auto it=tracks.find({re.time,re.id});if(it==tracks.end())continue;read(k,it->second);
            if(pars->size()<10||!std::all_of(pars->begin(),pars->end(),[](double v){return std::isfinite(v);})||pars->at(9)<=0||pars->at(8)<0)continue;
            sx=pars->at(0);sy=pars->at(4);x0=pars->at(2);y0=pars->at(6);
            if(std::abs(std::atan(sx))>angle||std::abs(std::atan(sy))>angle||std::sqrt(pars->at(8)/pars->at(9))>3)continue;
            // TrackFit3D uses -1 for a missing layer (as in CollectMuon).
            std::set<int> planes[2];for(int full:*kc){if(full<0)continue;int id=physicalID(full);if(id/100000<30&&connected(id)){ontrack.insert(id);planes[id/100000%2].insert(id/100000);}}
            if(planes[0].size()<7||planes[1].size()<7)continue;
            ++selectedTracks;
        }
        read(c,ordered?e:lookup.at(re.key()));require(re.key()==ce.key(),"Event join changed");
        size_t n=ci->size();require(eh->size()==n&&el->size()==n&&temp->size()==n,"Calibrated vector mismatch");
        if(k)require(xx->size()==n&&yy->size()==n&&zz->size()==n,"Position vector mismatch");
        std::map<int,int> counts;std::set<int> fulls;
        for(int full:*ci){require(fulls.insert(full).second,"Duplicate calibrated CellID");++counts[physicalID(full)];}
        for(size_t i=0;i<n;++i){
            int full=ci->at(i);cellid=physicalID(full);if(!selected(cellid))continue;
            auto it=raw.find(full);require(it!=raw.end()&&tags->at(it->second)==1,"Calibrated hit missing raw HitTag=1");
            auto ri=it->second;raw_hg=hg->at(ri);raw_lg=lg->at(ri);
            if(k){
                bool even=cellid/100000%2==0;
                double longitudinal=even?x0+sx*zz->at(i)-xx->at(i):y0+sy*zz->at(i)-yy->at(i);
                double transverse=even?y0+sy*zz->at(i)-yy->at(i):x0+sx*zz->at(i)-xx->at(i);
                if(!ontrack.count(cellid)||counts[cellid]!=1||raw_hg>=2600||!std::isfinite(longitudinal)||!std::isfinite(transverse)||std::abs(longitudinal)+std::abs(even?sx:sy)>=21.5||std::abs(transverse)+std::abs(even?sy:sx)>=2.25)continue;
            }
            run=re.run;memory=full%10000/100;saved_hg=eh->at(i);saved_lg=el->at(i);temperature=temp->at(i);sensor_temperature=sensorT(cellid);
            temperature_valid=std::isfinite(sensor_temperature)&&std::isfinite(temperature)&&temperature>0&&temperature<60;
            if(!temperature_valid)++missingSensors;
            path_factor=1/std::sqrt(1+sx*sx+sy*sy);hits.Fill();++signal;
        }
        if(scanned%10000==0)std::cout<<"scanned "<<scanned<<" / "<<take<<std::endl;
    }
    std::ofstream p(argv[6]);require(bool(p),"Cannot write pedestal table");p<<std::setprecision(17);
    p<<"run\tcellid\tmemory\ttemp_bin\ttemperature\tn\thg_mean\thg_sem\tlg_mean\tlg_sem\tmissing_temperature\n";
    for(const auto&[key,v]:ped){auto[rr,id,mem,tb]=key;p<<rr<<'\t'<<id<<'\t'<<mem<<'\t'<<tb<<'\t'<<(v.t.n?v.t.mean:NAN)<<'\t'<<v.h.n<<'\t'<<v.h.mean<<'\t'<<v.h.sem()<<'\t'<<v.l.mean<<'\t'<<v.l.sem()<<'\t'<<v.missing<<'\n';}
    std::ofstream a(argv[7]);require(bool(a),"Cannot write counters");
    a<<"events_total\tevents_scanned\tselected_tracks\tsignal_hits\tmissing_sensor_hits\tnonfinite_ped_hits\tmatching\n"<<total<<'\t'<<scanned<<'\t'<<selectedTracks<<'\t'<<signal<<'\t'<<missingSensors<<'\t'<<nonfinite<<'\t'<<(ordered?"verified_entry_order":"unique_event_key_or_raw_only")<<'\n';
    out->cd();hits.Write();require(!out->TestBit(TFile::kWriteError)&&bool(a)&&bool(p),"Write failure");
    return 0;
}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
