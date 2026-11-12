#include "../ssa_shower_shape/Common.hh"

// One union of physical CellIDs per sampled event. Never infer lost ADC values.
struct Signal {
  int truthN=0, rawN=0, taggedN=0, calibN=0;
  double truth=0, hg=0, lg=0, calib=0;
};
int main(int argc,char** argv) {
  try {
    require(argc==7,"Extract manifest beam_events strip_events output axis_layers energy_GeV");
    auto cached=openRoot(argv[2]), strips=openRoot(argv[3]);
    auto* beam=getTree(*cached,"events"); auto* strip=getTree(*strips,"events");
    int savedFile=0,savedRun=0,savedEvent=0,savedTime=0,stripFile=0,stripRun=0,stripEvent=0,stripTime=0;
    Long64_t savedEntry=0,stripEntry=0,cursor=0; bool savedValid=false,stripValid=false;
    double axisX=0,axisY=0;
    std::vector<int>* stripID=nullptr; std::vector<double>* stripE=nullptr;
    beam->SetBranchStatus("*",0); strip->SetBranchStatus("*",0);
    bind(beam,"fileIndex",&savedFile); bind(beam,"entry",&savedEntry); bind(beam,"run",&savedRun);
    bind(beam,"event",&savedEvent); bind(beam,"eventTime",&savedTime); bind(beam,"axisValid",&savedValid);
    bind(beam,"axisX",&axisX); bind(beam,"axisY",&axisY);
    bind(strip,"fileIndex",&stripFile); bind(strip,"entry",&stripEntry); bind(strip,"run",&stripRun);
    bind(strip,"event",&stripEvent); bind(strip,"eventTime",&stripTime); bind(strip,"axisValid",&stripValid);
    bind(strip,"cellID",&stripID); bind(strip,"energy",&stripE);
    require(strip->GetEntries()==beam->GetEntries(),"Beam and strip cache lengths differ");
    auto out=openRoot(argv[4],"RECREATE");
    TTree events("events","Matched MC truth, saved digitization, calibration; sampled beam cohort");
    int fileIndex=0,run=0,event=0,eventTime=0,invalidTruth=0,invalidRaw=0,invalidCalib=0;
    Long64_t entry=0; bool axisValid=false;
    std::vector<int> cellID,layer,truthMultiplicity,rawMultiplicity,taggedMultiplicity,calibMultiplicity;
    std::vector<double> x,y,truthEnergy,hgADC,lgADC,calibEnergy;
#define BR(n,t) events.Branch(#n,&n,#n "/" t)
    BR(fileIndex,"I"); BR(entry,"L"); BR(run,"I"); BR(event,"I"); BR(eventTime,"I"); BR(axisValid,"O");
    BR(axisX,"D"); BR(axisY,"D"); BR(invalidTruth,"I"); BR(invalidRaw,"I"); BR(invalidCalib,"I");
#undef BR
#define V(n) events.Branch(#n,&n)
    V(cellID); V(layer); V(x); V(y); V(truthMultiplicity); V(rawMultiplicity); V(taggedMultiplicity); V(calibMultiplicity);
    V(truthEnergy); V(hgADC); V(lgADC); V(calibEnergy);
#undef V
    std::ifstream manifest(argv[1]); require(bool(manifest),"Cannot open manifest"); std::string line;
    while(std::getline(manifest,line)) {
      if(line.empty())continue;
      std::istringstream row(line); std::string truthPath,rawPath,calibPath,indices;
      require(bool(std::getline(row,truthPath,'\t')) && bool(std::getline(row,rawPath,'\t')) &&
              bool(std::getline(row,calibPath,'\t')) && bool(std::getline(row,indices)),"Invalid manifest");
      auto tf=openRoot(truthPath),rf=openRoot(rawPath),cf=openRoot(calibPath);
      auto* tt=getTree(*tf,"MC_Truth"); auto* rt=getTree(*rf,"Raw_Hit"); Reader calib(getTree(*cf,"Calib_Hit"),false);
      require(tt->GetEntries()==rt->GetEntries() && tt->GetEntries()==calib.t->GetEntries(),"Stage entry counts differ");
      Truth primary(tt); std::vector<int>* tids=nullptr,*rids=nullptr,*tags=nullptr;
      std::vector<double>* te=nullptr,*hg=nullptr,*lg=nullptr;
      bind(tt,"CellID",&tids); bind(tt,"Hit_Energy",&te);
      int rawRun=0,rawTime=0,trigger=0,cycle=0; rt->SetBranchStatus("*",0);
      bind(rt,"Run_Num",&rawRun); bind(rt,"Event_Time",&rawTime); bind(rt,"TriggerID",&trigger); bind(rt,"CycleID",&cycle);
      bind(rt,"CellID",&rids); bind(rt,"HitTag",&tags); bind(rt,"HG_Charge",&hg); bind(rt,"LG_Charge",&lg);
      std::replace(indices.begin(),indices.end(),',',' '); std::istringstream selected(indices);
      while(selected>>entry) {
        require(cursor<beam->GetEntries() && beam->GetEntry(cursor)>0 && strip->GetEntry(cursor++)>0,"Cached event count mismatch");
        require(savedFile==fileIndex && savedEntry==entry && stripFile==fileIndex && stripEntry==entry,"Cached manifest identity mismatch");
        primary.read(entry); require(std::abs(primary.energy-std::stod(argv[6])*1000)<1e-6,"Wrong primary energy");
        require(rt->GetEntry(entry)>0,"Cannot read digitized entry"); calib.read(entry);
        require(std::tie(calib.run,calib.event,calib.time)==std::tie(savedRun,savedEvent,savedTime) &&
                std::tie(calib.run,calib.event,calib.time)==std::tie(stripRun,stripEvent,stripTime),"Cached calib identity mismatch");
        require(rawRun==calib.run && trigger==calib.event && rawTime==calib.time && cycle==entry && trigger==entry && rawTime==entry,
                "MC identity differs; truth entry matching is not justified");
        require(tids->size()==te->size() && rids->size()==tags->size() && rids->size()==hg->size() && rids->size()==lg->size(),"Vector length mismatch");
        auto hits=calib.hits(invalidCalib); const auto axis=entranceAxis(hits,std::stoi(argv[5]));
        require(axis.valid==savedValid && stripValid==savedValid,"Cached entrance selection changed");
        if(axis.valid)require(std::abs(axis.x-axisX)<1e-9 && std::abs(axis.y-axisY)<1e-9,"Cached entrance position changed");
        std::map<int,Signal> signals; invalidTruth=invalidRaw=0;
        for(size_t i=0;i<tids->size();++i) {
          require(std::isfinite(te->at(i)),"Nonfinite truth energy");
          if(!connected(tids->at(i))){++invalidTruth;continue;}
          auto& s=signals[physical(tids->at(i))]; ++s.truthN; s.truth+=te->at(i);
        }
        for(size_t i=0;i<rids->size();++i) {
          require(std::isfinite(hg->at(i)) && std::isfinite(lg->at(i)),"Nonfinite ADC");
          if(!connected(rids->at(i))){++invalidRaw;continue;}
          require(tags->at(i)==0 || tags->at(i)==1,"Unexpected HitTag");
          auto& s=signals[physical(rids->at(i))];
          if(!s.rawN){s.hg=hg->at(i);s.lg=lg->at(i);}else{s.hg=std::max(s.hg,hg->at(i));s.lg=std::max(s.lg,lg->at(i));}
          ++s.rawN; s.taggedN+=tags->at(i)!=0;
        }
        for(const auto& hit:hits){auto& s=signals[physical(hit.id)]; ++s.calibN; s.calib+=hit.e;}
        require(stripID->size()==stripE->size(),"Invalid strip cache");
        size_t nCalib=0;for(const auto& kv:signals)nCalib+=kv.second.calibN>0;
        require(nCalib==stripID->size(),"Cached physical strip count changed");
        for(size_t i=0;i<stripID->size();++i) {
          auto it=signals.find(stripID->at(i)); require(it!=signals.end() && it->second.calibN>0 &&
            std::abs(it->second.calib-stripE->at(i))<1e-9,"Cached strip signal changed");
        }
        cellID.clear();layer.clear();x.clear();y.clear();truthEnergy.clear();hgADC.clear();lgADC.clear();calibEnergy.clear();
        truthMultiplicity.clear();rawMultiplicity.clear();taggedMultiplicity.clear();calibMultiplicity.clear();
        for(const auto& kv:signals) {
          const int id=kv.first,l=id/100000; const auto& s=kv.second; auto* pos=EBUdecode(l,id%100000/10000,id%100);
          cellID.push_back(id);layer.push_back(l);x.push_back(pos[0]);y.push_back(pos[1]);
          truthMultiplicity.push_back(s.truthN);rawMultiplicity.push_back(s.rawN);taggedMultiplicity.push_back(s.taggedN);calibMultiplicity.push_back(s.calibN);
          truthEnergy.push_back(s.truth);hgADC.push_back(s.hg);lgADC.push_back(s.lg);calibEnergy.push_back(s.calib);
        }
        run=calib.run;event=calib.event;eventTime=calib.time;axisValid=savedValid;events.Fill();
      }
      ++fileIndex;
    }
    require(cursor==beam->GetEntries() && cursor>0,"Manifest does not cover beam cohort");
    out->cd();events.Write();
    TNamed("matching_rule","Physical CellID union; event identity and cached signals checked. Truth/calib signed energies summed; raw ADC is maximum per physical strip, not pedestal-subtracted. ADC plots exclude duplicates.").Write();
    TNamed("selection_rule","Cached beam axis-valid cohort; fixed event denominator in all stages. Raw_Hit is already threshold-selected; pre-threshold ADC and rejected fluctuations are unavailable.").Write();
    std::cout<<"Matched "<<cursor<<" sampled events\n";
  }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
