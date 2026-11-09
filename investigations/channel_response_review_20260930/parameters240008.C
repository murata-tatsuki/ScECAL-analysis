#include <TFile.h>
#include <TTree.h>
#include <set>
#include <iostream>
void parameters240008() {
 TFile f("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/all_hl_electron2023_hlratio_v3.root");
 auto* t=f.Get<TTree>("InterCalib");int cell=0;double slope=0,intercept=0,xmax=0;
 t->SetBranchAddress("CellID",&cell);t->SetBranchAddress("Slope",&slope);t->SetBranchAddress("Intercept",&intercept);t->SetBranchAddress("XMax",&xmax);
 std::set<int> ids={240008};
 std::cout<<"cell slope intercept xmax raw_high_switch"<<std::endl;
 for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);if(ids.count(cell))std::cout<<cell<<' '<<slope<<' '<<intercept<<' '<<xmax<<' '<<xmax-600<<std::endl;}
}
