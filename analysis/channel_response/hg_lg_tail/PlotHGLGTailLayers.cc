#include "HGLGTailLayerMaps.hh"
#include <iostream>
int main(int argc,char** argv)try{
 require(argc==2,"Usage: PlotHGLGTailLayers RESULT_DIRECTORY");
 fs::path out=argv[1];auto input=openFile((out/"tail_study.root").string());
 auto* data=object<TH2D>(*input,"data_low_tail_rate");auto* mc=object<TH2D>(*input,"mc_low_tail_rate");
 auto* beam=object<TParameter<double>>(*input,"beam_GeV");
 std::ifstream table(out/"channels.tsv");require(bool(table),"Cannot open channels.tsv");
 std::string line;std::getline(table,line);require(line.rfind("cell\tgood\tcommon\t",0)==0,"Unexpected channel table header");
 std::set<int> common,seen;
 while(std::getline(table,line)){
  std::istringstream row(line);int cell,good,shared;require(bool(row>>cell>>good>>shared),"Malformed channel row");
  require(connected(cell)&&physical(cell)==cell&&seen.insert(cell).second,"Invalid/duplicate channel");
  if(shared){require(good==1,"Common channel must be good");common.insert(cell);}
 }
 require(seen.size()==6300,"Incomplete channel table");
 require(common.size()==size_t(object<TParameter<int>>(*input,"common_channels")->GetVal()),"Common channel count mismatch");
 gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);
 drawTailLayerMaps(*data,*mc,common,out,beam->GetVal());
 std::cout<<"Saved 4 layer-map PNGs and hg_lg_tail_layer_maps.root in "<<out<<std::endl;
 return 0;
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<std::endl;return 1;}
