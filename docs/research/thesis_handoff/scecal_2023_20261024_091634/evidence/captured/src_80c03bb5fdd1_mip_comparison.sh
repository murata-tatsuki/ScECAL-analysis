#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

# Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/e-/20GeV
Decode_path=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV
outputDir=../results/mip/simulation/cosmicThreScan/100GeV

rm ../results/comparison.txt
touch ../results/comparison.txt

for i in $( ls ${Decode_path} ); do

  echo -n ${Decode_path}/${i} >> ../results/comparison.txt
  echo -n " " >> ../results/comparison.txt
  
done

var=`cat ../results/comparison.txt`
./MIP_comparison ../results/mip/simulation/cosmicThreScan/100GeV/comparison.root /home/murata_t/ScECAL_BeamTest/ECAL/results/mip_old/jiaxuan/all_auto_muon_v2_hist.root ${var}
