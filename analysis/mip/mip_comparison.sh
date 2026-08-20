#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

# Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/e-/20GeV
Decode_path=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/100GeV/sps
outputDir=../results/mip

rm ../tmp/mip.txt
touch ../tmp/mip.txt

for i in $( ls ${Decode_path} ); do

  echo -n ${Decode_path}/${i} >> ../tmp/mip.txt
  echo -n " " >> ../tmp/mip.txt
  
done

var=`cat ../tmp/mip.txt`
./MIP_calibration_simulation ../result/mip/MIP_simultaion_new_.root ${var}
