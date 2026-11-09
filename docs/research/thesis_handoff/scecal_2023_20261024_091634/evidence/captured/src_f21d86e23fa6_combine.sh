#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/mu-/100GeV
Calib_path=${data_dir}/${year}/sps/calib_simulation/e-/20GeV

rm ../results/combining.txt
touch ../results/combining.txt
ii=0
#for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do
for i in $( ls ${Calib_path} ); do
  echo -n ${Calib_path}/${i} >> ../results/combining.txt
  echo -n " " >> ../results/combining.txt
done

var=`cat ../results/combining.txt`
./combine_all_data ../results/mip/combine_e.root ${var}

