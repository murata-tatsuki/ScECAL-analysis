#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
Decode_path=${data_dir}/${year}/sps/decode/mu-/100GeV
SSA_path=${data_dir}/${year}/sps/ssa/mu-/100GeV

rm ../tmp/mip.txt
touch ../tmp/mip.txt
ii=0
#for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do
for i in $( ls ${Decode_path} ); do
  echo -n ${SSA_path}/${i} >> ../tmp/mip.txt
  echo -n " " >> ../tmp/mip.txt
  echo -n ${Decode_path}/${i} >> ../tmp/mip.txt
  echo -n " " >> ../tmp/mip.txt

#  echo ${ii} " " ${i}
#  let ii++

#  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done

var=`cat ../tmp/mip.txt`
./MIP_calibration_data ../result/mip/mip_data.root ${var}
