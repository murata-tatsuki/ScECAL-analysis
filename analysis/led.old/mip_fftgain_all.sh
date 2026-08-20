#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config
year=2023

#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
Decode_path=${data_dir}/${year}/sps/decode/mu-/100GeV
SSA_path=${data_dir}/${year}/sps/ssa/mu-/100GeV

# << COMMENTOUT
rm ../MIP_gain_allfilename.txt
touch ../MIP_gain_allfilename.txt
ii=0
#for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do
for i in $( ls ${Decode_path} ); do
  echo -n ${SSA_path}/${i} >> ../MIP_gain_allfilename.txt
  echo -n " " >> ../MIP_gain_allfilename.txt
  echo -n ${Decode_path}/${i} >> ../MIP_gain_allfilename.txt
  echo -n " " >> ../MIP_gain_allfilename.txt
done
# COMMENTOUT

var=`cat ../MIP_gain_allfilename.txt`
./selection_mip ../result/mip.root ${var}

# #Decode_path=${data_dir}/${year}/ps/simpleCalib/mu-/10GeV
# Decode_path=${data_dir}/${year}/ps/decode/mu-/10GeV
# SSA_path=${data_dir}/${year}/ps/ssa/mu-/10GeV

# for i in $( ls ${Decode_path} ); do
#   echo -n ${SSA_path}/${i} >> ../MIP_gain_allfilename.txt
#   echo -n " " >> ../MIP_gain_allfilename.txt
#   echo -n ${Decode_path}/${i} >> ../MIP_gain_allfilename.txt
#   echo -n " " >> ../MIP_gain_allfilename.txt

# #  echo ${ii} " " ${i}
# #  let ii++

# #  ./between_files_mip_analysis ../mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done


#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Decode_path=${data_dir}/${year}/sps/decode/auto_gain/mu-/100GeV
# SSA_path=${data_dir}/${year}/sps/ssa/auto_gain/mu-/100GeV

# for i in $( ls ${Decode_path} ); do
#   echo -n ${SSA_path}/${i} >> ../MIP_gain_allfilename.txt
#   echo -n " " >> ../MIP_gain_allfilename.txt
#   echo -n ${Decode_path}/${i} >> ../MIP_gain_allfilename.txt
#   echo -n " " >> ../MIP_gain_allfilename.txt

# #  echo ${ii} " " ${i}
# #  let ii++

# #  ./between_files_mip_analysis ../mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done

# var=`cat ../MIP_gain_allfilename.txt`
#var=`cat ../MIP_normalgain_allfilename_center.txt`
#./between_files_MIP_normal ../mip/MIP_normalgain_all_nocut_2chn.root ${var}
# ./between_files_MIP_newConvoluted ../mip/withPedestal/MIP_gain_all_binWidth5.root ${var}
# ./between_files_MIP_newConvoluted ../mip/withPedestal/test.root ${var}

