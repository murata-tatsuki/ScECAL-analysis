#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config


#Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-_noIntercept
Calib_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/decode/e-
Calib_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/decode/e-
#Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/tmp/${year}/${beam}/calib/${date_tag} ############ for edited calibration

fig_path=../result/threshold/figures/allGeVs
mkdir -p ${fig_path}/gaus/chips
mkdir -p ${fig_path}/exp/chips
mkdir -p ${fig_path}/landau/chips
mkdir -p ${fig_path}/gaus/layers
mkdir -p ${fig_path}/exp/layers
mkdir -p ${fig_path}/landau/layers

rm ../tmp/threshold.txt
touch ../tmp/threshold.txt

countGeV_sps=$(ls ${Calib_path_sps} -U1 | wc -l)
countGeV_ps=$(ls ${Calib_path_ps} -U1 | wc -l)
countGeV=$((${countGeV_sps}+${countGeV_ps}))
echo -n ${countGeV} >> ../tmp/threshold.txt
echo -n " " >> ../tmp/threshold.txt

for txt in $( ls ${Calib_path_sps} ); do
  echo -n ${txt} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
  # mkdir -p ${fig_path}/${txt}
done
for txt in $( ls ${Calib_path_ps} ); do
  echo -n ${txt} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
  # mkdir -p ${fig_path}/${txt}
done

for txt in $( ls ${Calib_path_sps} ); do
  count=$(ls ${Calib_path_sps}/${txt} -U1 | wc -l)
  echo -n ${count} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
done
for txt in $( ls ${Calib_path_ps} ); do
  count=$(ls ${Calib_path_ps}/${txt} -U1 | wc -l)
  echo -n ${count} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
done

for txt in $( ls ${Calib_path_sps}/*GeV/ECAL*.root ); do
  echo -n ${txt} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
#  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done
for txt in $( ls ${Calib_path_ps}/*GeV/ECAL*.root ); do
  echo -n ${txt} >> ../tmp/threshold.txt
  echo -n " " >> ../tmp/threshold.txt
done


var=`cat ../tmp/threshold.txt`
./threshold_finding ../result/threshold/threshold.root ${var} ${fig_path}


