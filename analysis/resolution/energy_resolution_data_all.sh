#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config


#Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-_noIntercept
Calib_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/calib/e-
Calib_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/calib/e-
#Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/tmp/${year}/${beam}/calib/${date_tag} ############ for edited calibration

cog_rage=20
fig_path=../result/resolution/figures/data/${cog_rage}mm

rm ../tmp/Eres_allfilename_normalgain.txt
touch ../tmp/Eres_allfilename_normalgain.txt

countGeV_sps=$(ls ${Calib_path_sps} -U1 | wc -l)
countGeV_ps=$(ls ${Calib_path_ps} -U1 | wc -l)
countGeV=$((${countGeV_sps}+${countGeV_ps}))
echo -n ${countGeV} >> ../tmp/Eres_allfilename_normalgain.txt
echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt

for txt in $( ls ${Calib_path_sps} ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
  mkdir -p ${fig_path}/${txt}
done
for txt in $( ls ${Calib_path_ps} ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
  mkdir -p ${fig_path}/${txt}
done

for txt in $( ls ${Calib_path_sps} ); do
  count=$(ls ${Calib_path_sps}/${txt} -U1 | wc -l)
  echo -n ${count} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
done
for txt in $( ls ${Calib_path_ps} ); do
  count=$(ls ${Calib_path_ps}/${txt} -U1 | wc -l)
  echo -n ${count} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
done

for txt in $( ls ${Calib_path_sps}/*GeV/ECAL*.root ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
#  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done
for txt in $( ls ${Calib_path_ps}/*GeV/ECAL*.root ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
done


var=`cat ../tmp/Eres_allfilename_normalgain.txt`
./between_files_energy_resolution ../result/resolution/Eres_data_v3_${cog_rage}.root ${var} ${cog_rage} 1 ${fig_path}




# Calib_path_2022=${data_dir}/2022/sps/calib/e+

# rm ../tmp/Eres_allfilename_normalgain.txt
# touch ../tmp/Eres_allfilename_normalgain.txt

# countGeV=$(ls ${Calib_path_2022} -U1 | wc -l)
# echo -n ${countGeV} >> ../tmp/Eres_allfilename_normalgain.txt
# echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt

# for txt in $( ls ${Calib_path_2022} ); do
#   echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
#   echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
# done

# for txt in $( ls ${Calib_path_2022} ); do
#   count=$(ls ${Calib_path_2022}/${txt} -U1 | wc -l)
#   echo -n ${count} >> ../tmp/Eres_allfilename_normalgain.txt
#   echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
# done

# for txt in $( ls ${Calib_path_2022}/*GeV/ECAL*.root ); do
#   echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
#   echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
# #  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done

# var=`cat ../tmp/Eres_allfilename_normalgain.txt`
# ./between_files_energy_resolution_2022 ../result/resolution/energy/Eres_normal_2022_all.root ${var}


