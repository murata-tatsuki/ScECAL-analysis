#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-
Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/results/${year}/${beam}/calib/${date_tag} ############ for edited calibration

Calib_path=/megraid01/users//ECAL_data/analysed/2023/ps

rm ../tmp/Eres_allfilename_normalgain.txt
touch ../tmp/Eres_allfilename_normalgain.txt

countGeV=$(ls ${Calib_path} -U1 | wc -l)
echo -n ${countGeV} >> ../tmp/Eres_allfilename_normalgain.txt
echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt

for txt in $( ls ${Calib_path} ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
done

for txt in $( ls ${Calib_path} ); do
  count=$(ls ${Calib_path}/${txt} -U1 | wc -l)
  echo ${txt} ${count}
  echo -n ${count} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
done

for txt in $( ls ${Calib_path}/*GeV/*.root ); do
  echo -n ${txt} >> ../tmp/Eres_allfilename_normalgain.txt
  echo -n " " >> ../tmp/Eres_allfilename_normalgain.txt
#  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done

var=`cat ../tmp/Eres_allfilename_normalgain.txt`
#./between_files_energy_resolution ../tmp/resolution/energy/Eres_all_normalgain.root ${var}
./between_files_energy_resolution ../result/result/resolution/Eres_all.root ${var}
