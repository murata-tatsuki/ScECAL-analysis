#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

# Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-
# Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/results/${year}/${beam}/calib/${date_tag} ############ for edited calibration
# Calib_path=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update/Result_MC/calib/e-
Calib_path=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/truth/e-/sps

cog_rage=200
fig_path=../result/resolution/figures/generation/${cog_rage}mm

mkdir -p ${fig_path}

rm ../tmp/resolution_generation.txt
touch ../tmp/resolution_generation.txt

countGeV=$(ls ${Calib_path} -U1 | wc -l)
echo -n ${countGeV} >> ../tmp/resolution_generation.txt
echo -n " " >> ../tmp/resolution_generation.txt

for txt in $( ls ${Calib_path} ); do
  echo -n ${txt} >> ../tmp/resolution_generation.txt
  echo -n " " >> ../tmp/resolution_generation.txt
  mkdir -p ${fig_path}/${txt}
done

for txt in $( ls ${Calib_path} ); do
  count=$(ls ${Calib_path}/${txt} -U1 | wc -l)
  echo ${txt} ${count}
  echo -n ${count} >> ../tmp/resolution_generation.txt
  # echo -n 1 >> ../tmp/resolution_generation.txt
  echo -n " " >> ../tmp/resolution_generation.txt
done

for txt in $( ls ${Calib_path}/*GeV/* ); do
  echo -n ${txt} >> ../tmp/resolution_generation.txt
  echo -n " " >> ../tmp/resolution_generation.txt
#  ./between_files_mip_analysis results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done


var=`cat ../tmp/resolution_generation.txt`
#./between_files_energy_resolution ../results/resolution/energy/Eres_all_normalgain.root ${var}
# ./between_files_energy_resolution_generation ../result/resolution/Eres_simulation_v0.root ${var}
./between_files_energy_resolution_generation ../result/resolution/Eres_generation_${cog_rage}.root ${var} ${cog_rage} 1 ${fig_path}
