#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config


#Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-
Calib_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/calib/e-
Calib_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/calib/e-
#Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/tmp/${year}/${beam}/calib/${date_tag} ############ for edited calibration

Calib_path_sim=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/sps

rm ../tmp/discrepancy.txt
touch ../tmp/discrepancy.txt

countGeV_data=$(ls ${Calib_path_sps}/40GeV -U1 | wc -l)
echo -n ${countGeV_data} >> ../tmp/discrepancy.txt
echo -n " " >> ../tmp/discrepancy.txt

countGeV_sim=$(ls ${Calib_path_sim}/40GeV -U1 | wc -l)
echo -n ${countGeV_sim} >> ../tmp/discrepancy.txt
echo -n " " >> ../tmp/discrepancy.txt

for txt in $( ls ${Calib_path_sps}/40GeV ); do
  echo -n ${Calib_path_sps}/40GeV/${txt} >> ../tmp/discrepancy.txt
  echo -n " " >> ../tmp/discrepancy.txt
done

for txt in $( ls ${Calib_path_sim}/40GeV ); do
  echo -n ${Calib_path_sim}/40GeV/${txt} >> ../tmp/discrepancy.txt
  echo -n " " >> ../tmp/discrepancy.txt
done

var=`cat ../tmp/discrepancy.txt`
./discrepancy ../result/discrepancy/edep.root ${var}




# Calib_path_2022=${data_dir}/2022/sps/calib/e+

# rm ../tmp/discrepancy.txt
# touch ../tmp/discrepancy.txt

# countGeV=$(ls ${Calib_path_2022} -U1 | wc -l)
# echo -n ${countGeV} >> ../tmp/discrepancy.txt
# echo -n " " >> ../tmp/discrepancy.txt

# for txt in $( ls ${Calib_path_2022} ); do
#   echo -n ${txt} >> ../tmp/discrepancy.txt
#   echo -n " " >> ../tmp/discrepancy.txt
# done

# for txt in $( ls ${Calib_path_2022} ); do
#   count=$(ls ${Calib_path_2022}/${txt} -U1 | wc -l)
#   echo -n ${count} >> ../tmp/discrepancy.txt
#   echo -n " " >> ../tmp/discrepancy.txt
# done

# for txt in $( ls ${Calib_path_2022}/*GeV/ECAL*.root ); do
#   echo -n ${txt} >> ../tmp/discrepancy.txt
#   echo -n " " >> ../tmp/discrepancy.txt
# #  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done

# var=`cat ../tmp/discrepancy.txt`
# ./between_files_energy_resolution_2022 ../result/resolution/energy/Eres_normal_2022_all.root ${var}


