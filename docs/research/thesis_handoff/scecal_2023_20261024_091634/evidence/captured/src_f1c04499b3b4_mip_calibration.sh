#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

# Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/mu-/100GeV
# Calib_path=${data_dir}/${year}/sps/calib_simulation/e-/20GeV
Decode_path=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV
outputDir=../results/mip/simulation/cosmicThreScan/100GeV/fitting/all_position

mkdir -p ${outputDir}

for i in $( ls ${Decode_path} ); do
  if [[ "${i}" == "mu-_100GeV_100000_thre_0.root" ]]; then
    continue
  fi
  # ./MIP_calibration_simulation ${outputDir}/${i} ${Decode_path}/${i}
done

var=$( ls ${Decode_path}/*thre_0.root )
echo ${var}


./MIP_calibration_simulation ${outputDir}/mu-_100GeV_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-100_-100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-100_-50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-100_0_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-100_100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-100_50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-50_-100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-50_-50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-50_0_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-50_100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos-50_50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos0_-100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos0_-50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos0_0_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos0_100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos0_50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos100_-100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos100_-50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos100_0_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos100_100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos100_50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos50_-100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos50_-50_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos50_0_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos50_100_thre0.root /megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/mu-/cosmicThreScan/100GeV/mu-_100GeV_pos50_50_thre0.root

# ./MIP_calibration_simulation ${outputDir}/${i} ${Decode_path}/${i}

