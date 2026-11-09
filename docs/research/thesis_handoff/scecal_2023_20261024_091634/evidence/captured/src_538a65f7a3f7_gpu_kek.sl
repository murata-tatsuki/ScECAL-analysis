#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch -a 0-20%10 example/analyzer_array.sl
#             sbatch -o test-%A_%a.out --error="test-%A_%a.err" gpu_kek.sl
#
###############################################################################

#SBATCH --partition=all
#SBATCH --ntasks=1       
#SBATCH --job-name=analyzer
#SBATCH --array=0
#SBATCH --output="analyzer-%A_%a.out"
#SBATCH --error="analyzer-%A_%a.err"

MY_HOST=`hostname`
MY_DATE=`date`
#
# Disable core dump. (Remove this line if you want to create core files)
#
ulimit -c 0

echo "Running on $MY_HOST at $MY_DATE"
echo "================================================================"

#
# Print envirionment
#
#echo "Running environment":
#env | sort
#echo "================================================================"

#
# Replace the start run number and configuration file by yours
#
#! /bin/bash

## 一番いいfitを選び出すもの

# source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/run/global_config

# #Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Decode_path=${data_dir}/${year}/sps/decode/mu-/100GeV
# SSA_path=${data_dir}/${year}/sps/ssa/mu-/100GeV

rm ../results/MIP_gain_allfilename.txt
touch ../results/MIP_gain_allfilename.txt
# ii=0
# #for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do

# Decode_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/decode/mipScan/layer_except_10_11_14_15_20_21_26_27
Decode_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/trackFit_layer16/mipScan/layer_except_10_11_14_15_20_21_26_27
SSA_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27

for i in $( ls ${Decode_path} ); do
  # rm ../results/MIP_gain_allfilename.txt
  # touch ../results/MIP_gain_allfilename.txt

  echo -n ${SSA_path}/${i} >> ../results/MIP_gain_allfilename.txt
  echo -n " " >> ../results/MIP_gain_allfilename.txt
  echo -n ${Decode_path}/${i} >> ../results/MIP_gain_allfilename.txt
  echo -n " " >> ../results/MIP_gain_allfilename.txt

  var=`cat ../results/MIP_gain_allfilename.txt`
  # echo ./selection ../results/mip/kek/fit/${i} ${var}
  # echo ./selection_1hit ../results/mip/kek/fit_1hit/${i} ${var}

# #  echo ${ii} " " ${i}
# #  let ii++

#  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done

var=`cat ../results/MIP_gain_allfilename.txt`
# ./selection ../results/mip/kek/fit/allRun.root ${var}
# ./selection_1hit ../results/mip/kek/fit_1hit/allRun.root ${var}
# ./selection_trackFit_ ../results/mip/kek/fit_trackFit_layer16/allRun.root ${var}
./selection_trackFit_ ../results/mip/kek/fit_trackFit_layer16/allRun_offtrack.root ${var}

# ./selection ../results/mip/kek/fit/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/decode/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root
# ./selection_1hit ../results/mip/kek/fit_1hit/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/decode/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root
# ./selection_trackFit_ ../results/mip/kek/fit_trackFit_layer16_/ECAL_Run171_20250605_162329_offtrack.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/trackFit_layer16_/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run171_20250605_162329.root

# ./selection ../results/mip/kek/fit/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/decode/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root
# ./selection_1hit ../results/mip/kek/fit_1hit/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/decode/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root
# ./selection_trackFit_ ../results/mip/kek/fit_trackFit/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ssa/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root /megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/trackFit_/mipScan/layer_except_10_11_14_15_20_21_26_27/ECAL_Run172_20250605_183026.root

# var=`cat ../results/MIP_gain_allfilename.txt`
# echo ./selection ../results/mip/kek/select_tempCorrect_25C.root ${var}




RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
