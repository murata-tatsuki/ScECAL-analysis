#!/bin/bash
###############################################################################
# A sample SLURM script to run analyzer on 1 node
#
#  usage:
#        (0) Edit this script
#        (1) Submit the job.
#            % cd $MEG2SYS/analyzer
#            % sbatch -a 0-20%10 example/analyzer_array.sl
#             sbatch -o test-%A_%a.out --error="test-%A_%a.err" gpu.sl
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

source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/run/global_config

#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
Decode_path=${data_dir}/${year}/sps/decode/mu-/100GeV
SSA_path=${data_dir}/${year}/sps/ssa/mu-/100GeV

rm ../results/MIP_gain_allfilename.txt
touch ../results/MIP_gain_allfilename.txt
ii=0
#for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do
for i in $( ls ${Decode_path} ); do
  echo -n ${SSA_path}/${i} >> ../results/MIP_gain_allfilename.txt
  echo -n " " >> ../results/MIP_gain_allfilename.txt
  echo -n ${Decode_path}/${i} >> ../results/MIP_gain_allfilename.txt
  echo -n " " >> ../results/MIP_gain_allfilename.txt

#  echo ${ii} " " ${i}
#  let ii++

#  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done


# #Decode_path=${data_dir}/${year}/ps/simpleCalib/mu-/10GeV
# Decode_path=${data_dir}/${year}/ps/decode/mu-/10GeV
# SSA_path=${data_dir}/${year}/ps/ssa/mu-/10GeV

# for i in $( ls ${Decode_path} ); do
#   echo -n ${SSA_path}/${i} >> ../results/MIP_gain_allfilename.txt
#   echo -n " " >> ../results/MIP_gain_allfilename.txt
#   echo -n ${Decode_path}/${i} >> ../results/MIP_gain_allfilename.txt
#   echo -n " " >> ../results/MIP_gain_allfilename.txt

# #  echo ${ii} " " ${i}
# #  let ii++

# #  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done


#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
# Decode_path=${data_dir}/${year}/sps/decode/auto_gain/mu-/100GeV
# SSA_path=${data_dir}/${year}/sps/ssa/auto_gain/mu-/100GeV

# for i in $( ls ${Decode_path} ); do
#   echo -n ${SSA_path}/${i} >> ../results/MIP_gain_allfilename.txt
#   echo -n " " >> ../results/MIP_gain_allfilename.txt
#   echo -n ${Decode_path}/${i} >> ../results/MIP_gain_allfilename.txt
#   echo -n " " >> ../results/MIP_gain_allfilename.txt

# #  echo ${ii} " " ${i}
# #  let ii++

# #  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
# done

var=`cat ../results/MIP_gain_allfilename.txt`
#var=`cat ../results/MIP_normalgain_allfilename_center.txt`
#./between_files_MIP_normal ../results/mip/MIP_normalgain_all_nocut_2chn.root ${var}
# ./between_files_MIP_newConvoluted ../results/mip/withPedestal/MIP_gain_all_binWidth5.root ${var}
# ./between_files_MIP_newConvoluted ../results/mip/withPedestal/test.root ${var}
./selection ../results/mip/withPedestal/select_tempCorrect_25C.root ${var}





RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
