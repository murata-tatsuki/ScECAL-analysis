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

data_files=(20230428_0829_Calib 20230503_0737_Calib 20230503_0835_Calib 20230504_1425_Calib 20230516_1400_Calib 20230517_1550_Calib 20230517_1606_Calib 20230518_1533_Calib 20230519_1140_Calib 20230520_1556_Calib 20230521_1455_Calib 20230522_1249_Calib 20230523_1435_Calib 20230524_1118_Calib)
dacList=(2100 2150 2200 2250 2300 2350 2400 2450 2500)
# dacList=(2500)

beam=sps
#beam=ps

rm ../results/led/LED_allfilename.txt
touch ../results/led/LED_allfilename.txt
i=0
for file in ${data_files[@]}; do
  if [[ i -ge 4 ]]
  then
    beam=ps
  fi
  #echo ${beam}
  for dac in ${dacList[@]}; do 
    rm ../results/led/LED_allfilename.txt
    touch ../results/led/LED_allfilename.txt
    for txt in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed_oldFormat/2023/${beam}/decode/LEDCali/${file}/*dac${dac}* ); do
    #  echo ${txt}
      echo -n ${txt} >> ../results/led/LED_allfilename.txt
      echo -n " " >> ../results/led/LED_allfilename.txt
    done
    # file_name=`basename $file`
    # prefix=${file_name%.dat}
    mkdir -p ../results/led/${file}
    outPitDir=../results/led/pic/${file}/${dac}
    mkdir -p ${outPitDir}

    var=`cat ../results/led/LED_allfilename.txt`
    # if [[ ${file} == *20230504* ]]
    # then
      ./between_files_LED ../results/led/${file}/${file}_LED_${dac}.root ${outPitDir} ${var}
      # echo ./between_files_LED ../results/led/${file}_LED_${dac}.root ${outPitDir} ${var}
    # fi

  done
  let i++
done




RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
