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


rm file_list/LED_allfilename.txt
touch file_list/LED_allfilename.txt


beam=sps
i=0
filenum=0
for file in ${data_files[@]}; do
  rm file_list/LED_allfilename.txt
  touch file_list/LED_allfilename.txt

  #echo ${file}
  # echo ${i}
  if [[ i -ge 4 ]]
  then
    beam=ps
  fi
  let i++

  LED_path=../results/led/${file}
  countLED=$(ls ${LED_path} -U1 | wc -l)
  if [ ${countLED} -ne 9 ]
  then
    echo ${LED_path} does not have enough number of files 
    continue
  fi

  let filenum++

  for dac in ${dacList[@]}; do 

    txt=../results/led/${file}/${file}_LED_${dac}.root
    echo ${txt}

    echo -n ${txt} >> file_list/LED_allfilename.txt
    echo -n " " >> file_list/LED_allfilename.txt
    # file_name=`basename $file`
    # prefix=${file_name%.dat}
    # outPitDir=../results/led/pic/${file}/${dac}
    # mkdir -p ${outPitDir}
    # var=`cat file_list/LED_allfilename.txt`
    # if [[ ${file} == *20230504* ]]
    # then
    # echo ./between_files_LED ../results/led/${file}_LED_all.root ${outPitDir} ${var}
    # fi
  done

  if [[ i -le 4 ]]
  then
    continue
  fi

  var=`cat file_list/LED_allfilename.txt`
  echo ./scanning ../results/led/${file}_LED.root ${var}
done

./scanning ../results/led/20230518_1533_Calib_LED_bestfit.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2100.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2150.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2200.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2250.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2300.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2350.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2400.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2450.root ../results/led/20230518_1533_Calib/20230518_1533_Calib_LED_2500.root




RET=$?
echo "Exit status: $RET"
echo "================================================================"
END_DATE=`date`
echo "End of run is $END_DATE ($SECONDS sec)"
