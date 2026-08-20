#! /bin/bash

## 一番いいfitを選び出すもの

data_files=(20230428_0829_Calib 20230503_0737_Calib 20230503_0835_Calib 20230504_1425_Calib 20230516_1400_Calib 20230517_1550_Calib 20230517_1606_Calib 20230518_1533_Calib 20230519_1140_Calib 20230520_1556_Calib 20230521_1455_Calib 20230522_1249_Calib 20230523_1435_Calib 20230524_1118_Calib)
dacList_=(2100 2150 2200 2250 2300 2350 2400 2450 2500)
dacList_20230503_0835_Calib=(2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100)
# dacList=(2500)

# rm -r jobs
mkdir -p jobs

beam=sps
#beam=ps
ledDir=LEDCali

rm ../tmp/LED_allfilename.txt
touch ../tmp/LED_allfilename.txt
i=0
for file in ${data_files[@]}; do
  if [[ i -ge 4 ]]
  then
    beam=ps
    ledDir=lightCalibData
  fi
  dacList=("${dacList_[@]}")
  if [[ i -eq 2 ]]
  then
    dacList=("${dacList_20230503_0835_Calib[@]}")
  fi
  for dac in ${dacList[@]}; do 
    rm ../tmp/LED_allfilename.txt
    touch ../tmp/LED_allfilename.txt
    for txt in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed/2023/${beam}/decode/${ledDir}/${file}/*dac${dac}* ); do
    #  echo ${txt}
      echo -n ${txt} >> ../tmp/LED_allfilename.txt
      echo -n " " >> ../tmp/LED_allfilename.txt
    done
    # file_name=`basename $file`
    # prefix=${file_name%.dat}
    mkdir -p ../result/gain/led/${file}
    outPitDir=../result/gain/figures/led/${file}/${dac}
    mkdir -p ${outPitDir}

    var=`cat ../tmp/LED_allfilename.txt`
    # if [[ ${file} == *20230504* ]]
    # then
      sbatch -o jobs/test-%A.out --error="jobs/test-%A.err" gpu.sl ./combine_files_LED ../result/gain/led/${file}/${file}_LED_${dac}.root ${outPitDir} ${var}
      # echo ./between_files_LED ../result/led/${file}_LED_${dac}.root ${outPitDir} ${var}
    # fi

  done
  let i++
done







# var=`cat file_list/LED_allfilename.txt`
# ./between_files_LED_scan ../result/LED_scan.root ${filenum} ${var}


