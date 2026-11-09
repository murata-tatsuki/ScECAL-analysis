#! /bin/bash

data_file=20230428_0829_Calib
#data_file=20230503_0737_Calib
#data_file=20230503_0835_Calib
#data_file=20230504_1425_Calib

#data_file=20230516_1400_Calib
#data_file=20230517_1550_Calib
#data_file=20230517_1606_Calib
#data_file=20230518_1533_Calib
#data_file=20230519_1140_Calib
#data_file=20230520_1556_Calib
#data_file=20230521_1455_Calib
#data_file=20230522_1249_Calib
#data_file=20230523_1435_Calib
#data_file=20230524_1118_Calib

data_files=(20230428_0829_Calib 20230503_0737_Calib 20230503_0835_Calib 20230504_1425_Calib 20230516_1400_Calib 20230517_1550_Calib 20230517_1606_Calib 20230518_1533_Calib 20230519_1140_Calib 20230520_1556_Calib 20230521_1455_Calib 20230522_1249_Calib 20230523_1435_Calib 20230524_1118_Calib)
dacList=(2100 2150 2200 2250 2300 2350 2400 2450 2500)
# dacList=(2500)

beam=sps
#beam=ps

rm ../results/led/LED_allfilename.txt
touch ../results/led/LED_allfilename.txt


#for txt in $( ls /home/jap/data2023/SPS/ECAL/decode/LEDCali/*dac* ); do
#for txt in $( ls /mnt2/jap/data2023/SPS/ECAL/decode/LEDCali/${data_file}/*dac* ); do
#for txt in $( ls /mnt2/jap/data2023/${beam}/ECAL/decode/LEDCali/${data_file}/*dac* ); do
i=0
for file in ${data_files[@]}; do
  # rm ../results/led/LED_allfilename.txt
  # touch ../results/led/LED_allfilename.txt

  #echo ${file}
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
       echo ./between_files_LED ../results/led/${file}/${file}_LED_${dac}.root ${outPitDir} ${var}
      # echo ./between_files_LED ../results/led/${file}_LED_${dac}.root ${outPitDir} ${var}
    # fi

  done
  let i++
done
#var=`cat ../results/led/LED_allfilename.txt`
#./between_files_LED ../results/led/${data_file}.root ${var}

#for txt in $( ls /megraid01/users/murata_t/scecal/ScECAL_BeamTest/ECAL/results/led/2023*_LED_all.root ); do
##  echo ${txt}
#  echo -n ${txt} >> ../results/led/LED_allfilename.txt
#  echo -n " " >> ../results/led/LED_allfilename.txt
#done

#var=`cat ../results/led/LED_allfilename.txt`
#./between_files_LED ../results/led/LED_all.root ${var}

