#！/bin/sh
path="$PWD"

inputfile="/megraid01/users/murata_t/scecal/ScECAL_CR/results/eventAlign/"
outputfile="/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/Double_readout/results/position_resolution/"

rm results/LED_allfilename.txt
touch results/LED_allfilename.txt



###### for file loop
#i=0
#for file in ${inputfile}*.root
#do
#  DIR_OUT1=${outputfile}$(basename ${file})
#  DIR_IN=${inputfile}$(basename ${file})
#  
#  echo ${file}
# 
#  if [ ${file} = "/megraid01/users/murata_t/scecal/ScECAL_CR/results/eventAlign/20210313_1600_cosmicRayTst.root" -o ${i} -ge 1 ] ; then
#    let i++
#    echo ${file}
#    continue
#  else
#    echo ${file}
#  fi
#  ./position_resolution ${DIR_OUT1} ${DIR_IN}
#  echo -n ${DIR_IN} >> results/position_resolution/allfilename.txt
#  echo -n " " >> results/position_resolution/allfilename.txt
#  
#  echo -n ${DIR_OUT1} >> results/position_resolution/alloutputname.txt
#  echo -n " " >> results/position_resolution/alloutputname.txt
#
#done

for txt in $( ls /megraid01/users/data_beamtest_SPS2022/ECAL_data/decode/LEDCalib/20221031_2214_Calib/*_dac* ); do   
#  echo ${txt}
  echo -n ${txt} >> results/LED_allfilename.txt
  echo -n " " >> results/LED_allfilename.txt
done

var=`cat results/LED_allfilename.txt`
./between_files_LED ../results/LED_all.root ${var}


