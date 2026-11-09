#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_Temperature/run/global_config

rm ../results/hl_allfilename.txt
touch ../results/hl_allfilename.txt

#for i in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed/2023/*/decode/*/*GeV/* ); do
#  echo -n ${i} >> ../results/hl_allfilename.txt
#  echo -n " " >> ../results/hl_allfilename.txt
#done

for i in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed/2023/*/decode/mu-/*GeV/* ); do
  echo -n ${i} >> ../results/hl_allfilename.txt
  echo -n " " >> ../results/hl_allfilename.txt
#  echo ${i}
#  let ii++
done

var=`cat ../results/hl_allfilename.txt`
#./between_files_hl_ratio ../results/hl_ratio/hl_all.root ${var}
./between_files_hl_ratio ../results/hl_ratio/hl_mu_all.root ${var}



#for i in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed/2023/*/decode/pi-/*GeV/* ); do
#  echo -n ${i} >> ../results/hl_allfilename.txt
#  echo -n " " >> ../results/hl_allfilename.txt
##  echo ${i}
##  let ii++
#done
#
#var=`cat ../results/hl_allfilename.txt`
##./between_files_hl_ratio ../results/hl_ratio/hl_all.root ${var}
#./between_files_hl_ratio ../results/hl_ratio/hl_pi_all.root ${var}
