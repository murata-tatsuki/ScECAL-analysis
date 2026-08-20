#! /bin/bash

## 一番いいfitを選び出すもの

data_files=(20230428_0829_Calib 20230503_0737_Calib 20230503_0835_Calib 20230504_1425_Calib 20230516_1400_Calib 20230517_1550_Calib 20230517_1606_Calib 20230518_1533_Calib 20230519_1140_Calib 20230520_1556_Calib 20230521_1455_Calib 20230522_1249_Calib 20230523_1435_Calib 20230524_1118_Calib)
dacList_=(2100 2150 2200 2250 2300 2350 2400 2450 2500)
dacList_20230503_0835_Calib=(2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100)
# dacList=(2500)

paralell_jobs () {
  beam=$1
  particle=$2

  energy=100GeV
  if [ $beam = 'ps' ];then
      energy=10GeV
  fi 
  rm ../tmp/beam_allfilename.txt
  touch ../tmp/beam_allfilename.txt

  # for txt in $( ls /megraid01/users/data_beamtest/ECAL_data/analysed/2023/${beam}/decode/${particle}/${energy}/* ); do
  #   # echo ${txt}
  #   echo -n ${txt} >> ../tmp/beam_allfilename.txt
  #   echo -n " " >> ../tmp/beam_allfilename.txt
  # done

  path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/${beam}/decode/${particle}

  for ene in $( ls ${path} ); do
    for txt in $( ls ${path}/${ene} );do

      # echo ${txt}
      echo -n ${path}/${ene}/${txt} >> ../tmp/beam_allfilename.txt
      echo -n " " >> ../tmp/beam_allfilename.txt

    done
  done

  mkdir -p ../result/gain/beam/${beam}/${particle}
  outPitDir=../result/gain/figures/beam/${beam}/all/${particle}
  mkdir -p ${outPitDir}/histo
  mkdir -p ${outPitDir}/temp/graph
  mkdir -p ${outPitDir}/temp/histo

  var=`cat ../tmp/beam_allfilename.txt`
  sbatch -o jobs/test-%A.out --error="jobs/test-%A.err" gpu.sl ./combine_files_beam ../result/gain/beam/${beam}_${particle}_all.root ${outPitDir} ${var}
}


mkdir -p jobs


paralell_jobs sps mu-
paralell_jobs ps mu-
# paralell_jobs sps e-
paralell_jobs ps e-


