#! /bin/bash

paralell_jobs () {
  track_path=$1
  beam=$2

  rm ../../tmp/efficiency_${beam}.txt
  touch ../../tmp/efficiency_${beam}.txt

  for dat_name in $(ls ${track_path})
  do
    prefix=`basename $dat_name`
    mkdir -p ../../result/trackFit/efficiency/${beam}
    outfile=../../result/trackFit/efficiency/${beam}/${prefix}
    prefix="${prefix:0:-5}"
    fig_path=../../result/trackFit/efficiency/figures/${beam}/${prefix}
    mkdir -p ${fig_path}

    echo -n ${track_path}/${dat_name} >> ../../tmp/efficiency_${beam}.txt
    echo -n " " >> ../../tmp/efficiency_${beam}.txt
  
    sbatch -o jobs/efficiency/test-%A.out --error="jobs/efficiency/test-%A.err" execute.sl ./channel_efficiecy ${outfile} ${track_path}/${dat_name} ${fig_path}
  done

  var=`cat ../../tmp/efficiency_${beam}.txt`
  outfile=../../result/trackFit/efficiency/${beam}/all.root
  fig_path=../../result/trackFit/efficiency/figures/${beam}/all
  mkdir -p ${fig_path}
  sbatch -o jobs/efficiency/test-%A.out --error="jobs/efficiency/test-%A.err" execute.sl ./channel_efficiecy ${outfile} ${var} ${fig_path}
}


rm -rf jobs/efficiency
mkdir -p jobs/efficiency


track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit/mu-/100GeV
paralell_jobs ${track_path} sps

track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit/mu-/10GeV/
paralell_jobs ${track_path} ps



