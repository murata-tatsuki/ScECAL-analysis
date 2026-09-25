#! /bin/bash

paralell_jobs () {
  track_path=$1
  beam=$2
  outPath=$3

  rm ../../tmp/efficiency_${beam}.txt
  touch ../../tmp/efficiency_${beam}.txt

  for dat_name in $(ls ${track_path}/*.root)
  do
    prefix=`basename $dat_name`
    mkdir -p ${outPath}/${beam}
    outfile=${outPath}/${beam}/${prefix}
    prefix="${prefix:0:-5}"
    fig_path=${outPath}/figures/${beam}/${prefix}
    mkdir -p ${fig_path}

    echo -n ${dat_name} >> ../../tmp/efficiency_${beam}.txt
    echo -n " " >> ../../tmp/efficiency_${beam}.txt
  
    sbatch -o ../jobs/efficiency/test-%A.out --error="../jobs/efficiency/test-%A.err" execute.sl ./channel_efficiecy ${outfile} ${dat_name} ${fig_path}
  done

  var=`cat ../../tmp/efficiency_${beam}.txt`
  outfile=${outPath}/${beam}/all.root
  fig_path=${outPath}/figures/${beam}/all
  mkdir -p ${fig_path}
  sbatch -o ../jobs/efficiency/test-%A.out --error="../jobs/efficiency/test-%A.err" execute.sl ./channel_efficiecy ${outfile} ${var} ${fig_path}
}


rm -rf ../jobs/efficiency
mkdir -p ../jobs/efficiency


# outPath=../../result/trackFit/efficiency

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit/mu-/100GeV
# paralell_jobs ${track_path} sps ${outPath}

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit/mu-/10GeV/
# paralell_jobs ${track_path} ps ${outPath}


# outPath=../../result/trackFit3D/efficiency

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit3D/mu-/100GeV
# paralell_jobs ${track_path} sps ${outPath}

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit3D/mu-/10GeV/
# paralell_jobs ${track_path} ps ${outPath}





## simulation


outPath_=../../result/trackFit3D/efficiency

track_path_=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/trackFit3D/mu-/mu_track/threshold/100GeV
paralell_jobs ${track_path_} simulation ${outPath_}

track_path_=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/trackFit3D/mu-/mu_track/threshold/10GeV
paralell_jobs ${track_path_} simulation ${outPath_}



