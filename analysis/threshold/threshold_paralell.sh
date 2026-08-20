#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

paralell_jobs () {
  Calib_path=$1
  fig_path=$2
  
  for txt in $( ls ${Calib_path} ); do
    rm ../tmp/threshold_${txt}.txt
    touch ../tmp/threshold_${txt}.txt
    mkdir -p ${fig_path}/${txt}
    mkdir -p ${fig_path}/${txt}/gaus/chips
    mkdir -p ${fig_path}/${txt}/exp/chips
    mkdir -p ${fig_path}/${txt}/landau/chips
    mkdir -p ${fig_path}/${txt}/gaus/layers
    mkdir -p ${fig_path}/${txt}/exp/layers
    mkdir -p ${fig_path}/${txt}/landau/layers

    count=$(ls ${Calib_path}/${txt} -U1 | wc -l)
    echo -n 1 >> ../tmp/threshold_${txt}.txt
    echo -n " " >> ../tmp/threshold_${txt}.txt
    echo -n ${txt} >> ../tmp/threshold_${txt}.txt
    echo -n " " >> ../tmp/threshold_${txt}.txt
    echo -n ${count} >> ../tmp/threshold_${txt}.txt
    echo -n " " >> ../tmp/threshold_${txt}.txt

    for txtRoot in $( ls ${Calib_path}/${txt}/ECAL*.root ); do
      echo -n ${txtRoot} >> ../tmp/threshold_${txt}.txt
      echo -n " " >> ../tmp/threshold_${txt}.txt
    #  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
    done

    var=`cat ../tmp/threshold_${txt}.txt`
    # ./threshold_finding ../result/threshold/threshold_${txt}.root ${var} ${fig_path}
    sbatch -o jobs/test-%A.out --error="jobs/test-%A.err" execute_paralell.sl ../result/threshold/threshold_${txt}.root ${var} ${fig_path}
  done
}




#Calib_path=${data_dir}/${data_year}/2023/sps/calib/e-_noIntercept
Calib_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/decode/e-
Calib_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/decode/e-
#Calib_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/tmp/${year}/${beam}/calib/${date_tag} ############ for edited calibration

fig_path=../result/threshold/figures


paralell_jobs ${Calib_path_sps} ${fig_path}
paralell_jobs ${Calib_path_ps} ${fig_path}
