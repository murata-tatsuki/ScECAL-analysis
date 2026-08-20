#! /bin/bash
source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

paralell_jobs () {
  decode_path=$1
  fig_path=$2
  cog_rage=$3
  root_path=$4
  preffix=$5
  suffix=$6
  mkdir -p ${root_path}
  
  for txt in $( ls ${decode_path} ); do
    rm ../tmp/analyze_${suffix}_${txt}.txt
    touch ../tmp/analyze_${suffix}_${txt}.txt
    mkdir -p ${fig_path}/${txt}/raw

    count=$(ls ${decode_path}/${txt} -U1 | wc -l)
    echo -n 1 >> ../tmp/analyze_${suffix}_${txt}.txt
    echo -n " " >> ../tmp/analyze_${suffix}_${txt}.txt
    echo -n ${txt} >> ../tmp/analyze_${suffix}_${txt}.txt
    echo -n " " >> ../tmp/analyze_${suffix}_${txt}.txt
    echo -n ${count} >> ../tmp/analyze_${suffix}_${txt}.txt
    echo -n " " >> ../tmp/analyze_${suffix}_${txt}.txt

    for txtRoot in $( ls ${decode_path}/${txt}/${preffix}*.root ); do
      echo -n ${txtRoot} >> ../tmp/analyze_${suffix}_${txt}.txt
      echo -n " " >> ../tmp/analyze_${suffix}_${txt}.txt
    #  ./between_files_mip_analysis tmp/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
    done

    var=`cat ../tmp/analyze_${suffix}_${txt}.txt`
    # ./threshold_finding ../result/threshold/analyze_${txt}.root ${var} ${fig_path}
    sbatch -o jobs/singleEnergy/test-%A.out --error="jobs/singleEnergy/test-%A.err" execute.sl ./SingleEnergyAnalysis ${root_path}/${txt}_${cog_rage}mm.root ${var} ${cog_rage} 1 ${fig_path}
  done
}


mkdir -p jobs/singleEnergy

# HLG=HighGain
HLG=LowGain


#decode_path=${data_dir}/${data_year}/2023/sps/decode/e-_noIntercept
decode_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/decode/e-
decode_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/decode/e-
#decode_path=/megraid01/users/murata_t/scecal/ScECAL_BeamTest/tmp/${year}/${beam}/decode/${date_tag} ############ for edited decoderation

cog_rage=200
fig_path=../result/discrepancy/figures/data/${HLG}/${cog_rage}mm
root_path=../result/discrepancy/data/${HLG}/${cog_rage}mm
mkdir -p ${fig_path}
mkdir -p ${root_path}

paralell_jobs ${decode_path_sps} ${fig_path} ${cog_rage} ${root_path} ECAL data
paralell_jobs ${decode_path_ps} ${fig_path} ${cog_rage} ${root_path} ECAL data



thre=threshold

decode_path_simulation=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/decode/e-/sps/${thre}
fig_path=../result/discrepancy/figures/simulation/${HLG}/${thre}/${cog_rage}mm
root_path=../result/discrepancy/simulation/${HLG}/${thre}/${cog_rage}mm
mkdir -p ${fig_path}
mkdir -p ${root_path}
paralell_jobs ${decode_path_simulation} ${fig_path} ${cog_rage} ${root_path} e- sim




# mu
cog_rage=200
fig_path=../result/discrepancy/figures/data/${HLG}/mu/${cog_rage}mm
root_path=../result/discrepancy/data/${HLG}/mu/${cog_rage}mm
mkdir -p ${fig_path}
mkdir -p ${root_path}

decode_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/decode/mu-
decode_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/decode/mu-

paralell_jobs ${decode_path_sps} ${fig_path} ${cog_rage} ${root_path} ECAL data
paralell_jobs ${decode_path_ps} ${fig_path} ${cog_rage} ${root_path} ECAL data