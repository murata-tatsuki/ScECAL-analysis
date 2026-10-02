#! /bin/bash
slurm_log_dir=$(realpath -m -- "$(dirname -- "${BASH_SOURCE[0]}")/../../log/trackFit/beamsize")
mkdir -p -- "$slurm_log_dir"


paralell_jobs () {
  track_path=$1
  beam=$2
  outPath=$3
  particle=$4

  rm ../../tmp/beamsize_${beam}.txt
  touch ../../tmp/beamsize_${beam}.txt
  
  for dat_name in $(ls ${track_path}/*.root)
  do
    prefix=`basename $dat_name`
    mkdir -p ${outPath}/${beam}/${particle}
    outfile=${outPath}/${beam}/${particle}/${prefix}
    prefix="${prefix:0:-5}"
    fig_path=${outPath}/figures/${beam}/${particle}/${prefix}
    mkdir -p ${fig_path}

    echo -n ${dat_name} >> ../../tmp/beamsize_${beam}.txt
    echo -n " " >> ../../tmp/beamsize_${beam}.txt
  
    sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" execute.sl ./beamsize_ssa ${outfile} ${dat_name} ${fig_path}
  done

  var=`cat ../../tmp/beamsize_${beam}.txt`
  energy=`basename $track_path`
  outfile=${outPath}/${beam}/${particle}/${energy}.root
  fig_path=${outPath}/figures/${beam}/${particle}/${energy}
  mkdir -p ${fig_path}
  # sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" execute.sl ./beamsize_ssa ${outfile} ${var} ${fig_path}
}


mkdir -p "${slurm_log_dir}"


outPath_=../../result/trackFit/beamsize_ssa

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/ssaCalib/mu-/100GeV
# paralell_jobs ${track_path_} sps ${outPath_} mu-

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/ssaCalib/mu-/10GeV/
# paralell_jobs ${track_path_} ps ${outPath_} mu-

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/ssaCalib/e-
# for enregy in $(ls ${track_path_})
# do
#   paralell_jobs ${track_path_}/${enregy} sps ${outPath_} e-
# done

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/ssaCalib/e-
# for enregy in $(ls ${track_path_})
# do
#   paralell_jobs ${track_path_}/${enregy} ps ${outPath_} e-
# done







## simulation

track_path_=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/ssaCalib/mu-/mu_track/threshold/100GeV
paralell_jobs ${track_path_} simulation ${outPath_} mu-

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/ssaCalib/mu-/10GeV/
# paralell_jobs ${track_path_} ps ${outPath_} mu-
