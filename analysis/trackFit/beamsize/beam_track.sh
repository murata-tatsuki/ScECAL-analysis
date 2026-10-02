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
  
    sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" execute.sl ./beamsize_track ${outfile} ${dat_name} ${fig_path}
  done

  # var=`cat ../../tmp/beamsize_${beam}.txt`
  # outfile=${outPath}/${beam}/all.root
  # fig_path=${outPath}/figures/${beam}/all
  # mkdir -p ${fig_path}
  # sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" execute.sl ./beamsize_track ${outfile} ${var} ${fig_path}
}


mkdir -p "${slurm_log_dir}"


outPath_=../../result/trackFit/beamsize_track

track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit/mu-/100GeV
paralell_jobs ${track_path_} sps ${outPath_} mu-

track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit/mu-/10GeV/
paralell_jobs ${track_path_} ps ${outPath_} mu-

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/ssa/e-
# for enregy in $(ls ${track_path_})
# do
#   paralell_jobs ${track_path_}/${enregy} sps ${outPath_} e-
# done

# track_path_=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/ssa/e-
# for enregy in $(ls ${track_path_})
# do
#   paralell_jobs ${track_path_}/${enregy} ps ${outPath_} e-
# done


# outPath=../../result/trackFit3D/efficiency

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/trackFit3D/mu-/100GeV
# paralell_jobs ${track_path} sps ${outPath}

# track_path=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/trackFit3D/mu-/10GeV/
# paralell_jobs ${track_path} ps ${outPath}



