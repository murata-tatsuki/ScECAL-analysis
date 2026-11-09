#! /bin/bash
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$script_dir" || exit 1
slurm_log_dir=$(realpath -m -- "${script_dir}/../log/mip_new")
result_dir=$(realpath -m -- "${script_dir}/../result/mip_new")
mkdir -p -- "$slurm_log_dir" "$result_dir/figures/chips" "$result_dir/figures/compare"

source /megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO/run/global_config

#Decode_path=${data_dir}/${year}/sps/simpleCalib/mu-/100GeV
Decode_path=${data_dir}/${year}/sps/decode/mu-/100GeV
SSA_path=${data_dir}/${year}/sps/ssa/mu-/100GeV

: > "${slurm_log_dir}/mip_inputs.txt"
ii=0
#for txt in $( ls /home/jap/data2023/ECAL/decode/LEDCali/*dac* ); do
for i in $( ls ${Decode_path} ); do
  echo -n ${SSA_path}/${i} >> "${slurm_log_dir}/mip_inputs.txt"
  echo -n " " >> "${slurm_log_dir}/mip_inputs.txt"
  echo -n ${Decode_path}/${i} >> "${slurm_log_dir}/mip_inputs.txt"
  echo -n " " >> "${slurm_log_dir}/mip_inputs.txt"

#  echo ${ii} " " ${i}
#  let ii++

#  ./between_files_mip_analysis ../results/mip/${i} ${SSA_path}/${i} ${Decode_path}/${i}
done

var=`cat "${slurm_log_dir}/mip_inputs.txt"`
# ./MIP_calibration ../result/mip_new/mip.root ${var}


sbatch -o "${slurm_log_dir}/test-%A.out" --error="${slurm_log_dir}/test-%A.err" --chdir="$script_dir" "$script_dir/execute.sl" "$script_dir/MIP_calibration" "$result_dir/mip.root" ${var}
