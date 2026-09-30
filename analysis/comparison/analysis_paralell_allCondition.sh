#!/bin/bash
set -euo pipefail
shopt -s nullglob
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"
result_base="$script_dir/../result/comparison"
log_dir="$script_dir/job"

paralell_jobs() {
  local calib_path=$1 fig_path=$2 cog_range=$3 channel_range=$4 root_path=$5 prefix=$6
  local condition="cog${cog_range}mm_channel${channel_range}mm"
  local energy_dir energy_name
  local -a files command
  if [[ ! -d "$calib_path" ]]; then
    echo "Calibration directory not found: $calib_path" >&2
    return 1
  fi
  for energy_dir in "$calib_path"/*/; do
    energy_name=$(basename "$energy_dir")
    files=("$energy_dir"/"$prefix"*.root)
    if ((${#files[@]} == 0)); then
      echo "No calibration ROOT files: $energy_dir" >&2
      continue
    fi
    command=(sbatch -o "$log_dir/test-%A.out"
      --error="$log_dir/test-%A.err"
      "$script_dir/execute_paralell.sl" "$script_dir/SingleEnergyAnalysis"
      "$root_path/${energy_name}_${condition}.root"
      1 "$energy_name" "${#files[@]}" "${files[@]}" "$cog_range" "$channel_range" 1 "$fig_path")
    if [[ ${DRY_RUN:-0} == 1 ]]; then
      printf '%q ' "${command[@]}"
      printf '\n'
    else
      mkdir -p "$root_path" "$fig_path/$energy_name/raw" "$log_dir"
      "${command[@]}"
    fi
  done
}

# x・y の半幅 [mm]。CoG は layer 9・10 の両方に適用、0 は CoG カットなし。
# compare_ds_allCondition.sh と同じ配列にする（小数も可）。
# 凡例は上から no cut → 20 → 10 → 5 mm（下ほど厳しいカット）。
cog_ranges=(0 20 10 5)
channel_ranges=(68.2 22.5)
# thres=(threshold threshold0 threshold50)
thres=(threshold)
calib_path_sps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/sps/calib/e-
calib_path_ps=/megraid01/users/data_beamtest/ECAL_data/analysed/2023/ps/calib/e-

for cog_range in "${cog_ranges[@]}"; do
  for channel_range in "${channel_ranges[@]}"; do
    condition="cog${cog_range}mm_channel${channel_range}mm"
    fig_path="$result_base/figures/data/$condition"
    root_path="$result_base/data/$condition"
    paralell_jobs "$calib_path_sps" "$fig_path" "$cog_range" "$channel_range" "$root_path" ECAL
    paralell_jobs "$calib_path_ps" "$fig_path" "$cog_range" "$channel_range" "$root_path" ECAL

    for thre in "${thres[@]}"; do
      calib_path_simulation="/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/sps/$thre"
      fig_path="$result_base/figures/simulation/$thre/$condition"
      root_path="$result_base/simulation/$thre/$condition"
      paralell_jobs "$calib_path_simulation" "$fig_path" "$cog_range" "$channel_range" "$root_path" e-
    done
  done
done
