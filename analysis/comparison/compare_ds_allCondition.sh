#!/bin/bash
set -euo pipefail
shopt -s nullglob
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir"
result_base="$script_dir/../result/comparison"

# 単一エネルギーの全ジョブが完了してから実行する。
thre=threshold
path_suffix=dataThre_cog_channel
# analysis_paralell_allCondition.sh と同じ配列にする。CoG 0 はカットなし。
# 凡例は上から no cut → 20 → 10 → 5 mm（下ほど厳しいカット）。
cog_ranges=(0 20 10 5)
channel_ranges=(68.2 22.5)
fig_path="$result_base/figures/comparison/$path_suffix"
files=()
reference_energies=()
count=0

for cog_range in "${cog_ranges[@]}"; do
  for channel_range in "${channel_ranges[@]}"; do
    condition="cog${cog_range}mm_channel${channel_range}mm"
    for sample in data "simulation/$thre"; do
      sample_files=("$result_base/$sample/$condition/"*"_${condition}.root")
      if ((${#sample_files[@]} == 0)); then
        echo "No input ROOT files: $result_base/$sample/$condition" >&2
        exit 1
      fi
      energies=()
      for file in "${sample_files[@]}"; do
        name=${file##*/}
        energies+=("${name%_${condition}.root}")
      done
      if ((count == 0)); then
        count=${#sample_files[@]}
        reference_energies=("${energies[@]}")
      elif [[ "${energies[*]}" != "${reference_energies[*]}" ]]; then
        echo "Energy files do not match: $sample, $condition" >&2
        exit 1
      fi
      files+=("${sample_files[@]}")
    done
  done
done

command=("$script_dir/MultiEnergyAnalysis" "$result_base/comparison_${path_suffix}.root"
  "$((2 * ${#cog_ranges[@]} * ${#channel_ranges[@]}))" "$count" "${files[@]}" "$fig_path")
if [[ ${DRY_RUN:-0} == 1 ]]; then
  printf '%q ' "${command[@]}"
  printf '\n'
else
  mkdir -p "$fig_path"
  "${command[@]}"
fi
