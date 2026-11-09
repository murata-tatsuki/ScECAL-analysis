#!/bin/bash
set -euo pipefail
shopt -s nullglob
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
result_base=${RESULT_BASE:-"$script_dir/../result/cog_channel_scan"}
read -r -a cog_ranges <<< "${COG_RANGES:-0 20 10 5}"
read -r -a channel_ranges <<< "${CHANNEL_RANGES:-68.2 22.5}"
read -r -a energies <<< "${ENERGIES:-0.5 1 2 3 4 5 10 20 30 40 50 60 70 80 100 120 150 200 250}"
# Only the current e_ssa / tail-excluded-data products are accepted.
for mode in dataThre_cog_channel dataThre; do
  selected_cogs=("${cog_ranges[@]}")
  if [[ $mode == dataThre ]]; then
    [[ " ${cog_ranges[*]} " == *' 0 '* ]] || continue
    selected_cogs=(0)
  fi
  files=()
  for cog in "${selected_cogs[@]}"; do
    for channel in "${channel_ranges[@]}"; do
      condition="cog${cog}mm_channel${channel}mm"
      for sample in data mc; do
        selection="data/exclude/channel${channel}mm/${cog}mm"
        [[ $sample != mc ]] || selection="mc/threshold/e_ssa/channel${channel}mm/${cog}mm"
        for energy in "${energies[@]}"; do
          file="$result_base/$selection/${energy}GeV_${condition}.root"
          [[ -f $file ]] || { echo "Missing current result: $file" >&2; exit 1; }
          files+=("$file")
        done
      done
    done
  done
  figures="$result_base/figures/comparison/$mode"
  command=("$script_dir/MultiEnergyAnalysis" "$result_base/comparison_${mode}.root"
    "$((2*${#selected_cogs[@]}*${#channel_ranges[@]}))" "${#energies[@]}" "${files[@]}" "$figures"
    --require-e-ssa-data-tail)
  if [[ ${DRY_RUN:-0} == 1 ]]; then
    printf '%q ' "${command[@]}"; printf '\n'
  else
    mkdir -p "$figures"
    {
      TZ=Asia/Tokyo date --iso-8601=seconds
      printf 'Comparison: %s; data tail excluded; MC e_ssa/threshold tail kept\n' "$mode"
      printf 'Command: '; printf '%q ' "${command[@]}"; printf '\n'
    } >> "$result_base/run.log"
    "${command[@]}"
  fi
done
