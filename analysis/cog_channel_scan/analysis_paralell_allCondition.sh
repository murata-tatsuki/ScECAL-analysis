#!/bin/bash
set -euo pipefail
shopt -s nullglob
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
result_base=${RESULT_BASE:-"$script_dir/../result/cog_channel_scan"}
log_dir=${LOG_DIR:-"$script_dir/../log/cog_channel_scan"}
tail_results=${TAIL_RESULTS:-"$script_dir/../result/channel_response/hg_lg_tail/threshold/e_ssa"}
mc_base=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold
data_base=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
read -r -a cog_ranges <<< "${COG_RANGES:-0 20 10 5}"
read -r -a channel_ranges <<< "${CHANNEL_RANGES:-68.2 22.5}"
read -r -a energies <<< "${ENERGIES:-0.5 1 2 3 4 5 10 20 30 40 50 60 70 80 100 120 150 200 250}"
read -r -a samples <<< "${SAMPLES:-data mc}"
max_parallel=${MAX_PARALLEL:-24}
[[ $max_parallel =~ ^[1-9][0-9]*$ ]] || { echo 'MAX_PARALLEL must be positive' >&2; exit 1; }
run_id=$(TZ=Asia/Tokyo date +%Y%m%dT%H%M%S)-$$
manifest="$log_dir/tasks_${run_id}.tsv"
commands=()
outputs=()
figure_dirs=()
for cog in "${cog_ranges[@]}"; do
  for channel in "${channel_ranges[@]}"; do
    condition="cog${cog}mm_channel${channel}mm"
    for sample in "${samples[@]}"; do
      for energy in "${energies[@]}"; do
        if [[ $sample == data ]]; then
          beam=sps
          [[ $energy == 0.5 || $energy == [1-5] ]] && beam=ps
          input_dir="$data_base/$beam/calib/e-/${energy}GeV"
          files=("$input_dir"/ECAL*.root)
          selection="data/exclude/channel${channel}mm/${cog}mm"
          label=data
          options=(--exclude-tail-events --tail-results "$tail_results")
          [[ -f "$tail_results/${energy}GeV/tail_study.root" ]] || { echo "Missing tail results: $energy GeV" >&2; exit 1; }
        elif [[ $sample == mc ]]; then
          input_dir="$mc_base/${energy}GeV"
          files=("$input_dir"/*.root)
          selection="mc/threshold/e_ssa/channel${channel}mm/${cog}mm"
          label=e_ssa
          options=(--keep-tail-events)
        else
          echo "Unknown sample: $sample" >&2; exit 1
        fi
        ((${#files[@]})) || { echo "No calibration inputs: $input_dir" >&2; exit 1; }
        output="$result_base/$selection/${energy}GeV_${condition}.root"
        figures="$result_base/figures/$selection"
        command=("$script_dir/SingleEnergyAnalysis" "$output" 1 "${energy}GeV" "${#files[@]}"
          "${files[@]}" "$cog" "$channel" 1 "$figures" --sample-label "$label"
          --skip-channel-plots "${options[@]}")
        if [[ ${DRY_RUN:-0} == 1 ]]; then
          printf '%q ' "${command[@]}"; printf '\n'
        else
          printf -v line '%s\t' "${command[@]}"
          commands+=("${line%$'\t'}")
          outputs+=("$output")
          figure_dirs+=("$figures/${energy}GeV")
        fi
      done
    done
  done
done
[[ ${DRY_RUN:-0} != 1 ]] || exit 0
((${#commands[@]})) || { echo 'No tasks' >&2; exit 1; }
active=$(squeue -h -u "$USER" -n cog_channel_scan -o '%i')
[[ -z $active ]] || { echo "cog_channel_scan jobs still active: $active" >&2; exit 1; }
mkdir -p "$log_dir" "$result_base/inputs"
for output in "${outputs[@]}"; do mkdir -p "${output%/*}"; done
mkdir -p "${figure_dirs[@]}"
printf '%s\n' "${commands[@]}" > "$manifest"
sha256sum "$script_dir/SingleEnergyAnalysis.cc" "$script_dir/MultiEnergyAnalysis.cc" \
  "$script_dir/../resolution/HLTailSelection.hh" "$script_dir/"*.sh "$script_dir/execute_paralell.sl" \
  "$script_dir/SingleEnergyAnalysis" "$script_dir/MultiEnergyAnalysis" > "$result_base/inputs/source_checksums.txt"
job=$(sbatch --parsable --array="0-$((${#commands[@]}-1))%$max_parallel" --mem=4G \
  --output="$log_dir/single-%A_%a.out" --error="$log_dir/single-%A_%a.err" \
  "$script_dir/execute_paralell.sl" --manifest "$manifest")
printf 'job_id\ttasks\tmanifest\n%s\t%s\t%s\n' "$job" "${#commands[@]}" "$manifest" > "$log_dir/submission_${run_id}.tsv"
{
  TZ=Asia/Tokyo date --iso-8601=seconds
  printf 'Command: bash %q\n' "$script_dir/analysis_paralell_allCondition.sh"
  printf 'MC calibration: %s\nData calibration: %s/{ps,sps}/calib/e-\n' "$mc_base" "$data_base"
  printf 'Data tail: exclude (z<-5, >=3 distinct channels); MC tail: keep\nTail results: %s\n' "$tail_results"
  printf 'CoG: %s\nChannel: %s\nEnergies: %s\nSamples: %s\n' "${cog_ranges[*]}" "${channel_ranges[*]}" "${energies[*]}" "${samples[*]}"
  printf 'Job: %s\nTask manifest: %s\n' "$job" "$manifest"
} >> "$result_base/run.log"
printf 'Submitted %s tasks as array %s (max concurrent: %s)\nManifest: %s\n' "${#commands[@]}" "$job" "$max_parallel" "$manifest"
