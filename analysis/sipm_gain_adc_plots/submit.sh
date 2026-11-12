#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
label=full_grouped
output_label=
sample=
kind=led
concurrent=8
partition=all
memory=4G
walltime=08:00:00
adc_range=data
spectral=1
fft_reconstruction=1
verify=metadata
execute=0
while (( $# )); do
  case $1 in
    --execute) execute=1; shift; continue ;;
    --help)
      echo 'submit.sh --sample NAME|all [--kind led|electron|muon|all] [--label full_grouped]'
      echo '          [--output-label NAME] [--max-concurrent 8] [--partition all]'
      echo '          [--memory 4G] [--walltime 08:00:00] [--adc-range data|fit|full]'
      echo '          [--fft-reconstruction 0|1] [--spectral-labels 0|1] [--verify-inputs metadata|sha256] [--execute]'
      echo 'Default is preparation only; --execute submits sample x 30 layer tasks.'
      exit 0 ;;
  esac
  (( $# >= 2 )) || { echo "Missing value: $1" >&2; exit 2; }
  case $1 in
    --label) label=$2 ;; --output-label) output_label=$2 ;; --sample) sample=$2 ;;
    --kind) kind=$2 ;; --max-concurrent) concurrent=$2 ;; --partition) partition=$2 ;;
    --memory) memory=$2 ;; --walltime) walltime=$2 ;; --adc-range) adc_range=$2 ;;
    --fft-reconstruction) fft_reconstruction=$2 ;; --spectral-labels) spectral=$2 ;; --verify-inputs) verify=$2 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift 2
done
output_label=${output_label:-$label}
for value in "$label" "$output_label" "$sample"; do
  [[ $value =~ ^[A-Za-z0-9_.-]+$ && $value != . && $value != .. ]] || { echo 'Invalid or missing label/sample' >&2; exit 2; }
done
[[ $concurrent =~ ^[1-9][0-9]*$ && $concurrent -le 128 ]]
[[ $kind == led || $kind == electron || $kind == muon || $kind == all ]]
[[ $adc_range == data || $adc_range == fit || $adc_range == full ]]
[[ $spectral == 0 || $spectral == 1 ]]
[[ $fft_reconstruction == 0 || $fft_reconstruction == 1 ]]
[[ $verify == metadata || $verify == sha256 ]]
[[ -x $project/analysis/sipm_gain_adc_plots/Plot ]] || { echo 'Build Plot first' >&2; exit 2; }
input=$project/analysis/result/sipm_gain/$label/data
samples=()
for k in led electron muon; do
  [[ $kind == all || $kind == "$k" ]] || continue
  for d in "$input/$k/"*; do
    [[ -d $d/collection ]] || continue
    name=${d##*/}
    [[ $sample == all || $sample == "$name" ]] || continue
    [[ $name =~ ^[A-Za-z0-9_.-]+$ ]]
    samples+=("$name")
  done
done
(( ${#samples[@]} )) || { echo 'No matching samples' >&2; exit 2; }
log=$project/analysis/log/sipm_gain_adc_plots/$output_label
mkdir -p "$log"
submission=$(mktemp -d "$log/submission.XXXXXX")
tasks=$submission/tasks.tsv
for name in "${samples[@]}"; do
  for (( layer=0; layer<30; ++layer )); do printf '%s\t%d\n' "$name" "$layer"; done
done > "$tasks"
total=$(( ${#samples[@]} * 30 ))
printf 'Samples: %d; tasks: %d; concurrent: %d\nTask manifest: %s\n' "${#samples[@]}" "$total" "$concurrent" "$tasks"
args=(--label "$label" --output-label "$output_label" --adc-range "$adc_range" --spectral-labels "$spectral" --fft-reconstruction "$fft_reconstruction" --verify-inputs "$verify")
printf '%q ' "${args[@]}" > "$submission/plot_options.txt"
printf '\n' >> "$submission/plot_options.txt"
if (( ! execute )); then
  echo 'Preparation only. Add --execute to submit; nothing has been queued.'
  exit 0
fi
previous=
for (( offset=0; offset<total; offset+=1000 )); do
  count=$(( total-offset )); (( count <= 1000 )) || count=1000
  cmd=(sbatch --parsable --job-name=sipm_adc_plots --partition="$partition"
       --array="0-$((count-1))%$concurrent" --cpus-per-task=1 --mem="$memory" --time="$walltime"
       --chdir="$project" --output="$submission/%A_%a.out" --error="$submission/%A_%a.err")
  [[ -z $previous ]] || cmd+=(--dependency="afterok:$previous")
  cmd+=("$project/analysis/sipm_gain_adc_plots/plot_job.sh" "$tasks" "$offset" "${args[@]}")
  printf '%q ' "${cmd[@]}" >> "$submission/commands.sh"
  printf '\n' >> "$submission/commands.sh"
  submitted=$("${cmd[@]}")
  previous=${submitted%%;*}
  printf '%s\t%d\t%d\n' "$previous" "$offset" "$count" >> "$submission/jobs.tsv"
  printf 'Submitted %s: task rows %d..%d\n' "$previous" "$offset" "$((offset+count-1))"
done
