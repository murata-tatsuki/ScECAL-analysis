#!/usr/bin/env bash
set -euo pipefail
stage=${1:?Usage: submit_stage.sh prepare|collect|merge|fit|compare|plots [label] [source-label]}
label=${2:-full_grouped_temperature_pool}
source_label=${3:-full_grouped}
[[ "$label" =~ ^[A-Za-z0-9_-][A-Za-z0-9_.-]*$ && "$source_label" =~ ^[A-Za-z0-9_-][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid label' >&2; exit 2; }
script_directory=$(cd -- "$(dirname -- "$0")" && pwd -P)
project_directory=$(cd -- "$script_directory/../.." && pwd -P)
cd "$project_directory"
if [[ "$stage" == fit ]]; then
  exec ./analysis/sipm_gain/Analyze submit --manifest "analysis/result/sipm_gain/$label/inputs/manifest.tsv" --label "$label" --stage fit --channel-groups 30 --max-concurrent 32 --memory 8G --walltime 24:00:00 --execute
fi
args=()
case "$stage" in
  prepare) ;;
  collect) args+=(--array=0-119%32) ;;
  merge) args+=(--array=0-3%4) ;;
  compare) ;;
  plots) args+=(--array=0-29%32) ;;
  *) echo 'Unknown stage' >&2; exit 2 ;;
esac
log_directory="$project_directory/analysis/log/sipm_gain_pool/$label/$stage"
mkdir -p "$log_directory"
sbatch --parsable --job-name="sipm_pool_$stage" --partition=all "${args[@]}" --cpus-per-task=1 --mem=8G --time=24:00:00 --chdir="$project_directory" --output="$log_directory/%A_%a.out" --error="$log_directory/%A_%a.err" "$script_directory/stage_job.sh" "$stage" "$label" "$source_label"
