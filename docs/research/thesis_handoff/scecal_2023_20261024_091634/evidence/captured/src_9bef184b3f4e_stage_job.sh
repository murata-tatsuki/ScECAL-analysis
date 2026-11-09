#!/usr/bin/env bash
set -euo pipefail
stage=${1:?stage required}
label=${2:?output label required}
source_label=${3:-full_grouped}
script_directory=$(cd -- "$(dirname -- "$0")" && pwd -P)
project_directory=$(cd -- "$script_directory/../.." && pwd -P)
cd "$project_directory"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0
case "$stage" in
  prepare)
    exec "$script_directory/Pool" --from-label "$source_label" --label "$label" --stage prepare --settings-policy common-within-beam ;;
  collect|merge)
    task=${SLURM_ARRAY_TASK_ID:?array task required}
    if [[ "$stage" == collect ]]; then group=$((task / 30)); layer=$((task % 30)); else group=$task; fi
    case "$group" in
      0) beam=ps; kind=electron ;;
      1) beam=ps; kind=muon ;;
      2) beam=sps; kind=electron ;;
      3) beam=sps; kind=muon ;;
      *) echo 'Invalid task group' >&2; exit 2 ;;
    esac
    args=(--from-label "$source_label" --label "$label" --stage "$stage" --beam "$beam" --kind "$kind")
    if [[ "$stage" == collect ]]; then args+=(--layer "$layer"); fi
    exec "$script_directory/Pool" "${args[@]}" ;;
  compare)
    exec ./analysis/sipm_gain/Analyze compare --label "$label" --stage numeric ;;
  plots)
    exec ./analysis/sipm_gain/Analyze compare --label "$label" --stage plots --layer "${SLURM_ARRAY_TASK_ID:?array task required}" ;;
  *) echo 'Unknown stage' >&2; exit 2 ;;
esac
