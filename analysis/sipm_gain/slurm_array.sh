#!/bin/bash
set -euo pipefail
[[ $# == 4 ]] || { echo 'Usage: slurm_array.sh RUNTIME PROJECT TASKS SHA256' >&2; exit 2; }
GAIN_RUNTIME=$1
export SIPM_GAIN_PROJECT=$2
GAIN_TASKS=$3
GAIN_SHA=$4
[[ "$GAIN_RUNTIME" = /* && -x "$GAIN_RUNTIME/Analyze" ]] || { echo 'Invalid frozen runtime' >&2; exit 2; }
GAIN_ACTUAL=$(sha256sum -- "$GAIN_TASKS")
[[ "${GAIN_ACTUAL%% *}" == "$GAIN_SHA" ]] || { echo 'Task manifest changed' >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0
exec "$GAIN_RUNTIME/Analyze" task --tasks "$GAIN_TASKS" --index "${SLURM_ARRAY_TASK_ID:?}"
