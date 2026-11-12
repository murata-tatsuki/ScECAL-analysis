#!/usr/bin/env bash
set -euo pipefail
if (( $# < 2 )); then
  echo 'Usage: plot_job.sh TASKS.tsv OFFSET [Plot options...]' >&2
  exit 2
fi
task_file=$1
task_offset=$2
shift 2
# Slurm copies this script into /var/spool/slurmd; BASH_SOURCE is not the
# repository path. submit.sh supplies the project directory via --chdir.
project=$(pwd -P)
[[ -x $project/analysis/sipm_gain_adc_plots/Plot ]] || {
  echo "ERROR: Plot not found under $project; submit with --chdir set to the project directory" >&2
  exit 2
}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0
: "${SLURM_ARRAY_TASK_ID:?This script runs as a Slurm array task}"
[[ $task_offset =~ ^[0-9]+$ && $SLURM_ARRAY_TASK_ID =~ ^[0-9]+$ ]]
line=$((task_offset + SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r sample layer < <(sed -n "${line}p" "$task_file")
[[ $sample =~ ^[A-Za-z0-9_.-]+$ && $layer =~ ^[0-9]+$ ]]
exec ./analysis/sipm_gain_adc_plots/Plot --project "$project" --sample "$sample" --layer "$layer" "$@"
