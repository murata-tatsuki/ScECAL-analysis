#!/bin/bash
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --job-name=cog_channel_scan
#SBATCH --output=/home/murata_t/ScECAL_BeamTest/analysis/log/cog_channel_scan/slurm-%j.out
#SBATCH --error=/home/murata_t/ScECAL_BeamTest/analysis/log/cog_channel_scan/slurm-%j.err
set -euo pipefail
ulimit -c 0
if [[ ${1:-} == --manifest ]]; then
  [[ $# == 2 && ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || { echo 'Invalid array manifest invocation' >&2; exit 1; }
  line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$2")
  [[ -n $line ]] || { echo 'Empty array task' >&2; exit 1; }
  IFS=$'\t' read -r -a command <<< "$line"
  printf 'Command: '; printf '%q ' "${command[@]}"; printf '\n'
  exec "${command[@]}"
fi
exec "$@"
