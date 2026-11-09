#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
LOG=$(realpath -m -- "$HERE/../../log/channel_response/muon_mip_response")
partition=all
memory=16G
time_limit=24:00:00
workers=4
args=()
while (($#));do case $1 in
 --partition)partition=$2;shift 2;; --mem)memory=$2;shift 2;; --time)time_limit=$2;shift 2;;
 --workers)workers=$2;args+=("$1" "$2");shift 2;;
 --help|-h)printf 'Usage: bash submit.sh [--partition all] [--mem 16G] [--time 24:00:00] [run.sh options]\n';bash "$HERE/run.sh" --help;exit 0;;
 *)args+=("$1");shift;;esac;done
for arg in "${args[@]}";do if [[ $arg == --dry-run ]];then exec bash "$HERE/run.sh" "${args[@]}";fi;done
bash "$HERE/run.sh" "${args[@]}" --dry-run
mkdir -p -- "$LOG"
job=$(sbatch --parsable --partition="$partition" --mem="$memory" --time="$time_limit" --cpus-per-task="$workers" --ntasks=1 --job-name=muon_mip_response --output="$LOG/slurm-%j.out" --error="$LOG/slurm-%j.err" "$HERE/../execute.sl" bash "$HERE/run.sh" "${args[@]}")
printf '%s\t%s\n' "$(date -Is)" "$job" >> "$LOG/submissions.tsv"
printf 'Submitted Slurm job %s\nLogs: %s/slurm-%s.{out,err}\n' "$job" "$LOG" "${job%%;*}"
