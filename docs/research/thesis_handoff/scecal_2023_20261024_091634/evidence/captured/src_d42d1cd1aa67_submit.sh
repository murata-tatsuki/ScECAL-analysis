#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
LOG=$(realpath -m -- "$HERE/../log/mip_calibration")
partition=all;memory=16G;time_limit=24:00:00;workers=4;args=()
while (($#));do case $1 in
 --partition)partition=$2;shift 2;; --mem)memory=$2;shift 2;; --time)time_limit=$2;shift 2;;
 --workers)workers=$2;args+=("$1" "$2");shift 2;;
 --help|-h)printf 'Usage: submit.sh [--partition all] [--mem 16G] [--time 24:00:00] [run.sh options]\n';bash "$HERE/run.sh" --help;exit 0;;
 *)args+=("$1");shift;;esac;done
for arg in "${args[@]}";do [[ $arg != --dry-run ]]||exec bash "$HERE/run.sh" "${args[@]}";done
bash "$HERE/run.sh" "${args[@]}" --dry-run
mkdir -p -- "$LOG"
job=$(sbatch --parsable --partition="$partition" --mem="$memory" --time="$time_limit" --cpus-per-task="$workers" --ntasks=1 --job-name=mip_calibration --output="$LOG/slurm-%j.out" --error="$LOG/slurm-%j.err" "$HERE/../channel_response/execute.sl" bash "$HERE/run.sh" "${args[@]}")
{ printf '%s\t%s\t' "$(date -Is)" "$job";printf '%q ' "$HERE/run.sh" "${args[@]}";printf '\n'; } >> "$LOG/submissions.tsv"
printf 'Submitted: %s\nLogs: %s/slurm-%s.{out,err}\n' "$job" "$LOG" "${job%%;*}"
