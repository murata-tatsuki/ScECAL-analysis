#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
LOG=$(realpath -m -- "$HERE/../../log/channel_response/time_stability")
partition=all;memory=8G;limit=02:00:00;cpus=1;dry=0;args=()
while (($#));do case $1 in
 --cpus)cpus=$2;shift 2;; --partition)partition=$2;shift 2;; --mem)memory=$2;shift 2;; --time)limit=$2;shift 2;;
 --dry-run)dry=1;shift;; --help|-h)printf 'submit.sh [--partition all] [--mem 8G] [--time 02:00:00] [--cpus 1] [run.sh options]\nMultiple fit processes with --cpus; parent is the sole ROOT/TSV writer.\n';exec bash "$HERE/run.sh" --help;;
 *)args+=("$1");shift;;esac;done
[[ $cpus =~ ^[0-9]+$ ]]&&((cpus>=1&&cpus<=33))||{ printf 'Invalid CPU count\n' >&2;exit 1; }
SLURM_CPUS_PER_TASK=$cpus bash "$HERE/run.sh" "${args[@]}" --dry-run
cmd=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task="$cpus" --mem="$memory" --time="$limit" --job-name=time_stability --chdir="$HERE" --output="$LOG/slurm-%j.out" --error="$LOG/slurm-%j.err" "$HERE/../execute.sl" bash "$HERE/run.sh" "${args[@]}")
if ((dry));then printf '%q ' "${cmd[@]}";printf '\n';exit 0;fi
mkdir -p -- "$LOG"
job=$("${cmd[@]}")
printf '%s\t%s\t' "$(date -Is)" "$job" >> "$LOG/submissions.tsv";printf '%q ' "${cmd[@]}" >> "$LOG/submissions.tsv";printf '\n' >> "$LOG/submissions.tsv"
printf 'Submitted %s; logs %s/slurm-%s.{out,err}\n' "$job" "$LOG" "${job%%;*}"
