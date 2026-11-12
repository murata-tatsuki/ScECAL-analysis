#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
LOG=$(realpath -m -- "$HERE/../../log/channel_response/strip_position_response")
partition=all;memory=8G;limit=04:00:00;workers=2;dry=0;args=()
while (($#));do case $1 in
 --partition)partition=$2;shift 2;; --mem)memory=$2;shift 2;; --time)limit=$2;shift 2;;
 --workers)workers=$2;args+=("$1" "$2");shift 2;; --dry-run)dry=1;shift;;
 --help|-h)echo 'submit.sh [--partition all] [--mem 8G] [--time 04:00:00] [run.sh options]';bash "$HERE/run.sh" --help;exit 0;;
 *)args+=("$1");shift;;esac;done
bash "$HERE/run.sh" "${args[@]}" --dry-run
cmd=(sbatch --parsable --partition="$partition" --mem="$memory" --time="$limit" --cpus-per-task="$workers" --ntasks=1 --job-name=strip_position_response --output="$LOG/slurm-%j.out" --error="$LOG/slurm-%j.err" "$HERE/../execute.sl" bash "$HERE/run.sh" "${args[@]}")
printf '%q ' "${cmd[@]}";printf '\n'
((dry==0))||exit 0
mkdir -p "$LOG"
job=$("${cmd[@]}")
printf '%s\t%s\t' "$(date -Is)" "$job" >> "$LOG/submissions.tsv";printf '%q ' "${cmd[@]}" >> "$LOG/submissions.tsv";printf '\n' >> "$LOG/submissions.tsv"
printf 'Submitted %s; logs %s\n' "$job" "$LOG"
