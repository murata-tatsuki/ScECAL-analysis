#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energies=100 partition=all label=beam_scan condition=point phase=scan dry=0
cuts=allChannel cut_option=0 args=()
while (($#)); do
 case "$1" in
 --energies) energies=$2;shift 2;;
 --partition) partition=$2;shift 2;;
 --label) label=$2;args+=("$1" "$2");shift 2;;
 --condition) condition=$2;args+=("$1" "$2");shift 2;;
 --phase) phase=$2;args+=("$1" "$2");shift 2;;
 --threshold-cut|--threshold-cuts)
   ((cut_option==0)) || { echo 'Specify threshold-cut OR threshold-cuts once' >&2;exit 2; }
   cuts=$2;cut_option=1;shift 2;;
 --dry-run) dry=1;shift;;
 --plots-only) args+=("$1");shift;;
 --help|-h)
   echo 'submit.sh [--energies 100|40,100] [--threshold-cuts allChannel,fitValid,precision,chi2_10,chi2_5,chi2_3] [--partition all] [--dry-run] [run.sh options]';exit 0;;
 --energy) echo 'Use --energies' >&2;exit 2;;
 *) (($#>=2)) || exit 2;args+=("$1" "$2");shift 2;;
 esac
done
[[ $label =~ ^[A-Za-z0-9_-]+$ && $condition =~ ^[A-Za-z0-9_-]+$ && $phase =~ ^(scan|validation|all)$ ]] || exit 2
[[ $energies != ,* && $energies != *, && $energies != *,,* && $cuts != ,* && $cuts != *, && $cuts != *,,* ]] || exit 2
IFS=, read -ra energy_list <<< "$energies"
IFS=, read -ra cut_list <<< "$cuts"
declare -A seen_energy=() seen_cut=()
for cut in "${cut_list[@]}"; do
 [[ $cut =~ ^[A-Za-z0-9_-]+$ && ! -v seen_cut[$cut] ]] || exit 2
 seen_cut[$cut]=1
done
for energy in "${energy_list[@]}"; do
 [[ $energy =~ ^[0-9]+(\.[0-9]+)?$ && ! -v seen_energy[$energy] ]] || exit 2
 seen_energy[$energy]=1
 for cut in "${cut_list[@]}"; do
  bash "$HERE/run.sh" --energy "$energy" --threshold-cut "$cut" "${args[@]}" --dry-run >/dev/null
 done
done
logs=$(realpath -m "$HERE/../log/strip_stage_response/$label/$phase/$condition")
if ((!dry)); then
 mkdir -p "$logs";record="$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv"
 printf 'energy\tthreshold_cut\tjob_id\n' > "$record"
fi
for energy in "${energy_list[@]}"; do
 for cut in "${cut_list[@]}"; do
  if ((!dry)); then mkdir -p "$logs/$cut";fi
  cmd=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=2 --mem=8G --time=04:00:00 --export=ALL --job-name="strip_stages_${condition}_${cut}_${energy}" --output="$logs/$cut/${energy}GeV-%j.out" --error="$logs/$cut/${energy}GeV-%j.err" "$HERE/../channel_response/execute.sl" bash "$HERE/run.sh" --energy "$energy" --threshold-cut "$cut" "${args[@]}")
  printf '%q ' "${cmd[@]}";printf '\n'
  if ((!dry)); then job=$("${cmd[@]}");printf '%s\t%s\t%s\n' "$energy" "$cut" "$job" >> "$record";fi
 done
done
