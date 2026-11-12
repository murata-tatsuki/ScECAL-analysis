#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energies=100 partition=all label=beam_scan condition=default phase=scan dry=0
args=()
while (($#)); do
 case "$1" in
 --energies) energies=$2;shift 2;;
 --partition) partition=$2;shift 2;;
 --label) label=$2;args+=("$1" "$2");shift 2;;
 --condition) condition=$2;args+=("$1" "$2");shift 2;;
 --phase) phase=$2;args+=("$1" "$2");shift 2;;
 --dry-run) dry=1;shift;;
 --plots-only|--threshold-diagnostics-only|--skip-threshold-diagnostics) args+=("$1");shift;;
 --help|-h) echo 'submit.sh [--energies 100|20,100|all] [--partition all] [--dry-run] [run.sh options]';exit 0;;
 --energy) echo 'Use --energies' >&2;exit 2;;
 *) (($#>=2)) || exit 2;args+=("$1" "$2");shift 2;;
 esac
done
[[ $label =~ ^[A-Za-z0-9_-]+$ && $condition =~ ^[A-Za-z0-9_-]+$ && $phase =~ ^(scan|validation|all)$ ]] || exit 2
[[ $energies != all ]] || energies=0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250
IFS=, read -ra list <<< "$energies"
declare -A seen=()
for energy in "${list[@]}"; do
 [[ $energy =~ ^(0\.5|1|2|3|4|5|10|20|30|40|50|60|70|80|100|120|150|200|250)$ && ! -v seen[$energy] ]] || exit 2
 seen[$energy]=1
 # Read-only preflight checks inputs/options and refuses output collisions.
 bash "$HERE/run.sh" --energy "$energy" "${args[@]}" --dry-run >/dev/null
done
logs=$(realpath -m "$HERE/../log/beam_quantitative/$label/$phase/$condition")
if ((!dry)); then mkdir -p "$logs";record="$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv";fi
for energy in "${list[@]}"; do
 cmd=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=2 --mem=8G --time=04:00:00 --export=ALL --job-name="beam_quant_${energy}" --output="$logs/${energy}GeV-%j.out" --error="$logs/${energy}GeV-%j.err" "$HERE/../channel_response/execute.sl" bash "$HERE/run.sh" --energy "$energy" "${args[@]}")
 printf '%q ' "${cmd[@]}";printf '\n'
 if ((!dry)); then job=$("${cmd[@]}");printf '%s\t%s\n' "$energy" "$job" >> "$record";fi
done
