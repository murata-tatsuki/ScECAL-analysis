#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energies=100 partition=all label=baseline dry=0
args=()
while (($#)); do
 case "$1" in
 --energies) energies=$2;shift 2;;
 --partition) partition=$2;shift 2;;
 --label) label=$2;args+=("$1" "$2");shift 2;;
 --dry-run) dry=1;shift;;
 --plots-only) args+=("$1");shift;;
 --help|-h) echo 'submit.sh [--energies 100|20,100|all] [--partition all] [--dry-run] [run.sh options]';exit 0;;
 --energy) echo 'Use --energies in submit.sh' >&2;exit 2;;
 *) (($#>=2)) || exit 2;args+=("$1" "$2");shift 2;;
 esac
done
[[ $label =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
[[ $energies != all ]] || energies=0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250
IFS=, read -ra list <<< "$energies"
declare -A seen=()
for energy in "${list[@]}"; do
 [[ $energy =~ ^(0\.5|1|2|3|4|5|10|20|30|40|50|60|70|80|100|120|150|200|250)$ && ! -v seen[$energy] ]] || exit 2
 seen[$energy]=1
done
logs=$(realpath -m "$HERE/../log/ssa_shower_shape/$label")
if ((!dry)); then mkdir -p "$logs";record="$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv";fi
for energy in "${list[@]}"; do
 cmd=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1 --export=ALL --job-name="ssa_shape_${energy}" --output="$logs/${energy}GeV-%j.out" --error="$logs/${energy}GeV-%j.err" "$HERE/../channel_response/execute.sl" bash "$HERE/run.sh" --energy "$energy" "${args[@]}")
 printf '%q ' "${cmd[@]}";printf '\n'
 if ((!dry)); then job=$("${cmd[@]}");printf '%s\t%s\n' "$energy" "$job" >> "$record";fi
done
