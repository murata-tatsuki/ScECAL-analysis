#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energies=5,20,100 partition=all dry=0 label=initial
submission_dir=$(pwd -P)
args=()
while (($#));do case $1 in
 --energies) energies=$2;shift 2;; --partition) partition=$2;shift 2;;
 --dry-run) dry=1;shift;; --label) label=$2;args+=("$1" "$2");shift 2;;
 --plots-only|--intervals-only|--maps-only|--style-only) args+=("$1");shift;; --help|-h) echo 'submit.sh [--energies all|5,20,100] [--partition all] [--dry-run] [run.sh options, including --track-cut off|on|both|v2 and track thresholds]';exit;;
 --energy) echo 'Use --energies here' >&2;exit 2;;
 *) (($#>=2))||exit 2;args+=("$1" "$2");shift 2;;esac;done
[[ $label =~ ^[A-Za-z0-9_-]+$ ]]||exit 2
[[ $energies != all ]]||energies=0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250
IFS=, read -ra list <<< "$energies"
declare -A seen=()
for energy in "${list[@]}";do [[ $energy =~ ^(0\.5|1|2|3|4|5|10|20|30|40|50|60|70|80|100|120|150|200|250)$ && ! -v seen[$energy] ]]||exit 2;seen[$energy]=1;done
logs=$(realpath -m "$HERE/../log/shower_leak/$label")
if ((!dry));then mkdir -p "$logs";record=$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv;fi
for energy in "${list[@]}";do
 cmd=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1 --export=ALL --chdir="$submission_dir" --job-name="shower_leak_${energy}" --output="$logs/${energy}GeV-%j.out" --error="$logs/${energy}GeV-%j.err" "$HERE/../channel_response/execute.sl" bash "$HERE/run.sh" --energy "$energy" "${args[@]}")
 printf '%q ' "${cmd[@]}";printf '\n'
 if ((!dry));then job=$("${cmd[@]}");printf '%s\t%s\n' "$energy" "$job" >> "$record";fi
done
