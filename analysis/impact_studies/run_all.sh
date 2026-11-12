#!/bin/bash
# Direct (non-batch) full-energy run, followed by an exact-list summary.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "$HERE/energies.sh"
beam=auto dry=0
run_args=() summary_args=()
while (($#));do
 case $1 in
  -h|--help)
   echo 'Usage: bash run_all.sh [run.sh options]'
   echo 'Runs all 19 energies sequentially and summarizes after all succeed.'
   echo 'For parallel cluster jobs: bash submit.sh --energy all --summarize'
   exit 0;;
  --dry-run) dry=1;shift;;
  --no-plots) run_args+=("$1");shift;;
  --*)
   (($#>=2))||{ echo "Missing value for $1" >&2;exit 1; }
   case $1 in
    --energy) [[ $2 == all ]]||{ echo 'run_all.sh only accepts --energy all' >&2;exit 1; };;
    --beam) beam=$2;run_args+=("$1" "$2");;
    --study|--mc-tag|--label|--result-dir) run_args+=("$1" "$2");summary_args+=("$1" "$2");;
    *) run_args+=("$1" "$2");;
   esac
   shift 2;;
  *) echo "Unexpected argument: $1" >&2;exit 1;;
 esac
done
case $beam in auto) energies_csv=$all_energies;; ps) energies_csv=$ps_energies;; sps) energies_csv=$sps_energies;; *) echo 'Invalid beam' >&2;exit 1;; esac
IFS=, read -r -a energies <<< "$energies_csv"
# Check every input before any processing starts.
for energy in "${energies[@]}";do
 bash "$HERE/run.sh" "${run_args[@]}" --energy "$energy" --dry-run
done
printf 'After all %s energies succeed: ' "${#energies[@]}"
printf '%q ' bash "$HERE/summarize.sh" "${summary_args[@]}" --energies "$energies_csv"
printf '\n'
((!dry))||exit 0
for energy in "${energies[@]}";do
 bash "$HERE/run.sh" "${run_args[@]}" --energy "$energy"
done
bash "$HERE/summarize.sh" "${summary_args[@]}" --energies "$energies_csv"
