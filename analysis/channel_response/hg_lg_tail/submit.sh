#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
submission_dir=$(pwd -P)
energies=100 queue=${LSF_QUEUE:-s} logs=$HERE/../../log/channel_response/hg_lg_tail dry=0 time_limit='' memory=''
args=()
while (($#));do
 case $1 in
  -h|--help) cat <<HELP
Usage: bash submit.sh [--energy VALUE|all] [--energies CSV] [run.sh options]
Default: 100 GeV only. One LSF job per energy, all channels/events by default.
--energy all includes PS 0.5,1,2,3,4,5 and SPS 10,20,30,40,50,60,70,80,100,120,150,200,250.
--queue NAME (default: LSF_QUEUE or s), --log-dir PATH (default: $HERE/../../log/channel_response/hg_lg_tail),
--mem MB (LSF -M/-R), --time MINUTES|HH:MM (LSF -W), --dry-run.
--dry-run validates inputs and prints bsub commands; it creates/submits nothing.
Status: bjobs
HELP
   exit 0;;
  --dry-run) dry=1;shift;;
  --*) (($#>=2))||exit 1;case $1 in --energy|--energies) energies=$2;;--queue) queue=$2;;--partition) echo 'Use --queue for KEKCC LSF' >&2;exit 1;;--log-dir) logs=$2;;--mem) memory=$2;;--time) time_limit=$2;;*) args+=("$1" "$2");;esac;shift 2;;
  *) echo "Unexpected argument $1" >&2;exit 1;;
 esac
done
[[ -z $memory || $memory =~ ^[1-9][0-9]*$ ]]||{ echo '--mem requires a positive integer in MB' >&2;exit 1; }
[[ -z $time_limit || $time_limit =~ ^[1-9][0-9]*$ || $time_limit =~ ^[0-9]+:[0-5][0-9]$ ]]||{ echo '--time requires LSF minutes or HH:MM' >&2;exit 1; }
[[ $energies != all ]]||energies=0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250
[[ $energies =~ ^[0-9]+([.][0-9]+)?(,[0-9]+([.][0-9]+)?)*$ ]]||{ echo 'Invalid energy list' >&2;exit 1; }
IFS=, read -ra values <<< "$energies"
declare -A seen=()
for energy in "${values[@]}";do [[ ! -v seen[$energy] ]]||{ echo 'Duplicate energy' >&2;exit 1; };seen[$energy]=1;bash "$HERE/run.sh" --energy "$energy" "${args[@]}" --dry-run;done
logs=$(realpath -m -- "$logs")
make_command(){
 command=(bsub -q "$queue" -n 1 -cwd "$submission_dir" -J "hg_lg_tail_${1}GeV" -o "$logs/${1}GeV-%J.out" -e "$logs/${1}GeV-%J.err")
 [[ -z $memory ]]||command+=(-M "$memory" -R "rusage[mem=$memory]")
 [[ -z $time_limit ]]||command+=(-W "$time_limit")
 command+=(bash "$HERE/run.sh" --energy "$1" "${args[@]}")
}
for energy in "${values[@]}";do make_command "$energy";printf '%q ' "${command[@]}";printf '\n';done
((!dry))||exit 0
command -v bsub >/dev/null
( flock 9; make -C "$HERE" all ) 9>"$HERE/.build.lock"
mkdir -p -- "$logs"
record=$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv
printf 'energy_GeV\tjob_id\n' > "$record"
job_id_pattern='Job <([0-9]+)>'
for energy in "${values[@]}";do
 make_command "$energy"
 reply=$("${command[@]}")||{ echo "Submission failed; earlier jobs: $record" >&2;exit 1; }
 [[ $reply =~ $job_id_pattern ]]||{ echo "Unexpected bsub reply: $reply; earlier jobs: $record" >&2;exit 1; }
 id=${BASH_REMATCH[1]}
 printf '%s\t%s\n' "$energy" "$id" >> "$record"
 printf 'Submitted %s GeV: %s\n' "$energy" "$id"
done
printf 'Submission record: %s\n' "$record"
