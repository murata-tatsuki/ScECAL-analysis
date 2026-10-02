#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
submission_dir=$(pwd -P)
energies=100 partition=all logs=$HERE/../../log/channel_response/calibration_residual dry=0 time_limit='' memory=''
args=()
while (($#));do
 case $1 in
  -h|--help) cat <<HELP
Usage: bash submit.sh [--energy VALUE|all] [--energies CSV] [run.sh options]
Default: 100 GeV only. One Slurm job per energy, all channels/events by default; run and temperature kept separate.
--energy all includes PS 0.5,1,2,3,4,5 and SPS 10,20,30,40,50,60,70,80,100,120,150,200,250.
--partition NAME (all), --log-dir PATH (default: $HERE/../../log/channel_response/calibration_residual),
--mem SIZE, --time LIMIT, --dry-run.
--dry-run validates inputs and prints sbatch commands; it creates/submits nothing.
Status: bjobs or squeue -u \$USER
HELP
   exit 0;;
  --dry-run) dry=1;shift;;
  --*) (($#>=2))||exit 1;case $1 in --energy|--energies) energies=$2;;--partition) partition=$2;;--log-dir) logs=$2;;--mem) memory=$2;;--time) time_limit=$2;;*) args+=("$1" "$2");;esac;shift 2;;
  *) echo "Unexpected argument $1" >&2;exit 1;;
 esac
done
[[ $energies != all ]]||energies=0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250
[[ $energies =~ ^[0-9]+([.][0-9]+)?(,[0-9]+([.][0-9]+)?)*$ ]]||{ echo 'Invalid energy list' >&2;exit 1; }
IFS=, read -ra values <<< "$energies"
declare -A seen=()
for energy in "${values[@]}";do [[ ! -v seen[$energy] ]]||{ echo 'Duplicate energy' >&2;exit 1; };seen[$energy]=1;bash "$HERE/run.sh" --energy "$energy" "${args[@]}" --dry-run;done
logs=$(realpath -m -- "$logs")
make_command(){
 command=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1 --export=ALL --chdir="$submission_dir" --job-name="calibration_residual_${1}GeV" --output="$logs/${1}GeV-%j.out" --error="$logs/${1}GeV-%j.err")
 [[ -z $memory ]]||command+=(--mem="$memory")
 [[ -z $time_limit ]]||command+=(--time="$time_limit")
 command+=("$HERE/../execute.sl" bash "$HERE/run.sh" --energy "$1" "${args[@]}")
}
for energy in "${values[@]}";do make_command "$energy";printf '%q ' "${command[@]}";printf '\n';done
((!dry))||exit 0
command -v sbatch >/dev/null
( flock 9; make -C "$HERE" all ) 9>"$HERE/.build.lock"
mkdir -p -- "$logs"
record=$logs/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv
printf 'energy_GeV\tjob_id\n' > "$record"
for energy in "${values[@]}";do
 make_command "$energy"
 reply=$("${command[@]}")||{ echo "Submission failed; earlier jobs: $record" >&2;exit 1; }
 id=${reply%%;*};[[ $id =~ ^[0-9]+$ ]]||{ echo "Unexpected sbatch reply: $reply" >&2;exit 1; }
 printf '%s\t%s\n' "$energy" "$reply" >> "$record"
 printf 'Submitted %s GeV: %s\n' "$energy" "$reply"
done
printf 'Submission record: %s\n' "$record"
