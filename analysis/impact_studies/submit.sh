#!/bin/bash
# Submit one Slurm job per beam energy; run.sh defaults to all channels/events.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
submission_dir=$(pwd -P)
source "$HERE/energies.sh"
energies_csv=$all_energies
beam=auto
partition=all
log_dir=$HERE/../log/impact_studies
time_limit=''
memory=''
dry_run=0
energy_option=''
run_args=()
summary_args=()
summarize=0

usage() {
    cat <<EOF
Usage: bash submit.sh [submission options] [run.sh options]
  --energy VALUE|all   Submit one energy, or all energies (also the default)
  --energies LIST      Comma-separated energies, one job each
  --beam NAME          auto (default), ps, or sps; selects default energies too
                       PS: $ps_energies GeV
                       SPS: $sps_energies GeV
                       auto: submit both groups; select data beam by energy
  --partition NAME     Slurm partition (default: all)
  --time LIMIT         Optional Slurm time limit; otherwise partition default
  --mem SIZE           Optional memory per job, e.g. 8G; otherwise cluster default
  --log-dir PATH       Job logs and submission records (default: $HERE/../log/impact_studies)
  --dry-run            Check all inputs and show commands without submitting/building
  --summarize          Submit an energy-scan job after all analysis jobs succeed
  -h, --help           Show this help

Other options are passed to run.sh, including --mc-tag, --label, --plot-channels,
--max-events, --result-dir, --study and --no-plots.
Without selection options, each job processes all connected channels in
layers 0-29 and all events for both data and MC, then writes PNG figures.
Paths relative to the submission directory keep the same meaning in each job.

Examples:
  bash submit.sh --energy 100 --dry-run
  bash submit.sh --energy 100
  bash submit.sh --beam ps
  bash submit.sh --energies 30,50,100,150,250
  bash submit.sh

Status: bjobs (Slurm compatibility command on this cluster), or squeue -u \$USER
EOF
}
fail() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
while (($#)); do
    case $1 in
        -h|--help) usage; exit 0 ;;
        --dry-run) dry_run=1; shift ;;
        --summarize) summarize=1; shift ;;
        --no-plots) run_args+=("$1"); shift ;;
        --*)
            (($# >= 2)) || fail "Missing value for $1"
            case $1 in
                --energy|--energies)
                    [[ -z $energy_option ]] || fail 'Use --energy or --energies only once'
                    energy_option=$1
                    energies_csv=$2
                    [[ $1 != --energy || $2 != *,* ]] || fail 'Use --energies for multiple energies' ;;
                --partition) partition=$2 ;;
                --beam) beam=$2 ;;
                --time) time_limit=$2 ;;
                --mem) memory=$2 ;;
                --log-dir) log_dir=$2 ;;
                --study|--mc-tag|--label|--result-dir)
                    run_args+=("$1" "$2"); summary_args+=("$1" "$2") ;;
                *) run_args+=("$1" "$2") ;;
            esac
            shift 2 ;;
        *) fail "Unexpected argument: $1" ;;
    esac
done
if [[ $energies_csv == all ]];then
    energies_csv=$all_energies
    energy_option=''
fi
case $beam in
    auto) ;;
    ps|sps)
        if [[ -z $energy_option ]]; then
            if [[ $beam == ps ]]; then energies_csv=$ps_energies; else energies_csv=$sps_energies; fi
        fi
        run_args+=(--beam "$beam") ;;
    *) fail 'beam must be auto, ps, or sps' ;;
esac
[[ $energies_csv =~ ^[A-Za-z0-9_.-]+(,[A-Za-z0-9_.-]+)*$ ]] || fail 'Invalid energy list'
[[ $partition =~ ^[A-Za-z0-9_.-]+$ ]] || fail 'Invalid partition'
IFS=, read -r -a energies <<< "$energies_csv"
declare -A seen=()
for energy in "${energies[@]}"; do
    [[ ! -v seen[$energy] ]] || fail "Duplicate energy: $energy"
    seen[$energy]=1
done
log_dir=$(realpath -m -- "$log_dir")

# Validate every requested energy before submitting even the first job.
# This checks input pairs, calibration constants, options and output collisions.
for energy in "${energies[@]}"; do
    printf 'Checking %s GeV...\n' "$energy"
    bash "$HERE/run.sh" --energy "$energy" "${run_args[@]}" --dry-run
done

command_for_energy() {
    local energy=$1
    command=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1
        --export=ALL --chdir="$submission_dir" --job-name="hl_intercept_${energy}GeV"
        --output="$log_dir/${energy}GeV-%j.out" --error="$log_dir/${energy}GeV-%j.err")
    [[ -z $time_limit ]] || command+=(--time="$time_limit")
    [[ -z $memory ]] || command+=(--mem="$memory")
    command+=("$HERE/execute.sl" bash "$HERE/run.sh" --energy "$energy" "${run_args[@]}")
}
for energy in "${energies[@]}"; do
    command_for_energy "$energy"
    printf '%q ' "${command[@]}"
    printf '\n'
done
summary_command() {
    local dependencies=$1
    command=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1
        --export=ALL --chdir="$submission_dir" --job-name=hl_intercept_energy_scan
        --dependency="afterok:$dependencies" --kill-on-invalid-dep=yes
        --output="$log_dir/energy_scan-%j.out" --error="$log_dir/energy_scan-%j.err")
    [[ -z $time_limit ]] || command+=(--time="$time_limit")
    [[ -z $memory ]] || command+=(--mem="$memory")
    command+=("$HERE/execute.sl" bash "$HERE/summarize.sh" "${summary_args[@]}" --energies "$energies_csv")
}
if ((summarize && dry_run));then
    summary_command ENERGY_JOB_IDS
    printf '%q ' "${command[@]}"; printf '\n'
fi
((!dry_run)) || exit 0
command -v sbatch >/dev/null || fail 'sbatch is not available; run this on the cluster login host'

# Compile before jobs start so normal parallel jobs see up-to-date executables.
make -C "$HERE" all
mkdir -p -- "$log_dir"
submission_log=$log_dir/submissions_$(date +%Y%m%d_%H%M%S)_$$.tsv
printf 'energy_GeV\tjob_id\n' > "$submission_log"
printf 'Submission record: %s\n' "$submission_log"
job_ids=()
for energy in "${energies[@]}"; do
    command_for_energy "$energy"
    if job_id=$("${command[@]}"); then
        dependency_id=${job_id%%;*}
        [[ $dependency_id =~ ^[0-9]+$ ]] || fail "Unexpected sbatch job ID: $job_id"
        job_ids+=("$dependency_id")
        printf '%s\t%s\n' "$energy" "$job_id" >> "$submission_log"
        printf 'Submitted %s GeV: job %s\n' "$energy" "$job_id"
    else
        fail "Submission failed for $energy GeV; earlier submissions are listed in $submission_log"
    fi
done
if ((summarize));then
    dependencies=$(IFS=:;printf '%s' "${job_ids[*]}")
    summary_command "$dependencies"
    if job_id=$("${command[@]}");then
        printf 'energy_scan\t%s\n' "$job_id" >> "$submission_log"
        printf 'Submitted energy summary after all analysis jobs succeed: %s\n' "$job_id"
    else
        fail "Summary submission failed; analysis jobs are listed in $submission_log"
    fi
fi
printf 'Check jobs with: bjobs\n'
