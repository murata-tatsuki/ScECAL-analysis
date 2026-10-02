#!/bin/bash
# Fill one sample pair per energy, then render a throttled energy/layer array.
set -euo pipefail
HERE=$(dirname -- "$(realpath -e -- "${BASH_SOURCE[0]}")")
submission_dir=$(pwd -P)
ps_energies=0.5,1,2,3,4,5
sps_energies=10,20,30,40,50,60,70,80,100,120,150,200,250
energies_csv=$ps_energies,$sps_energies
beam=auto partition=all log_dir=$HERE/../../log/channel_response/adc_energy
layers=0-29 channels=all plot_concurrency=32
time_limit='' memory='' energy_option=''
dry_run=0 no_plots=0 plots_only=0
run_args=()
fail(){ printf 'ERROR: %s\n' "$*" >&2;exit 1; }
usage(){ cat <<HELP
Usage: bash submit.sh [submission options] [run.sh options]
  --energy VALUE|all    One energy or all (default: all PS+SPS energies)
  --energies LIST       Comma-separated energies
  --beam auto|ps|sps    Select data beam and default energy list
  --plot-concurrency N Maximum simultaneous layer plot tasks (default: 32)
  --no-plots           Submit histogram aggregation only
  --plots-only         Plot completed sample ROOTs without re-reading events
  --save-canvases       Also save per-layer comparison ROOTs (default: PNG only)
  --scratch-root PATH  Node-local plotting workspace/cache (default: /tmp)
  --partition NAME     Default: all
  --time LIMIT         Optional time limit for each task
  --mem SIZE           Optional memory per task
  --log-dir PATH       Logs and submission/task records (default: $HERE/../../log/channel_response/adc_energy)
  --dry-run            Validate inputs and print commands; create/submit nothing

Other run.sh options, including --label, --channels, --layers, --max-events,
and input/output paths, are forwarded. Defaults analyse all events/channels.
Each energy has one fill job. After ALL fill jobs succeed, one array renders
one energy/layer per task with a global concurrency cap for this submission.
--plots-only skips fill jobs and uses completed data/MC histogram ROOTs.
Examples:
  bash submit.sh --energy 100 --label fast_v1
  bash submit.sh --energy all --plot-concurrency 32 --label fast_v1
  bash submit.sh --plots-only --energy 100
Status: bjobs or squeue -u \$USER
HELP
}
while (($#));do
 case $1 in
  -h|--help) usage;exit 0;;
  --dry-run) dry_run=1;shift;;
  --no-plots) no_plots=1;shift;;
  --plots-only) plots_only=1;shift;;
  --save-canvases) run_args+=("$1");shift;;
  --*)
   (($#>=2))||fail "Missing value for $1"
   case $1 in
    --energy|--energies)
     [[ -z $energy_option ]]||fail 'Specify energies only once'
     energy_option=$1;energies_csv=$2
     [[ $1 != --energy || $2 != *,* ]]||fail 'Use --energies for a list';;
    --beam) beam=$2;;
    --partition) partition=$2;;
    --time) time_limit=$2;;
    --mem) memory=$2;;
    --log-dir) log_dir=$2;;
    --plot-concurrency) plot_concurrency=$2;;
    --layers) layers=$2;run_args+=("$1" "$2");;
    --channels) channels=$2;run_args+=("$1" "$2");;
    *) run_args+=("$1" "$2");;
   esac
   shift 2;;
  *) fail "Unexpected argument $1";;
 esac
done
((!no_plots || !plots_only))||fail '--no-plots and --plots-only cannot be combined'
[[ $plot_concurrency =~ ^[1-9][0-9]*$ ]]&&((plot_concurrency<=128))||fail 'plot-concurrency must be 1..128'
if [[ $energies_csv == all ]];then energies_csv=$ps_energies,$sps_energies;energy_option='';fi
case $beam in
 auto) ;;
 ps|sps)
  if [[ -z $energy_option ]];then
   if [[ $beam == ps ]];then energies_csv=$ps_energies;else energies_csv=$sps_energies;fi
  fi
  run_args+=(--beam "$beam");;
 *) fail 'beam must be auto, ps, or sps';;
esac
[[ $energies_csv =~ ^[A-Za-z0-9_.-]+(,[A-Za-z0-9_.-]+)*$ ]]||fail 'Invalid energy list'
[[ $partition =~ ^[A-Za-z0-9_.-]+$ ]]||fail 'Invalid partition'
IFS=, read -r -a energies <<< "$energies_csv"
declare -A seen=() selected_layers=()
for energy in "${energies[@]}";do
 [[ ! -v seen[$energy] ]]||fail "Duplicate energy $energy";seen[$energy]=1
 printf 'Checking %s GeV...\n' "$energy"
 mode=--no-plots;((!plots_only))||mode=--plots-only
 bash "$HERE/run.sh" --energy "$energy" "${run_args[@]}" "$mode" --dry-run
done
[[ $layers =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]]||fail 'Invalid layers'
IFS=, read -r -a ranges <<< "$layers"
for range in "${ranges[@]}";do
 lo=${range%%-*};hi=${range##*-};lo=$((10#$lo));hi=$((10#$hi))
 ((lo<=hi && hi<32))||fail 'Layers must be within 0..31'
 for ((l=lo;l<=hi;++l));do selected_layers[$l]=1;done
done
if [[ $channels != all ]];then
 declare -A channel_layers=()
 IFS=, read -r -a cells <<< "$channels"
 for cell in "${cells[@]}";do channel_layers[$((10#$cell/100000))]=1;done
 for l in "${!selected_layers[@]}";do [[ -v channel_layers[$l] ]]||unset 'selected_layers[$l]';done
fi
plot_layers=()
for ((l=0;l<32;++l));do [[ ! -v selected_layers[$l] ]]||plot_layers+=("$l");done
((${#plot_layers[@]}))||fail 'No selected layers contain requested channels'
log_dir=$(realpath -m -- "$log_dir")
record_tag=$(date +%Y%m%d_%H%M%S)_$$
tasks=$log_dir/plot_tasks_$record_tag.tsv
task_count=$((${#energies[@]}*${#plot_layers[@]}))
base_command(){
 command=(sbatch --parsable --partition="$partition" --ntasks=1 --cpus-per-task=1 --export=ALL --chdir="$submission_dir")
 [[ -z $time_limit ]]||command+=(--time="$time_limit")
 [[ -z $memory ]]||command+=(--mem="$memory")
}
fill_command(){
 local energy=$1
 base_command
 command+=(--job-name="channel_fill_${energy}GeV" --output="$log_dir/${energy}GeV-fill-%j.out" --error="$log_dir/${energy}GeV-fill-%j.err"
  "$HERE/../execute.sl" bash "$HERE/run.sh" --energy "$energy" "${run_args[@]}" --no-plots)
}
plot_command(){
 local dependency=$1
 base_command
 command+=(--job-name=channel_plot --array="0-$((task_count-1))%$plot_concurrency"
  --output="$log_dir/plot-%A_%a.out" --error="$log_dir/plot-%A_%a.err")
 [[ -z $dependency ]]||command+=(--dependency="afterok:$dependency" --kill-on-invalid-dep=yes)
 command+=("$HERE/../execute.sl" bash "$HERE/plot_array.sh" "$tasks" "${run_args[@]}")
}
print_command(){ printf '%q ' "$@";printf '\n'; }
if ((!plots_only));then
 for energy in "${energies[@]}";do fill_command "$energy";print_command "${command[@]}";done
fi
if ((!no_plots));then
 dependency='';((plots_only))||dependency=FILL_JOB_IDS
 plot_command "$dependency";print_command "${command[@]}"
 printf 'Plot tasks: %s energies x %s layers = %s tasks, concurrency %s\n' "${#energies[@]}" "${#plot_layers[@]}" "$task_count" "$plot_concurrency"
fi
((!dry_run))||exit 0
command -v sbatch >/dev/null||fail 'sbatch is not available'
make -C "$HERE" all
mkdir -p -- "$log_dir"
record=$log_dir/submissions_$record_tag.tsv
printf 'stage\tenergy_GeV\tjob_id\n' > "$record"
submit(){
 local stage=$1 energy=$2 reply
 if ! reply=$("${command[@]}");then fail "Submission failed; earlier jobs are recorded in $record";fi
 job_id=${reply%%;*}
 [[ $job_id =~ ^[0-9]+$ ]]||fail "Unexpected sbatch reply $reply; see $record"
 printf '%s\t%s\t%s\n' "$stage" "$energy" "$reply" >> "$record"
 printf 'Submitted %s %s: %s\n' "$stage" "$energy" "$reply"
}
ids=()
if ((!plots_only));then
 for energy in "${energies[@]}";do fill_command "$energy";submit fill "$energy";ids+=("$job_id");done
fi
if ((!no_plots));then
 for energy in "${energies[@]}";do for l in "${plot_layers[@]}";do printf '%s\t%s\n' "$energy" "$l";done;done > "$tasks"
 dependency=$(IFS=:;printf '%s' "${ids[*]}")
 plot_command "$dependency";submit plot_array all
fi
printf 'Submission record: %s\nCheck jobs with: bjobs\n' "$record"
