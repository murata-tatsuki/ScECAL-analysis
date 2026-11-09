#!/usr/bin/env bash
# Select calibrated samples and comparison legends in a sourced shell config.
set -euo pipefail
shopt -s nullglob

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
dry_run=0
tail_event_mode=exclude
tail_results_dir="$script_dir/../result/channel_response/hg_lg_tail/threshold/e_ssa"
usage() {
  echo "Usage: bash run_samples.sh [--dry-run] [--exclude-tail-events|--keep-tail-events] [--tail-results DIR] single|multi config.sh"
  echo "Default: exclude tail-rich events and skip MultiEnergy channel plots. Multi requires SingleEnergy results with the same selection."
}
tail_mode_set=0
while (( $# )); do
  case $1 in
    --dry-run) dry_run=1; shift ;;
    --exclude-tail-events|--keep-tail-events)
      (( tail_mode_set == 0 )) || { echo 'Specify one tail selection option' >&2; exit 1; }
      tail_mode_set=1
      if [[ $1 == --keep-tail-events ]]; then tail_event_mode=keep; else tail_event_mode=exclude; fi
      shift ;;
    --tail-results)
      (( $# >= 2 )) || { echo '--tail-results requires a directory' >&2; exit 1; }
      tail_results_dir=$(realpath -m -- "$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
done
if [[ $# != 2 || ( $1 != single && $1 != multi ) ]]; then
  usage >&2
  exit 1
fi
readonly tail_event_mode
stage=$1
config=$(realpath -- "$2")
cd -- "$script_dir"

# Config may override these arrays, e.g. single_runner=(sbatch ... execute_paralell.sl).
single_runner=()
multi_runner=()
# Skip expensive per-channel comparisons by default. Configs may opt in with multi_options=().
multi_options=(--skip-channel-plots)
multi_output=''
multi_figures=''
labels=()
comparison_files=()
files_per_sample=0
single_jobs=0
multi_jobs=0
declare -A multi_outputs=() multi_figure_paths=()
declare -A single_outputs=()

fail() { echo "run_samples: $*" >&2; exit 1; }
run_command() {
  if (( dry_run )); then
    printf '%q ' "$@"
    printf '\n'
  else
    "$@"
  fi
}

# single_sample OUTPUT_ROOT_DIR FIGURE_DIR COG_MM ONLY_BEST FILE_GLOB CALIB_DIR...
# Each CALIB_DIR contains energy directories such as 0.5GeV/, 10GeV/, 100GeV/.
# Multiple directories (e.g. PS and SPS) are merged by beam energy into one sample.
single_sample() {
  [[ $stage == single ]] || return 0
  (( $# >= 6 )) || fail 'single_sample requires output/figure dirs, CoG, only-best, file glob and calibration dirs'
  local root_dir=$1 figure_dir=$2 cog=$3 only_best=$4 pattern=$5
  shift 5
  [[ -n $root_dir && -n $figure_dir ]] || fail 'Empty SingleEnergy output directory'
  [[ $cog =~ ^[1-9][0-9]*$ ]] || fail "Invalid CoG range: $cog"
  [[ $only_best == 0 || $only_best == 1 ]] || fail 'ONLY_BEST must be 0 or 1'
  local base directory energy file output
  local -a energies=() files=() tail_options=("--${tail_event_mode}-tail-events")
  if [[ $tail_event_mode == exclude ]]; then tail_options+=(--tail-results "$tail_results_dir"); fi
  local -A seen_energies=() seen_files=()
  for base in "$@"; do
    [[ -d $base ]] || fail "Calibration directory not found: $base"
    for directory in "$base"/*GeV; do
      [[ -d $directory ]] || continue
      energy=${directory##*/}
      energy=${energy%GeV}
      # SingleEnergy's existing parser supports integer energies and 0.5 GeV.
      [[ $energy =~ ^[1-9][0-9]*$ || $energy == 0.5 ]] || fail "Unsupported energy directory: $directory"
      if [[ ! ${seen_energies[$energy]+yes} ]]; then
        energies+=("$energy")
        seen_energies[$energy]=1
      fi
    done
  done
  (( ${#energies[@]} )) || fail 'No energy directories found'
  mapfile -t energies < <(printf '%s\n' "${energies[@]}" | sort -g)
  for energy in "${energies[@]}"; do
    files=()
    seen_files=()
    for base in "$@"; do
      # Match basenames without word splitting or evaluating the configured glob.
      for file in "$base/${energy}GeV/"*; do
        [[ -f $file && ${file##*/} == $pattern ]] || continue
        file=$(realpath -- "$file")
        [[ ! ${seen_files[$file]+yes} ]] || fail "Duplicate calibrated input: $file"
        seen_files[$file]=1
        files+=("$file")
      done
    done
    (( ${#files[@]} )) || fail "No files matching $pattern for ${energy}GeV"
    output=$(realpath -m -- "$root_dir/${energy}GeV_${cog}mm.root")
    [[ ! ${single_outputs[$output]+yes} ]] || fail "Repeated SingleEnergy output: $output"
    single_outputs[$output]=1
    if (( ! dry_run )); then mkdir -p -- "$root_dir" "$figure_dir/${energy}GeV/raw"; fi
    run_command "${single_runner[@]}" "$script_dir/SingleEnergyAnalysis" \
      "$output" 1 "${energy}GeV" "${#files[@]}" "${files[@]}" "$cog" "$only_best" "$figure_dir" "${tail_options[@]}"
    single_jobs=$((single_jobs + 1))
  done
}

# multi_sample LEGEND ROOT_FILE... (shell globs or explicit lists are accepted).
multi_sample() {
  [[ $stage == multi ]] || return 0
  (( $# >= 2 )) || fail 'multi_sample requires a legend and at least one ROOT file'
  local label=$1 file
  shift
  [[ -n $label && $label != --* ]] || fail 'Legend must be non-empty and must not start with --'
  if (( ${#labels[@]} == 0 )); then files_per_sample=$#; fi
  (( $# == files_per_sample )) || fail "Different number of energy files for: $label"
  for file in "$@"; do
    [[ -f $file ]] || fail "SingleEnergy ROOT file not found: $file"
  done
  labels+=("$label")
  comparison_files+=("$@")
}

# Finish the current group and reset its samples for the next comparison.
# multi_run OUTPUT_ROOT FIGURE_DIR
multi_run() {
  [[ $stage == multi ]] || return 0
  (( $# == 2 )) || fail 'multi_run requires output ROOT and figure directory'
  local output=$1 figures=$2 file canonical_output canonical_figures
  (( ${#labels[@]} >= 2 )) || fail 'Select at least two multi_sample entries per comparison'
  [[ -n $output && -n $figures ]] || fail 'Set comparison output and figure directory'
  canonical_output=$(realpath -m -- "$output")
  canonical_figures=$(realpath -m -- "$figures")
  [[ ! ${multi_outputs[$canonical_output]+yes} ]] || fail "Repeated comparison output: $output"
  [[ ! ${multi_figure_paths[$canonical_figures]+yes} ]] || fail "Repeated comparison figure directory: $figures"
  for file in "${comparison_files[@]}"; do
    [[ $canonical_output != "$(realpath -- "$file")" ]] || fail 'Comparison output would overwrite an input'
  done
  multi_outputs[$canonical_output]=1
  multi_figure_paths[$canonical_figures]=1
  if (( ! dry_run )); then mkdir -p -- "$(dirname -- "$output")" "$figures"; fi
  local -a selection_options=("--${tail_event_mode}-tail-events")
  local option
  for option in "${multi_options[@]}"; do
    if [[ $option == --compare-tail-selections ]]; then
      (( tail_mode_set == 0 )) || fail '--compare-tail-selections cannot be combined with --exclude-tail-events/--keep-tail-events'
      selection_options=()
    fi
  done
  run_command "${multi_runner[@]}" "$script_dir/MultiEnergyAnalysis" \
    "$output" "${#labels[@]}" "$files_per_sample" "${comparison_files[@]}" "$figures" \
    --labels "${labels[@]}" "${multi_options[@]}" "${selection_options[@]}"
  multi_jobs=$((multi_jobs + 1))
  labels=()
  comparison_files=()
  files_per_sample=0
}

# This is a user-authored shell configuration, with paths relative to resolution/.
source "$config"

if [[ $stage == single ]]; then
  (( single_jobs )) || fail 'No single_sample entries selected'
else
  # Backward compatibility: configs with one implicit comparison still work.
  if (( ${#labels[@]} )); then
    [[ -n $multi_output && -n $multi_figures ]] || fail 'Call multi_run for the pending samples, or set multi_output and multi_figures'
    multi_run "$multi_output" "$multi_figures"
  fi
  (( multi_jobs )) || fail 'No multi_sample entries selected'
fi
