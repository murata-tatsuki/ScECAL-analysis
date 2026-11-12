# Compare calibrated SSA MC with and without momentum spread.
# Usage: bash run_samples.sh single|multi momentum_spread_samples.sh
# Existing data / spread MC: MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh multi momentum_spread_samples.sh
# Data SingleEnergy: MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh single momentum_spread_samples.sh
# Data front-cut comparison: MOMENTUM_SPREAD_COMPARE=data_front bash run_samples.sh --front-position-range 20 multi momentum_spread_samples.sh
# Relative paths start at analysis/resolution/.

comparison_target=${MOMENTUM_SPREAD_COMPARE:-mc}
[[ $comparison_target == mc || $comparison_target == data || $comparison_target == data_front ]] || fail 'MOMENTUM_SPREAD_COMPARE must be mc, data or data_front'
if [[ $comparison_target == data || $comparison_target == data_front ]]; then
  [[ $tail_event_mode == exclude ]] || fail 'Data mode requires --exclude-tail-events'
fi
if [[ $comparison_target == data_front ]]; then
  [[ $stage == multi && -n $front_position_range_mm ]] || fail 'data_front requires multi and --front-position-range MM'
fi

mc_base=/gpfs/group/ilc/users/murata/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-
calibration_old="$mc_base/e_ssa/threshold"
calibration_spread="$mc_base/momentum_spread/default"
tail_results_old=../result/channel_response/hg_lg_tail/threshold/e_ssa_kekcc
tail_results_spread=../result/channel_response/hg_lg_tail/momentum_spread/default
result_base=../result/resolution
figure_base="$result_base/figures"
data_calibration=/gpfs/group/ilc/users/murata/data_beamtest/ECAL_data/analysed/2023
front_suffix=''
if [[ -n $front_position_range_mm ]]; then front_suffix="_front${front_position_range_mm}mm"; fi
# Keep the existing job-log condition names; result directories describe the selection explicitly.
condition="$tail_event_mode$front_suffix"
selection="tail_${tail_event_mode}_front_off"
if [[ -n $front_suffix ]]; then selection="tail_${tail_event_mode}_front${front_position_range_mm}mm"; fi
data_single_base="$result_base/single/data/beam2023_reprocessed"
data_single="$data_single_base/$selection"
data_figures="$figure_base/single/data/beam2023_reprocessed/$selection"
old_single="$result_base/single/mc/e_ssa_threshold/$selection"
old_figures="$figure_base/single/mc/e_ssa_threshold/$selection"
spread_single="$result_base/single/mc/spread_ps0p3_sps2/$selection"
spread_figures="$figure_base/single/mc/spread_ps0p3_sps2/$selection"
# The original data inputs remain an unknown-provenance sample unless explicitly replaced.
data_results=${MOMENTUM_SPREAD_DATA_RESULTS:-$result_base/unknown/data}
if [[ -n $front_suffix ]]; then data_results=${MOMENTUM_SPREAD_DATA_RESULTS:-$data_single}; fi
if [[ -n $front_suffix && $comparison_target == data && $stage == multi ]]; then
  [[ $(realpath -m -- "$data_results") == "$(realpath -m -- "$data_single")" ]] || fail 'Front-position comparison requires the matching data SingleEnergy directory; omit MOMENTUM_SPREAD_DATA_RESULTS'
fi
comparison_dir="$result_base/multi/mc_vs_mc/e_ssa_threshold_vs_spread_ps0p3_sps2/$selection"
comparison_root="$comparison_dir/comparison.root"
comparison_figures="$figure_base/multi/mc_vs_mc/e_ssa_threshold_vs_spread_ps0p3_sps2/$selection"
comparison_record="$comparison_dir/inputs"
if [[ $comparison_target == data_front ]]; then
  comparison_dir="$result_base/multi/data_vs_data/front_off_vs_front${front_position_range_mm}/tail_$tail_event_mode"
  comparison_root="$comparison_dir/comparison.root"
  comparison_figures="$figure_base/multi/data_vs_data/front_off_vs_front${front_position_range_mm}/tail_$tail_event_mode"
  comparison_record="$comparison_dir/inputs"
elif [[ $comparison_target == data ]]; then
  if [[ $(realpath -m -- "$data_results") == "$(realpath -m -- "$data_single")" ]]; then
    comparison_dir="$result_base/multi/data_vs_mc/beam2023_reprocessed_vs_spread_ps0p3_sps2/$selection"
    comparison_root="$comparison_dir/comparison.root"
    comparison_figures="$figure_base/multi/data_vs_mc/beam2023_reprocessed_vs_spread_ps0p3_sps2/$selection"
    comparison_record="$comparison_dir/inputs"
  else
    comparison_dir="$result_base/unknown/custom/momentum_spread/comparison"
    comparison_root="$comparison_dir/data_$condition.root"
    comparison_figures="$result_base/unknown/figures/custom/momentum_spread/comparison/data_$condition"
    comparison_record="$comparison_dir/inputs/$condition"
  fi
fi
cog_ranges=(200 20 10 5)
single_options=(--skip-channel-plots)

# Check the selected all-energy tail results before submitting any SingleEnergy jobs.
# The executable also verifies the exact input paths and event identities.
if [[ $stage == single && $tail_event_mode == exclude ]] && (( ! dry_run )); then
  calibration_dirs=("$calibration_old" "$calibration_spread")
  tail_bases=("$tail_results_old" "$tail_results_spread")
  if [[ $comparison_target == data ]]; then
    calibration_dirs=("$data_calibration/ps/calib/e-" "$data_calibration/sps/calib/e-")
    tail_bases=("$tail_results_old")
  fi
  for calibration in "${calibration_dirs[@]}"; do
    for directory in "$calibration"/*GeV; do
      [[ -d $directory ]] || continue
      energy=${directory##*/}
      for tail_base in "${tail_bases[@]}"; do
        for input in tail_study.root inputs/data.tsv inputs/mc.tsv; do
          [[ -f "$tail_base/$energy/$input" ]] || fail "Missing tail result: $tail_base/$energy/$input; run hg_lg_tail on all events for the selected samples first"
        done
      done
    done
  done
fi

# LSF inherits the submitting ROOT/GCC environment. Finish single jobs before multi.
log_dir=$(realpath -m -- ../log/resolution/momentum_spread)
record_dir="$result_base/single/mc/inputs/momentum_spread/$selection"
if [[ $stage == multi ]]; then record_dir=$comparison_record; fi
if [[ $comparison_target == data ]]; then
  log_dir="$log_dir/data"
  if [[ $stage == single ]]; then record_dir="$data_single/inputs"; fi
elif [[ $comparison_target == data_front ]]; then
  log_dir="$log_dir/data_front"
fi
if (( ! dry_run )); then
  mkdir -p -- "$log_dir/$condition/$stage" "$record_dir"
  {
    date --iso-8601=seconds
    printf 'Command: MOMENTUM_SPREAD_COMPARE=%s ' "$comparison_target"
    if [[ $comparison_target == data && $stage == multi ]]; then
      printf 'MOMENTUM_SPREAD_DATA_RESULTS=%q ' "$data_results"
    fi
    printf 'bash run_samples.sh --%s-tail-events ' "$tail_event_mode"
    if [[ -n $front_position_range_mm ]]; then
      printf '%s %q ' '--front-position-range' "$front_position_range_mm"
    fi
    printf '%s momentum_spread_samples.sh\n' "$stage"
    printf 'Working directory: %s\n' "$PWD"
    if [[ $comparison_target == data && $stage == single ]]; then
      printf 'Data inputs: %s\nData inputs: %s\n' "$data_calibration/ps/calib/e-" "$data_calibration/sps/calib/e-"
      printf 'Data tail results: %s\n' "$tail_results_old"
      printf 'ROOT base: %s\nFigure base: %s\n' "$data_single" "$data_figures"
    elif [[ $comparison_target == data_front ]]; then
      printf 'Data without front cut: %s\nData with front cut: %s\n' "$data_single_base/tail_${tail_event_mode}_front_off" "$data_single"
      printf '%s\n' 'Selection provenance: reuse the selections recorded in the existing data SingleEnergy ROOT files'
      printf 'Comparison ROOT: %s\nComparison figures: %s\n' "$comparison_root" "$comparison_figures"
    elif [[ $comparison_target == data ]]; then
      printf 'Data SingleEnergy: %s\nMC SingleEnergy: %s\n' "$data_results" "$spread_single"
      printf '%s\n' 'Selection provenance: reuse the selections recorded in the existing data and MC SingleEnergy ROOT files'
      printf 'Comparison ROOT: %s\nComparison figures: %s\n' "$comparison_root" "$comparison_figures"
    else
      printf 'Old MC: %s\nMomentum-spread MC: %s\n' "$calibration_old" "$calibration_spread"
      printf 'Old MC tail results: %s\nMomentum-spread MC tail results: %s\n' "$tail_results_old" "$tail_results_spread"
      printf 'Old MC SingleEnergy: %s\nMomentum-spread MC SingleEnergy: %s\n' "$old_single" "$spread_single"
      printf 'Comparison ROOT: %s\nComparison figures: %s\n' "$comparison_root" "$comparison_figures"
    fi
    printf 'CoG [mm]: %s\n' "${cog_ranges[*]}"
    printf 'HL-tail events: %s; bad-channel exclusion: 1; channel plots: skipped\n' "$tail_event_mode"
    printf 'Front position half-width [mm]: %s; centre (0,0); CoG selection retained\n' "${front_position_range_mm:-disabled}"
    printf 'LSF queue: %s\n' "${LSF_QUEUE:-s}"
  } > "$record_dir/$stage.run.log"
fi
single_runner=(bsub -q "${LSF_QUEUE:-s}" -n 1 -cwd "$script_dir" -J resolution_momentum_single
  -o "$log_dir/$condition/single/%J.out" -e "$log_dir/$condition/single/%J.err")
multi_runner=(bsub -q "${LSF_QUEUE:-s}" -n 1 -cwd "$script_dir" -J resolution_momentum_multi
  -o "$log_dir/$condition/multi/%J.out" -e "$log_dir/$condition/multi/%J.err")

if [[ $comparison_target == data_front ]]; then
  for cog in "${cog_ranges[@]}"; do
    multi_sample 'Data without front cut' "$data_single_base/tail_${tail_event_mode}_front_off/${cog}mm/"*.root
    multi_sample 'Data with front cut' "$data_single/${cog}mm/"*.root
  done
  multi_run "$comparison_root" "$comparison_figures"
  return 0
fi

if [[ $comparison_target == data ]]; then
  for cog in "${cog_ranges[@]}"; do
    tail_results_dir=$tail_results_old
    single_sample "$data_single/${cog}mm" "$data_figures/${cog}mm" \
      "$cog" 1 'ECAL*.root' "$data_calibration/ps/calib/e-" "$data_calibration/sps/calib/e-"
    multi_sample 'Data' "$data_results/${cog}mm/"*.root
    multi_sample 'MC with momentum spread' "$spread_single/${cog}mm/"*.root
  done
  multi_run "$comparison_root" "$comparison_figures"
  return 0
fi

for cog in "${cog_ranges[@]}"; do
  tail_results_dir=$tail_results_old
  single_sample "$old_single/${cog}mm" "$old_figures/${cog}mm" \
    "$cog" 1 '*.root' "$calibration_old"
  tail_results_dir=$tail_results_spread
  single_sample "$spread_single/${cog}mm" "$spread_figures/${cog}mm" \
    "$cog" 1 '*.root' "$calibration_spread"

  multi_sample 'MC no momentum spread' "$old_single/${cog}mm/"*.root
  multi_sample 'MC with momentum spread' "$spread_single/${cog}mm/"*.root
done
multi_run "$comparison_root" "$comparison_figures"
