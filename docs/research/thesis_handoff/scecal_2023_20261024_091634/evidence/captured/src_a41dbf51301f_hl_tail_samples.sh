slurm_log_dir=$(realpath -m -- "$(dirname -- "${BASH_SOURCE[0]}")/../log/resolution")

# Data / current SSA MC, using the same tail-rich flags as hg_lg_tail.
# Run with: bash run_samples.sh [--keep-tail-events] single|multi hl_tail_samples.sh
# Paths are relative to analysis/resolution/. Both modes live under custom/hl_tail.
data_base=/home/murata_t/data_beamtest/ECAL_data/analysed/2023
mc_calibration=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold
result_base=../result/resolution/custom/hl_tail
figure_base=../result/resolution/figures/custom/hl_tail
cog_ranges=(200 20 10 5)

if (( ! dry_run )); then mkdir -p "${slurm_log_dir}/hl_tail/$tail_event_mode/$stage"; fi
single_runner=(sbatch -o "${slurm_log_dir}/hl_tail/$tail_event_mode/single/%j.out" --error="${slurm_log_dir}/hl_tail/$tail_event_mode/single/%j.err" execute_paralell.sl)
multi_runner=(sbatch -o "${slurm_log_dir}/hl_tail/$tail_event_mode/multi/%j.out" --error="${slurm_log_dir}/hl_tail/$tail_event_mode/multi/%j.err" execute_paralell.sl)
# Finish all single jobs before submitting multi. Empty runner arrays run locally.

for cog in "${cog_ranges[@]}"; do
  single_sample "$result_base/data/$tail_event_mode/${cog}mm" "$figure_base/data/$tail_event_mode/${cog}mm" "$cog" 1 'ECAL*.root' \
    "$data_base/ps/calib/e-" "$data_base/sps/calib/e-"
  single_sample "$result_base/mc/$tail_event_mode/${cog}mm" "$figure_base/mc/$tail_event_mode/${cog}mm" "$cog" 1 '*.root' "$mc_calibration"
done

# MultiEnergy: choose data (default) or MC, with four CoGs and both tail selections.
# The selected target's exclude and keep SingleEnergy jobs must have finished first.
multi_target=${HL_TAIL_MULTI_TARGET:-data}
[[ $multi_target == data || $multi_target == mc ]] || fail "HL_TAIL_MULTI_TARGET must be data or mc: $multi_target"
multi_options+=(--compare-tail-selections)
tail_comparison_root=../result/resolution/custom/hl_tail
tail_comparison_figures=../result/resolution/figures/custom/hl_tail
read -r -a excluded_energies <<< "${HL_TAIL_MULTI_EXCLUDE_ENERGIES:-}"
for cog in "${cog_ranges[@]}"; do
  for mode in exclude keep; do
    input_dir="$tail_comparison_root/$multi_target/$mode/${cog}mm"
    input_files=()
    for root_file in "$input_dir/"*.root; do
      energy_name=${root_file##*/}
      energy_name=${energy_name%%GeV_*}
      skip_energy=0
      for excluded_energy in "${excluded_energies[@]}"; do
        if [[ $energy_name == "$excluded_energy" ]]; then skip_energy=1; break; fi
      done
      if (( ! skip_energy )); then input_files+=("$root_file"); fi
    done
    multi_sample "${multi_target^^} tail $mode" "${input_files[@]}"
  done
done

multi_run "$tail_comparison_root/comparison/${multi_target}_all_cog_tail.root" "$tail_comparison_figures/comparison/${multi_target}_all_cog_tail"
