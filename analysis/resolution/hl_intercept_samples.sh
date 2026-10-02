slurm_log_dir=$(realpath -m -- "$(dirname -- "${BASH_SOURCE[0]}")/../log/resolution")

# Compare SSA simulation with / without the HG–LG intercalibration intercept.
# All relative paths are interpreted from analysis/resolution/.

calibration_a=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold
calibration_b=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa_noHLintercept/threshold
result_base=../result/resolution/custom/hl_intercept
figure_base=../result/resolution/figures/custom/hl_intercept

# Submit one Slurm job per sample / CoG / energy, as in analysis_paralell.sh.
# execute_paralell.sl sets partition=all, ntasks=1 and job-name=analyzer.
if [[ $dry_run == 0 ]]; then
  if [[ $stage == single ]]; then
    mkdir -p "${slurm_log_dir}/singleEnergy"
  elif [[ $stage == multi ]]; then
    mkdir -p "${slurm_log_dir}/compare"
  fi
fi
single_runner=(sbatch -o "${slurm_log_dir}/singleEnergy/test-%A.out" --error="${slurm_log_dir}/singleEnergy/test-%A.err" execute_paralell.sl)
multi_runner=(sbatch -o "${slurm_log_dir}/compare/test-%A.out" --error="${slurm_log_dir}/compare/test-%A.err" execute_paralell.sl)
# Both stages use Slurm. Submission returns before analysis finishes.
# Submit the multi stage after all SingleEnergy jobs have completed successfully.

# CoG ranges used by both SingleEnergy and MultiEnergy.
cog_ranges=(200 20 10 5)

# Each call selects the calibration input(s), ROOT/figure destinations and cuts.
# You can pass PS and SPS calibration directories together after '*.root'.
for cog in "${cog_ranges[@]}"; do
  single_sample "$result_base/mc/after/${cog}mm" "$figure_base/mc/after/${cog}mm" "$cog" 1 '*.root' "$calibration_a"
  single_sample "$result_base/mc/before/${cog}mm" "$figure_base/mc/before/${cog}mm" "$cog" 1 '*.root' "$calibration_b"
done

# Channel canvases are skipped by default; resolution and summaries are saved.
# To also generate channel canvases, uncomment BEFORE adding other multi_options:
# multi_options=()

# One comparison: four CoG selections times after/before = eight series.
# The first condition is solid; the second is dotted. CoG is added to captions.
for comparison_cog in "${cog_ranges[@]}"; do
  for condition in after before; do
    input_dir="$result_base/mc/$condition/${comparison_cog}mm"
    # Read old results only when this selection has no new directory yet.
    if [[ ! -d $input_dir && -d $result_base/$condition/${comparison_cog}mm ]]; then
      input_dir="$result_base/$condition/${comparison_cog}mm"
    fi
    multi_sample "MC $condition" "$input_dir/"*.root
  done
done
multi_run "$result_base/mc/comparison/all.root" "$figure_base/mc/comparison/all"
