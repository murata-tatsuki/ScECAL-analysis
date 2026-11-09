# Edit these input/output paths. Relative paths start at analysis/resolution/.
calibration_after=/path/to/calibration_after
calibration_before=/path/to/calibration_before
result_base=../result/resolution/custom/my_correction
figure_base=../result/resolution/figures/custom/my_correction
sample_type=mc  # Set data or mc to match the calibration inputs.
cog_ranges=(200 20 10 5)

# Empty runner arrays execute locally. For Slurm logs, use analysis/log/resolution/.
# Set single_runner/multi_runner as in hl_intercept_samples.sh.
# Finish all SingleEnergy jobs before starting MultiEnergy.
for cog in "${cog_ranges[@]}"; do
  single_sample "$result_base/$sample_type/after/${cog}mm" "$figure_base/$sample_type/after/${cog}mm" "$cog" 1 '*.root' "$calibration_after"
  single_sample "$result_base/$sample_type/before/${cog}mm" "$figure_base/$sample_type/before/${cog}mm" "$cog" 1 '*.root' "$calibration_before"
done

# Channel canvases are skipped by default; resolution and summaries are saved.
# To also generate channel canvases, uncomment BEFORE adding other multi_options:
# multi_options=()
for cog in "${cog_ranges[@]}"; do
  multi_sample "After" "$result_base/$sample_type/after/${cog}mm/"*.root
  multi_sample "Before" "$result_base/$sample_type/before/${cog}mm/"*.root
done
multi_run "$result_base/comparison/all.root" "$figure_base/comparison/all"
