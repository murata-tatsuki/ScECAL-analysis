# Edit these input/output paths. Relative paths start at analysis/resolution/.
calibration_after=/path/to/calibration_after
calibration_before=/path/to/calibration_before
result_base=../result/resolution/single
figure_base=../result/resolution/figures/single
sample_type=mc  # Set data or mc to match the calibration inputs.
selection="tail_${tail_event_mode}_front_off"
if [[ -n $front_position_range_mm ]]; then selection="tail_${tail_event_mode}_front${front_position_range_mm}mm"; fi
comparison_type="${sample_type}_vs_${sample_type}"
cog_ranges=(200 20 10 5)

# Empty runner arrays execute locally. For Slurm logs, use analysis/log/resolution/.
# Set single_runner/multi_runner as in hl_intercept_samples.sh.
# Finish all SingleEnergy jobs before starting MultiEnergy.
for cog in "${cog_ranges[@]}"; do
  single_sample "$result_base/$sample_type/my_correction_after/$selection/${cog}mm" "$figure_base/$sample_type/my_correction_after/$selection/${cog}mm" "$cog" 1 '*.root' "$calibration_after"
  single_sample "$result_base/$sample_type/my_correction_before/$selection/${cog}mm" "$figure_base/$sample_type/my_correction_before/$selection/${cog}mm" "$cog" 1 '*.root' "$calibration_before"
done

# Channel canvases are skipped by default; resolution and summaries are saved.
# To also generate channel canvases, uncomment BEFORE adding other multi_options:
# multi_options=()
for cog in "${cog_ranges[@]}"; do
  multi_sample "After" "$result_base/$sample_type/my_correction_after/$selection/${cog}mm/"*.root
  multi_sample "Before" "$result_base/$sample_type/my_correction_before/$selection/${cog}mm/"*.root
done
multi_run "../result/resolution/multi/$comparison_type/my_correction_after_vs_before/$selection/comparison.root" \
  "../result/resolution/figures/multi/$comparison_type/my_correction_after_vs_before/$selection"
