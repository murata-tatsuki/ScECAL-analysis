# Select each hit by its calibrated strip-center coordinates [mm].
# Usage: bash run_samples.sh single strip_samples.example.sh
# After all SingleEnergy jobs finish: bash run_samples.sh multi strip_samples.example.sh
# Edit the half-widths here, or choose one with run_samples.sh --strip-range MM.
strip_ranges_mm=("${strip_range_mm:-22.5}" "68.2")
if [[ -n $strip_range_mm ]]; then strip_ranges_mm=("$strip_range_mm"); fi
result_base=../result/resolution/single/data/beam2023_reprocessed
figure_base=../result/resolution/figures/single/data/beam2023_reprocessed
data_base=/gpfs/group/ilc/users/murata/data_beamtest/ECAL_data/analysed/2023
selection="tail_${tail_event_mode}_front_off"
if [[ -n $front_position_range_mm ]]; then selection="tail_${tail_event_mode}_front${front_position_range_mm}mm"; fi
comparison_name=''

# Channel canvases can be enabled with single_options=().
single_options=(--skip-channel-plots)
# Empty runner arrays execute locally; use the existing Slurm wrapper if needed.
# single_runner=(sbatch execute_paralell.sl)

for strip_range_mm in "${strip_ranges_mm[@]}"; do
  condition="${selection}_strip${strip_range_mm}mm"
  comparison_name+="${comparison_name:+_vs_}strip${strip_range_mm}mm"
  single_sample "$result_base/$condition/200mm" \
    "$figure_base/$condition/200mm" 200 1 'ECAL*.root' \
    "$data_base/ps/calib/e-" "$data_base/sps/calib/e-"
  multi_sample "Data, strip |X|,|Y| < ${strip_range_mm} mm" \
    "$result_base/$condition/200mm/"*.root
done

# One range produces only SingleEnergy outputs; compare at least two ranges for multi.
if [[ $stage == multi ]]; then
  (( ${#strip_ranges_mm[@]} >= 2 )) || fail 'Select at least two strip_ranges_mm for MultiEnergy'
  multi_run "../result/resolution/multi/data_vs_data/$comparison_name/$selection/comparison.root" \
    "../result/resolution/figures/multi/data_vs_data/$comparison_name/$selection"
fi
