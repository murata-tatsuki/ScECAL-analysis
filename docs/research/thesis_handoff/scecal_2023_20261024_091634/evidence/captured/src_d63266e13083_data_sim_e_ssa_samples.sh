# Compare existing data and SSA-hitmap MC before-intercept SingleEnergy results.
# Usage: bash run_samples.sh --keep-tail-events multi data_sim_e_ssa_samples.sh
# Relative paths start at analysis/resolution/. MultiEnergy runs locally.
[[ $stage == multi ]] || fail 'This config only compares existing SingleEnergy results; use stage multi'
[[ $tail_event_mode == keep ]] || fail 'These legacy inputs retain HL-tail events; specify --keep-tail-events'

data_base=../result/resolution/default/data
ssa_base=../result/resolution/custom/hl_intercept/mc/before
comparison_base=../result/resolution/default/comparison/data_sim_e-_ssa
comparison_figures=../result/resolution/figures/default/comparison/data_sim_e-_ssa

# Use the same 18 energies in all eight series. The saved data fit at
# 0.5 GeV / 20 mm has a negative mean, so 0.5 GeV is excluded throughout.
energies=(1 2 3 4 5 10 20 30 40 50 60 70 80 100 120 150 200 250)
# CoG fixes the color/marker; data is solid/filled, sim dotted/open.
for cog in 200 20 10 5; do
  data_files=()
  sim_files=()
  for energy in "${energies[@]}"; do
    data_files+=("$data_base/${cog}mm/${energy}GeV_${cog}mm.root")
    sim_files+=("$ssa_base/${cog}mm/${energy}GeV_${cog}mm.root")
  done
  multi_sample 'data' "${data_files[@]}"
  multi_sample 'sim' "${sim_files[@]}"
done

if (( ! dry_run )); then
  mkdir -p -- "$comparison_base"
  {
    date --iso-8601=seconds
    printf '%s\n' 'Command: bash run_samples.sh --keep-tail-events multi data_sim_e_ssa_samples.sh'
    printf 'Working directory: %s\n' "$PWD"
    printf '%s\n' 'HL-tail events: keep (legacy SingleEnergy inputs); no new event selection'
    printf '%s\n' 'Samples: data / sim; CoG: 200 (nocut), 20, 10, 5 mm'
    printf 'ROOT: %s/all.root\nFigures: %s\n' "$comparison_base" "$comparison_figures"
    printf '%s\n' 'Excluded: 0.5 GeV in all series (invalid saved data fit at CoG 20 mm; user approved)'
    printf 'Input: %s\n' "${comparison_files[@]}"
  } > "$comparison_base/run.log"
fi

multi_run "$comparison_base/all.root" "$comparison_figures"
