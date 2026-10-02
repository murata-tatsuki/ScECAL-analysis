# Compare existing Gaussian-beam and SSA-hitmap SingleEnergy results.
# Usage: bash run_samples.sh --keep-tail-events multi sim_sps_e_ssa_samples.sh
# Relative paths start at analysis/resolution/. MultiEnergy runs locally.
[[ $stage == multi ]] || fail 'This config only compares existing SingleEnergy results; use stage multi'
[[ $tail_event_mode == keep ]] || fail 'These legacy inputs retain HL-tail events; specify --keep-tail-events'

oval_base=../result/resolution/default/simulation/threshold
ssa_base=../result/resolution/custom/hl_intercept/mc/before
comparison_base=../result/resolution/default/comparison/sim_sps_e-_ssa
comparison_figures=../result/resolution/figures/default/comparison/sim_sps_e-_ssa

# CoG fixes the color/marker; oval is solid/filled, SSA hitmap dotted/open.
for cog in 200 20 10 5; do
  multi_sample 'oval' "$oval_base/${cog}mm/"*.root
  multi_sample 'SSA hitmap' "$ssa_base/${cog}mm/"*.root
done

if (( ! dry_run )); then
  mkdir -p -- "$comparison_base"
  {
    date --iso-8601=seconds
    printf '%s\n' 'Command: bash run_samples.sh --keep-tail-events multi sim_sps_e_ssa_samples.sh'
    printf 'Working directory: %s\n' "$PWD"
    printf '%s\n' 'HL-tail events: keep (legacy SingleEnergy inputs); no new event selection'
    printf '%s\n' 'Samples: oval / SSA hitmap; CoG: 200 (nocut), 20, 10, 5 mm'
    printf 'ROOT: %s/all.root\nFigures: %s\n' "$comparison_base" "$comparison_figures"
    printf 'Input: %s\n' "${comparison_files[@]}"
  } > "$comparison_base/run.log"
fi

multi_run "$comparison_base/all.root" "$comparison_figures"
