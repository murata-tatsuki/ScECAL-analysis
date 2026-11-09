#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930"
input="$base/20261002"
out="$base/resolution_impact"
requested=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold
mkdir -p "$out/inputs"
exec > >(tee "$out/run.log") 2>&1
printf 'Verification date: %s\nCommand: bash %q\nRequested MC: %s\nResolved MC: %s\n' "$(TZ=Asia/Tokyo date --iso-8601=seconds)" "$HERE/resolution_impact.sh" "$requested" "$(realpath "$requested")"
for energy in 5 20 100; do
 dir="$input/${energy}GeV"
 while IFS=$'\t' read -r raw cal; do
  [[ $(realpath "$cal") == "$(realpath "$requested")/"* ]] || { echo 'MC input differs';exit 1; }
 done < "$dir/baseline/inputs/mc.tsv"
 while IFS=$'\t' read -r raw cal; do stat -c '%y\t%s\t%n' "$raw" "$cal";done < <(cat "$dir/baseline/inputs/"{data,mc}.tsv) | diff - "$dir/baseline/inputs/file_stats_after.tsv"
 awk '/Analysis_edit\// {print}' "$dir/baseline/inputs/checksums.txt" | sha256sum -c -
 sha256sum -c "$dir/inputs/tail_checksums.txt"
 sha256sum "$input/sensitivity/comparison/20mm/${energy}GeV/"{events,summary,layer_covariance,covariance_summary}.tsv > "$out/inputs/${energy}GeV_checksums.txt"
 printf '%s GeV: source metadata, parameters, and tail flags verified\n' "$energy"
done
make -C "$HERE" ResolutionImpact
"$HERE/ResolutionImpact" "$input" "$out"
"$HERE/ResolutionImpact" "$input" "$out" --fit-scan
sha256sum "$HERE/ResolutionImpact.cc" "$HERE/resolution_impact.sh" > "$out/inputs/code_checksums.txt"
