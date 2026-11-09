#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
energy=${1:-100};label=${2:-}
case "$energy" in 100|20|5);;*) echo 'Use 5, 20 or 100' >&2;exit 1;;esac
[[ $label != *[^a-zA-Z0-9_-]* ]] || { echo 'Invalid result label' >&2;exit 1; }
params=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
args=(--hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root")
result_base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930${label:+/$label}"
out="$result_base/${energy}GeV"
tail="$PROJECT/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/${energy}GeV"
[[ -f $out/baseline/study.root && -f $tail/tail_study.root ]]
mkdir -p "$out/truth" "$out/inputs"
make -C "$HERE" TruthResponse Sensitivity
"$PROJECT/investigations/response_diagnostics/AuditMIP" "${args[@]}" --out "$out/inputs/mip_audit.tsv"
# Reuse production tail flags only if all selected input files predate that result.
while IFS=$'\t' read -r raw cal; do
 [[ ! $raw -nt $tail/tail_study.root && ! $cal -nt $tail/tail_study.root ]] || { echo "Stale tail flags: $cal" >&2;exit 1; }
done < <(cat "$out/baseline/inputs/"{data,mc}.tsv)
sha256sum "$tail/tail_study.root" "$tail/inputs/"{data,mc}.tsv > "$out/inputs/tail_checksums.txt"
printf '%s\n' "$tail" > "$out/inputs/tail_source.txt"
"$HERE/TruthResponse" --manifest "$out/baseline/inputs/mc.tsv" --audit "$out/inputs/mip_audit.tsv" "${args[@]}" --out "$out/truth"
"$HERE/Sensitivity" --input "$out/baseline" --out "$result_base/sensitivity" --audit "$out/inputs/mip_audit.tsv" "${args[@]}" --mask "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" --tail-results "$tail"
sha256sum "$HERE/"*.cc "$HERE/"*.sh "$PROJECT/investigations/raw_calibration/InspectRawCalibration.cc" "$PROJECT/investigations/response_diagnostics/InspectResponse.cc" "$PROJECT/analysis/resolution/HLTailSelection.hh" "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" > "$out/baseline/inputs/investigation_code_checksums.txt"
while IFS=$'\t' read -r a b; do stat -c '%y\t%s\t%n' "$a" "$b"; done < <(cat "$out/baseline/inputs/"{data,mc}.tsv) > "$out/baseline/inputs/file_stats_after.tsv"
cmp "$out/baseline/inputs/file_stats_before.tsv" "$out/baseline/inputs/file_stats_after.tsv"
