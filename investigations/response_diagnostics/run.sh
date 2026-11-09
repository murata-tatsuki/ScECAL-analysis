#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
root=$(cd -- "$HERE/../.." && pwd -P)
base=${1:-$root/analysis/result/impact_studies/hl_intercept/threshold/production_v1}
out=${2:-$root/investigations/result/response_diagnostics/threshold/production_v1}
params=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
args=(--hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root")
make -C "$HERE" all
mkdir -p "$out/inputs"
printf 'input_label_directory=%s\n' "$base" > "$out/inputs/configuration.txt"
sha256sum "$HERE/"*.cc "$HERE/"*.sh "$params/all_auto_muon_v4_trackfit.root" "$params/all_hl_electron2023_hlratio_v3.root" "$params/pedestal2023_SPS.root" >> "$out/inputs/configuration.txt"
cp "$root/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$out/inputs/comparison_source.cc"
"$HERE/SurveyImpact" --base "$base" --out "$out"
"$HERE/AuditMIP" "${args[@]}" --out "$out/mip_constants_audit.tsv"
bash "$HERE/inspect.sh" "$base" "$out"
for energy in 1 5 100 200;do "$HERE/InspectTopology" --input "$base/${energy}GeV" --out "$out/${energy}GeV" --comparison-source "$out/inputs/comparison_source.cc";done
for energy in 1 100;do "$HERE/CheckTruthResponse" --input "$base/${energy}GeV" --audit "$out/mip_constants_audit.tsv" "${args[@]}" --out "$out/${energy}GeV" --events 10000;done
"$HERE/PlotResponse" "$out"
bash "$HERE/summarize.sh" "$out" "$base"
bash "$HERE/verify_results.sh" "$out"
