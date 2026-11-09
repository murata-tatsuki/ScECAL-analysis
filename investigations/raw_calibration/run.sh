#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
energy=${1:-100}
[[ $energy =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo 'Specify one beam energy, not all' >&2; exit 1; }
out=${2:-$PROJECT/investigations/result/raw_calibration/${energy}GeV}
params=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
[[ ! -e $out/channels.tsv ]] || { echo "Results already exist: $out" >&2; exit 1; }
make -C "$HERE"
mkdir -p "$out/inputs"
sha256sum "$HERE/"*.cc "$HERE/run.sh" "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$params/pedestal2023_SPS.root" "$params/all_hl_electron2023_hlratio_v3.root" "$params/all_auto_muon_v4_trackfit.root" "$PROJECT/analysis/result/threshold/threshold.root" > "$out/inputs/checksums.txt"
"$HERE/InspectRawCalibration" --input "$PROJECT/analysis/result/impact_studies/hl_intercept/threshold/production_v1/${energy}GeV" --out "$out" --hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root" --threshold "$PROJECT/analysis/result/threshold/threshold.root" --comparison-source "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" 2>&1 | tee "$out/run.log"

bash "$HERE/tails.sh" "$energy" "$out"
"$HERE/SummarizeRaw" "$out"
