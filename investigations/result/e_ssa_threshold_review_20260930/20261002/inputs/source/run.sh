#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
energy=${1:-100}
label=${2:-}
case "$energy" in 100|20) beam=sps;; 5) beam=ps;; *) echo 'Representative energies only: 5, 20 or 100' >&2; exit 1;; esac
[[ $label != *[^a-zA-Z0-9_-]* ]] || { echo 'Invalid result label' >&2; exit 1; }
MC=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
DATA=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930${label:+/$label}"
out="$base/${energy}GeV"
[[ ! -e $out/baseline/study.root ]] || { echo "Output already exists: $out" >&2;exit 1; }
mkdir -p "$out/baseline/inputs" "$out/raw" "$out/calibrated" "$out/truth"
params="$MC/Analysis_edit/share"
args=(--hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root")
shopt -s nullglob
files=("$DATA/$beam/decode/e-/${energy}GeV/"*.root)
(( ${#files[@]} > 0 )) || { echo 'No data files' >&2; exit 1; }
for a in "${files[@]}"; do b="$DATA/$beam/calib/e-/${energy}GeV/${a##*/}";[[ -f $b ]];printf '%s\t%s\n' "$a" "$b";done > "$out/baseline/inputs/data.tsv"
files=("$MC/Result_MC/decode/e-/e_ssa/threshold/${energy}GeV/"*.root)
(( ${#files[@]} >= 6 )) || { echo 'Fewer than six MC files' >&2; exit 1; }
for i in 0 1 2 3 4 5;do idx=$((i*(${#files[@]}-1)/5));a=${files[$idx]};b="$MC/Result_MC/calib/e-/e_ssa/threshold/${energy}GeV/${a##*/}";[[ -f $a && -f $b ]];printf '%s\t%s\n' "$a" "$b";done > "$out/baseline/inputs/mc.tsv"
{
 date --iso-8601=seconds
 printf 'Command: bash %q %q %q\n' "$HERE/run.sh" "$energy" "$label"
 printf 'MC available files: %s; sampled: 6; data: 2000 events/file in ten spread blocks\n' "${#files[@]}"
 printf 'MC source: %s/Result_MC/{decode,calib}/e-/e_ssa/threshold/%sGeV\n' "$MC" "$energy"
} > "$out/run.log"
while IFS=$'\t' read -r a b; do stat -c '%y\t%s\t%n' "$a" "$b"; done < <(cat "$out/baseline/inputs/"{data,mc}.tsv) > "$out/baseline/inputs/file_stats_before.tsv"
sha256sum "$HERE/run.sh" "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$MC/Analysis_edit/src/"{MCDigi,Extract}.cxx "$params/"{all_hl_electron2023_hlratio_v3.root,all_auto_muon_v4_trackfit.root,pedestal2023_SPS.root} > "$out/baseline/inputs/checksums.txt"
make -C "$HERE" PrepareReference
"$HERE/PrepareReference" --mc-manifest "$out/baseline/inputs/mc.tsv" "${args[@]}" --energy "$energy" --out "$out/baseline" > "$out/baseline/run.log" 2>&1
echo "$energy GeV: paired MC/calibration baseline complete"
"$PROJECT/investigations/raw_calibration/InspectRawCalibration" --input "$out/baseline" --out "$out/raw" "${args[@]}" --threshold "$PROJECT/analysis/result/threshold/threshold.root" --comparison-source "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" > "$out/raw/run.log" 2>&1
echo "$energy GeV: decode/calibration audit complete"
"$PROJECT/investigations/response_diagnostics/InspectResponse" --input "$out/baseline" --out "$out/calibrated" "${args[@]}" --data-per-file 2000 --mc-per-file 1000 --comparison-source "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" > "$out/calibrated/run.log" 2>&1
echo "$energy GeV: calibrated response complete"
while IFS=$'\t' read -r a b; do stat -c '%y\t%s\t%n' "$a" "$b"; done < <(cat "$out/baseline/inputs/"{data,mc}.tsv) > "$out/baseline/inputs/file_stats_after.tsv"
cmp "$out/baseline/inputs/file_stats_before.tsv" "$out/baseline/inputs/file_stats_after.tsv"
if [[ ${REVIEW_CORE_ONLY:-0} != 1 ]]; then bash "$HERE/extra.sh" "$energy" "$label"; fi
