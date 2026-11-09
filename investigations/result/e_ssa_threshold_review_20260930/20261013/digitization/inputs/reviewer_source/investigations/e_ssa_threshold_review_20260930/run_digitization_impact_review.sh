#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
[[ $# -le 1 && $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||{ echo 'Usage: bash run_digitization_impact_review.sh [LABEL]' >&2;exit 1; }
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
OUT=$BASE/digitization
MANIFEST=$OUT/inputs/manifest.tsv
[[ -f $MANIFEST ]]||{ echo "Missing frozen digitization manifest: $MANIFEST" >&2;exit 1; }
exec 8>"$OUT/inputs/comparison.lock"
flock -n 8||{ echo 'Another process is writing the paired comparison' >&2;exit 1; }
sha256sum -c "$OUT/inputs/source_checksums.txt" > "$OUT/inputs/frozen_production_source_verification.log"
INPUTS=("$MANIFEST")
while IFS=$'\t' read -r task energy file_index old_raw old_cal truth new_raw new_cal;do
 [[ $task == task ]]&&continue
 for path in "$old_raw" "$old_cal" "$truth" "$new_raw" "$new_cal";do [[ -f $path ]]||{ echo "Missing pilot input: $path" >&2;exit 1; };INPUTS+=("$path");done
 dest=${new_raw%/decode/*}
 [[ -f $dest/${new_raw##*/}.COMPLETE ]]||{ echo "Production pilot generation incomplete: $new_raw" >&2;exit 1; }
done < "$MANIFEST"
for energy in 5 20 100;do INPUTS+=("$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv" "$BASE/${energy}GeV/baseline/reference.tsv");done
PARAM_BASE=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
PARAM_FILES=("$PARAM_BASE/all_hl_electron2023_hlratio_v3.root" "$PARAM_BASE/all_auto_muon_v4_trackfit.root" "$PARAM_BASE/pedestal2023_SPS.root" "$PARAM_BASE/20230504_1425_Calib_LED_all.root" "$PROJECT/analysis/result/threshold/threshold.root")
INPUTS+=("${PARAM_FILES[@]}")
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_stats_before.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_sha256_before.txt"
SNAP=$OUT/inputs/reviewer_source
mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan" "$OUT/inputs/reviewer_parameters"
cp -- "$HERE/DigitizationImpactReview.cc" "$HERE/run_digitization_impact_review.sh" "$HERE/ResolutionImpact.cc" "$HERE/CurrentParameters.hh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
cp -- "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
cp -- "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
cp -- "${PARAM_FILES[@]}" "$OUT/inputs/reviewer_parameters/"
sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* "$OUT/inputs/reviewer_parameters/"* > "$OUT/inputs/reviewer_source_parameter_checksums.txt"
printf 'bash %q %q\n' "$HERE/run_digitization_impact_review.sh" "$LABEL" > "$OUT/inputs/comparison_command.sh"
mkdir -p -- "$PROJECT/.scratch"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/digitization_impact_review.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++17 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/DigitizationImpactReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/DigitizationImpactReview" > "$OUT/comparison_build.log" 2>&1
"$TASK_TMP/DigitizationImpactReview" "$BASE" "$OUT" > "$OUT/comparison_run.log" 2>&1
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_stats_after.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_sha256_after.txt"
cmp -s -- "$OUT/inputs/comparison_input_stats_before.tsv" "$OUT/inputs/comparison_input_stats_after.tsv"
cmp -s -- "$OUT/inputs/comparison_input_sha256_before.txt" "$OUT/inputs/comparison_input_sha256_after.txt"
touch "$OUT/COMPARISON_COMPLETE"
printf 'Completed two-file production digitization impact comparison: %s\n' "$OUT"
