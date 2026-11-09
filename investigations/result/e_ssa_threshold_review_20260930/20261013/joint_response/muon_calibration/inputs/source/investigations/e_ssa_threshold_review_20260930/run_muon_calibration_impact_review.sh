#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
[[ $LABEL =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid label' >&2;exit 1; }
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
OUT=$BASE/joint_response/muon_calibration
OLD=$PROJECT/analysis/result/channel_response/muon_mip_response/default
NEW=$PROJECT/analysis/result/channel_response/muon_mip_response/threshold_MPVwoCT
[[ -e $BASE/digitization/COMPARISON_COMPLETE && -e $OLD/COMPLETE && -e $NEW/COMPLETE && ! -e $OUT/COMPLETE ]] || { echo 'Inputs incomplete or diagnostic already complete' >&2;exit 1; }
mkdir -p -- "$OUT/inputs" "$PROJECT/.scratch"
exec 8>"$OUT/inputs/run.lock"
flock -n 8 || { echo 'Diagnostic running' >&2;exit 1; }
INPUTS=("$OLD/fits.tsv" "$OLD/plateau_stability.tsv" "$NEW/fits.tsv" "$NEW/plateau_stability.tsv" "$BASE/digitization/inputs/manifest.tsv")
SMALL_INPUTS=("${INPUTS[@]}")
for energy in 5 20 100;do
 INPUTS+=("$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv" "$BASE/digitization/comparison/20mm/${energy}GeV/paired_events.tsv" "$BASE/${energy}GeV/baseline/inputs/data.tsv")
 SMALL_INPUTS+=("$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv" "$BASE/digitization/comparison/20mm/${energy}GeV/paired_events.tsv" "$BASE/${energy}GeV/baseline/inputs/data.tsv")
 while IFS=$'\t' read -r raw calib;do INPUTS+=("$calib");done < "$BASE/${energy}GeV/baseline/inputs/data.tsv"
done
while IFS=$'\t' read -r task energy file_index oldraw oldcal truth newraw newcal;do [[ $task == task ]] && continue;INPUTS+=("$newcal");done < "$BASE/digitization/inputs/manifest.tsv"
INPUTS+=("$BASE/digitization/inputs/reviewer_parameters/"*.root)
SMALL_INPUTS+=("$BASE/digitization/inputs/reviewer_parameters/"*.root)
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/input_stats_before.tsv"
sha256sum -- "${SMALL_INPUTS[@]}" > "$OUT/inputs/input_sha256_before.txt"
SNAP=$OUT/inputs/source
mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies"
for f in MuonCalibrationImpactReview.cc DigitizationImpactReview.cc ResolutionImpact.cc CurrentParameters.hh run_muon_calibration_impact_review.sh;do cp -- "$HERE/$f" "$SNAP/investigations/e_ssa_threshold_review_20260930/";done
cp -- "$BASE/digitization/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* > "$OUT/inputs/source_sha256.txt"
printf 'bash %q %q\n' "$HERE/run_muon_calibration_impact_review.sh" "$LABEL" > "$OUT/inputs/command.sh"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/muon_calibration_impact.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++20 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/MuonCalibrationImpactReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/MuonCalibrationImpactReview" > "$OUT/build.log" 2>&1
"$TASK_TMP/MuonCalibrationImpactReview" "$BASE" "$OLD" "$NEW" "$OUT" > "$OUT/run.log" 2>&1
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/input_stats_after.tsv"
sha256sum -- "${SMALL_INPUTS[@]}" > "$OUT/inputs/input_sha256_after.txt"
cmp -s "$OUT/inputs/input_stats_before.tsv" "$OUT/inputs/input_stats_after.tsv"
cmp -s "$OUT/inputs/input_sha256_before.txt" "$OUT/inputs/input_sha256_after.txt"
sha256sum -c "$OUT/inputs/source_sha256.txt" > "$OUT/inputs/source_verification.log"
touch "$OUT/VERIFIED"
echo "Independent muon data calibration diagnostic complete: $OUT"
