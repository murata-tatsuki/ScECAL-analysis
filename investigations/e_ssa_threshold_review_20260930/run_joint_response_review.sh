#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}; [[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||exit 1
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
PILOT=$BASE/digitization; OUT=$BASE/joint_response/digitization
[[ -f $OUT/inputs/manifest.tsv && ! -e $OUT/COMPARISON_COMPLETE ]]||{ echo 'Joint manifest incomplete or comparison exists' >&2;exit 1; }
exec 8>"$OUT/inputs/comparison.lock"; flock -n 8||{ echo 'Joint review already running' >&2;exit 1; }
sha256sum -c "$OUT/inputs/source_parameter_checksums.txt" > "$OUT/inputs/source_verification.log"
INPUTS=("$OUT/inputs/manifest.tsv" "$PILOT/inputs/manifest.tsv")
while IFS=$'\t' read -r task case_name energy file_index old_raw old_cal truth raw cal;do
 [[ $task == task ]]&&continue
 [[ -f $OUT/mc/$case_name/${energy}GeV/${raw##*/}.COMPLETE ]]||{ echo "Joint output incomplete: $raw" >&2;exit 1; }
 INPUTS+=("$old_raw" "$old_cal" "$truth" "$raw" "$cal")
done < "$OUT/inputs/manifest.tsv"
while IFS=$'\t' read -r task energy file_index old_raw old_cal truth raw cal;do [[ $task == task ]]&&continue;INPUTS+=("$raw" "$cal");done < "$PILOT/inputs/manifest.tsv"
for energy in 5 20 100;do INPUTS+=("$PILOT/comparison/20mm/${energy}GeV/paired_events.tsv" "$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv" "$BASE/${energy}GeV/baseline/reference.tsv");done
INPUTS+=("$PILOT/inputs/reviewer_parameters/"*.root)
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_stats_before.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_sha256_before.txt"
SNAP=$OUT/inputs/reviewer_source
mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan"
cp -- "$HERE/"{JointResponseReview,DigitizationImpactReview,ResolutionImpact}.cc "$HERE/CurrentParameters.hh" "$HERE/run_joint_response_review.sh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
cp -- "$PILOT/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
cp -- "$PILOT/inputs/reviewer_source/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* > "$OUT/inputs/reviewer_code_checksums.txt"
printf 'bash %q %q\n' "$HERE/run_joint_response_review.sh" "$LABEL" > "$OUT/inputs/review_command.sh"
mkdir -p -- "$PROJECT/.scratch"; TASK_TMP=$(mktemp -d "$PROJECT/.scratch/joint_response_review.XXXXXXXX"); trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"; read -r -a ROOT_LIBS <<< "$(root-config --libs)"
for driver in DigitizationImpactReview JointResponseReview;do
 g++ -std=c++20 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/$driver.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/$driver" > "$OUT/logs/build_review_$driver.log" 2>&1
done
for case_name in muon_norm norm_threshold10 norm_pixels19600 norm_pixels19600_threshold10;do
 CASE_OUT=$OUT/comparison/$case_name; mkdir -p -- "$CASE_OUT/inputs"
 ln -s -- "$PILOT/inputs/reviewer_parameters" "$CASE_OUT/inputs/reviewer_parameters"
 ln -s -- "$SNAP" "$CASE_OUT/inputs/reviewer_source"
 awk -F '\t' -v name="$case_name" 'BEGIN{OFS="\t";print "task","energy_GeV","file_index","old_decode","old_calib","truth","new_decode","new_calib"} NR>1&&$2==name{print $1,$3,$4,$5,$6,$7,$8,$9}' "$OUT/inputs/manifest.tsv" > "$CASE_OUT/inputs/manifest.tsv"
 awk -F '\t' -v name="$case_name" 'NR==1||$1==name' "$OUT/inputs/cases.tsv" > "$CASE_OUT/inputs/model_condition.tsv"
 args=("$BASE" "$CASE_OUT"); [[ $case_name != muon_norm ]]&&args+=(100)
 "$TASK_TMP/DigitizationImpactReview" "${args[@]}" > "$OUT/logs/review_$case_name.log" 2>&1
 touch "$CASE_OUT/COMPARISON_COMPLETE"
done
"$TASK_TMP/JointResponseReview" "$BASE" "$OUT" > "$OUT/logs/review_summary.log" 2>&1
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_stats_after.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_sha256_after.txt"
cmp -s "$OUT/inputs/comparison_stats_before.tsv" "$OUT/inputs/comparison_stats_after.tsv"
cmp -s "$OUT/inputs/comparison_sha256_before.txt" "$OUT/inputs/comparison_sha256_after.txt"
touch "$OUT/COMPARISON_COMPLETE"
printf 'Completed joint ADC/hit response review: %s\n' "$OUT"
