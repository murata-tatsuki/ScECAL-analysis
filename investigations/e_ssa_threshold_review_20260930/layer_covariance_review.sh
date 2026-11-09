#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
[[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]] || exit 1
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
OUT=$BASE/joint_response/layer_covariance

if [[ ${2:-} == --plots-only ]]; then
 [[ -f $OUT/VERIFIED ]] || { echo 'Numerical results not verified' >&2; exit 1; }
 exec 8>"$OUT/inputs/run.lock"
 flock -n 8 || { echo 'Review already running' >&2; exit 1; }
 PLOT=$OUT/inputs/plot_source
 mkdir -p -- "$PLOT/investigations/e_ssa_threshold_review_20260930" "$PLOT/analysis/impact_studies" "$PROJECT/.scratch"
 cp -- "$HERE/"{LayerCovarianceReview,DigitizationImpactReview,ResolutionImpact}.cc "$HERE/CurrentParameters.hh" "$HERE/layer_covariance_review.sh" "$PLOT/investigations/e_ssa_threshold_review_20260930/"
 cp -- "$BASE/digitization/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$PLOT/analysis/impact_studies/"
 sha256sum -- "$OUT/"*.tsv "$OUT/layer_covariance.root" > "$OUT/inputs/replot_before.txt"
 TASK_TMP=$(mktemp -d "$PROJECT/.scratch/layer_covariance_plot.XXXXXXXX")
 trap 'rm -rf -- "$TASK_TMP"' EXIT
 read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
 read -r -a ROOT_LIBS <<< "$(root-config --libs)"
 g++ -std=c++20 -O2 "${ROOT_CFLAGS[@]}" "$PLOT/investigations/e_ssa_threshold_review_20260930/LayerCovarianceReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/plot" > "$OUT/replot_build.log" 2>&1
 "$TASK_TMP/plot" --plots-only "$OUT" > "$OUT/replot.log" 2>&1
 sha256sum -- "$OUT/"*.tsv "$OUT/layer_covariance.root" > "$OUT/inputs/replot_after.txt"
 cmp -- "$OUT/inputs/replot_before.txt" "$OUT/inputs/replot_after.txt"
 sha256sum -- "$PLOT/investigations/e_ssa_threshold_review_20260930/"* "$PLOT/analysis/impact_studies/"* > "$OUT/inputs/plot_source_checksums.txt"
 printf 'bash %q %q --plots-only\n' "$HERE/layer_covariance_review.sh" "$LABEL" > "$OUT/inputs/replot_command.sh"
 exit 0
fi

[[ -f $BASE/joint_response/digitization/COMPARISON_COMPLETE && ! -e $OUT/COMPLETE ]] || { echo 'Prerequisite incomplete or output complete' >&2; exit 1; }
mkdir -p -- "$OUT/inputs" "$PROJECT/.scratch"
exec 8>"$OUT/inputs/run.lock"
flock -n 8 || { echo 'Review already running' >&2; exit 1; }
SNAP=$OUT/inputs/source
mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies"
cp -- "$HERE/"{LayerCovarianceReview,DigitizationImpactReview,ResolutionImpact}.cc "$HERE/CurrentParameters.hh" "$HERE/layer_covariance_review.sh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
cp -- "$BASE/digitization/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
INPUTS=("$BASE/digitization/inputs/manifest.tsv" "$BASE/joint_response/digitization/inputs/manifest.tsv" "$BASE/digitization/inputs/reviewer_source/analysis/cog_channel_scan/SingleEnergyAnalysis.cc")
for energy in 5 20 100; do INPUTS+=("$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv"); done
while IFS=$'\t' read -r task energy fi oldraw oldcal truth raw cal; do [[ $task == task ]] || INPUTS+=("$cal"); done < "$BASE/digitization/inputs/manifest.tsv"
while IFS=$'\t' read -r task name energy fi oldraw oldcal truth raw cal; do [[ $name != norm_pixels19600_threshold10 ]] || INPUTS+=("$cal"); done < "$BASE/joint_response/digitization/inputs/manifest.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/checksums_before.txt"
sha256sum -- "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* > "$OUT/inputs/source_checksums.txt"
printf 'bash %q %q\n' "$HERE/layer_covariance_review.sh" "$LABEL" > "$OUT/inputs/command.sh"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/layer_covariance.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++20 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/LayerCovarianceReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/LayerCovarianceReview" > "$OUT/build.log" 2>&1
"$TASK_TMP/LayerCovarianceReview" "$BASE" "$OUT" > "$OUT/run.log" 2>&1
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/checksums_after.txt"
cmp -- "$OUT/inputs/checksums_before.txt" "$OUT/inputs/checksums_after.txt"
touch "$OUT/VERIFIED"
printf 'Completed layer covariance: %s\n' "$OUT"
