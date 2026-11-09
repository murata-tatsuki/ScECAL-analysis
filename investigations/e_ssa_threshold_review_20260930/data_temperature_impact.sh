#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
[[ $# -le 1 && $LABEL =~ ^[A-Za-z0-9_-]+$ ]] || { echo 'Usage: data_temperature_impact.sh [LABEL]' >&2;exit 1; }
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
OUT=$BASE/digitization/temperature_estimate/data_calibration
[[ ! -e $OUT/summary.tsv ]] || { echo 'Scientific output already exists' >&2;exit 1; }
mkdir -p -- "$OUT/inputs/source/investigations/e_ssa_threshold_review_20260930" "$OUT/inputs/source/analysis/impact_studies" "$OUT/inputs/source/analysis/cog_channel_scan" "$PROJECT/.scratch"
exec 8>"$OUT/inputs/run.lock"
flock -n 8 || { echo 'Another reviewer is writing' >&2;exit 1; }
SRC=$OUT/inputs/source
cp -- "$HERE/"{DataTemperatureImpact.cc,data_temperature_impact.sh,DigitizationImpactReview.cc,ResolutionImpact.cc,CurrentParameters.hh} "$SRC/investigations/e_ssa_threshold_review_20260930/"
cp -- "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$SRC/analysis/impact_studies/"
cp -- "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SRC/analysis/cog_channel_scan/"
INPUTS=("$BASE/digitization/inputs/reviewer_parameters/"*.root)
for beam in 5 20 100;do
 INPUTS+=("$BASE/integrated/comparison/20mm/${beam}GeV/events.tsv" "$BASE/${beam}GeV/baseline/inputs/data.tsv")
 while IFS=$'\t' read -r raw cal;do INPUTS+=("$raw" "$cal");done < "$BASE/${beam}GeV/baseline/inputs/data.tsv"
done
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/input_stats_before.tsv"
sha256sum "$SRC/investigations/e_ssa_threshold_review_20260930/"* "$SRC/analysis/impact_studies/"* "$SRC/analysis/cog_channel_scan/"* > "$OUT/inputs/source_checksums.txt"
printf 'bash %q %q\n' "$HERE/data_temperature_impact.sh" "$LABEL" > "$OUT/inputs/command.sh"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/data_temperature_impact.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++17 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SRC/investigations/e_ssa_threshold_review_20260930/DataTemperatureImpact.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/DataTemperatureImpact" > "$OUT/build.log" 2>&1
"$TASK_TMP/DataTemperatureImpact" "$BASE" "$OUT" > "$OUT/run.log" 2>&1
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/input_stats_after.tsv"
cmp -s -- "$OUT/inputs/input_stats_before.tsv" "$OUT/inputs/input_stats_after.tsv"
touch "$OUT/COMPLETE"
printf 'Completed %s\n' "$OUT"
