#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT=$(cd -- "$HERE/../.." && pwd)
BASE=${1:-"$PROJECT/investigations/result/e_ssa_threshold_review_20260930/20261013"}
ENERGY=${2:-all}
MODE=${3:-}
[[ -z "$MODE" || "$MODE" == --plots-only ]] || { echo 'Usage: raw_adc_review.sh [RESULT_BASE] [5|20|100|all] [--plots-only]' >&2; exit 2; }
OUT="$BASE/raw_adc"
if [[ -z "$MODE" ]]; then
  for BEAM in 5 20 100; do
    [[ "$ENERGY" == all || "$ENERGY" == "$BEAM" ]] || continue
    [[ ! -e "$OUT/comparison/20mm/${BEAM}GeV/histograms.root" ]] || { echo "Existing scientific output for ${BEAM}GeV; use --plots-only or a fresh result base." >&2; exit 3; }
  done
fi
mkdir -p -- "$PROJECT/.scratch/raw_adc_build"
BUILD=$(mktemp -d "$PROJECT/.scratch/raw_adc_build/build.XXXXXX")
trap 'rm -rf -- "$BUILD"; rmdir -- "$PROJECT/.scratch/raw_adc_build" 2>/dev/null || true' EXIT
if [[ ! -x "$HERE/RawADCReview" || "$HERE/RawADCReview.cc" -nt "$HERE/RawADCReview" || "$HERE/CurrentParameters.hh" -nt "$HERE/RawADCReview" || "$PROJECT/analysis/impact_studies/StudyCommon.hh" -nt "$HERE/RawADCReview" ]]; then
  g++ -O2 -std=c++17 $(root-config --cflags) "$HERE/RawADCReview.cc" -o "$BUILD/RawADCReview" $(root-config --libs)
  mv -- "$BUILD/RawADCReview" "$HERE/RawADCReview"
fi
if [[ "$MODE" == --plots-only ]]; then
  mkdir -p -- "$OUT/inputs/plot_source"
  cp -- "$HERE/RawADCReview.cc" "$HERE/raw_adc_review.sh" "$OUT/inputs/plot_source/"
  sha256sum "$HERE/RawADCReview" "$OUT/inputs/plot_source/"* > "$OUT/inputs/plot_source_checksums.tsv"
  rg --files "$OUT" | rg '\.(root|tsv)$' | rg -v '/inputs/' | sort | xargs -d '\n' sha256sum > "$OUT/inputs/scientific_checksums.sha256"
  "$HERE/RawADCReview" --plots-only "$OUT" "$ENERGY" > "$OUT/repaint.log" 2>&1
  sha256sum -c "$OUT/inputs/scientific_checksums.sha256" > "$OUT/inputs/repaint_validation.log"
else
  mkdir -p -- "$OUT/inputs/source"
  cp -- "$HERE/RawADCReview.cc" "$HERE/CurrentParameters.hh" "$HERE/raw_adc_review.sh" "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$OUT/inputs/source/"
  sha256sum "$HERE/RawADCReview" "$OUT/inputs/source/"* > "$OUT/inputs/source_checksums.tsv"
  "$HERE/RawADCReview" "$BASE" "$OUT" "$ENERGY" > "$OUT/execute.log" 2>&1
fi
