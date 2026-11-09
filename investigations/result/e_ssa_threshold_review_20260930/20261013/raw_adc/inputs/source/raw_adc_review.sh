#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT=$(cd -- "$HERE/../.." && pwd)
BASE=${1:-"$PROJECT/investigations/result/e_ssa_threshold_review_20260930/20261013"}
ENERGY=${2:-all}
OUT="$BASE/raw_adc"
mkdir -p -- "$PROJECT/.scratch/raw_adc_build"
BUILD=$(mktemp -d "$PROJECT/.scratch/raw_adc_build/build.XXXXXX")
trap 'rm -rf -- "$BUILD"; rmdir -- "$PROJECT/.scratch/raw_adc_build" 2>/dev/null || true' EXIT
if [[ ! -x "$HERE/RawADCReview" || "$HERE/RawADCReview.cc" -nt "$HERE/RawADCReview" || "$HERE/CurrentParameters.hh" -nt "$HERE/RawADCReview" || "$PROJECT/analysis/impact_studies/StudyCommon.hh" -nt "$HERE/RawADCReview" ]]; then
  g++ -O2 -std=c++17 $(root-config --cflags) "$HERE/RawADCReview.cc" -o "$BUILD/RawADCReview" $(root-config --libs)
  mv -- "$BUILD/RawADCReview" "$HERE/RawADCReview"
fi
mkdir -p -- "$OUT/inputs/source"
cp -- "$HERE/RawADCReview.cc" "$HERE/CurrentParameters.hh" "$HERE/raw_adc_review.sh" "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$OUT/inputs/source/"
sha256sum "$HERE/RawADCReview" "$OUT/inputs/source/"* > "$OUT/inputs/source_checksums.tsv"
"$HERE/RawADCReview" "$BASE" "$OUT" "$ENERGY" > "$OUT/execute.log" 2>&1
