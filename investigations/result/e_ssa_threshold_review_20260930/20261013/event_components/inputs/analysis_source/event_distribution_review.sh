#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
label=${1:-20261013}
[[ $label != *[^a-zA-Z0-9_-]* && -n $label ]] || { echo 'Invalid label' >&2; exit 1; }
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$label"
input=${2:-"$base/integrated"}
input=$(realpath -- "$input")
out="$base/event_components"
[[ ! -e "$out/distribution_summary.tsv" ]] || { echo "Output exists: $out" >&2; exit 1; }
for energy in 5 20 100; do [[ -s "$input/comparison/20mm/${energy}GeV/events.tsv" ]]; done
mkdir -p "$out/inputs/analysis_source"
exec > >(tee "$out/run.log") 2>&1
printf 'Started: %s\nCommand: bash %q %q %q\nInput: %s\nOutput: %s\n' "$(date --iso-8601=seconds)" "$HERE/event_distribution_review.sh" "$label" "$input" "$input" "$out"
build="$PROJECT/.scratch/event_distribution_review_build"
mkdir -p "$build"
temporary=$(mktemp "$build/EventDistributionReview.XXXXXX")
cleanup() { [[ -z ${temporary:-} ]] || rm -f -- "$temporary"; rmdir "$build" 2>/dev/null || :; }
trap cleanup EXIT HUP INT TERM
c++ -std=c++17 -O2 $(root-config --cflags) "$HERE/EventDistributionReview.cc" $(root-config --libs) -lASImage -o "$temporary"
chmod 755 "$temporary"
mv -f -- "$temporary" "$HERE/EventDistributionReview"
temporary=
cleanup
cp -p -- "$HERE/EventDistributionReview.cc" "$HERE/event_distribution_review.sh" "$PROJECT/analysis/channel_response/adc_energy/FastPng.hh" "$out/inputs/analysis_source/"
if [[ -f "$input/inputs/checksums.txt" ]]; then cp -p -- "$input/inputs/checksums.txt" "$out/inputs/integrated_review_checksums.txt"; fi
inputs=("$input/covariance.tsv")
for energy in 5 20 100; do inputs+=("$input/comparison/20mm/${energy}GeV/events.tsv"); done
sha256sum "$HERE/EventDistributionReview.cc" "$HERE/event_distribution_review.sh" "$HERE/EventDistributionReview" "$PROJECT/analysis/channel_response/adc_energy/FastPng.hh" "${inputs[@]}" > "$out/inputs/checksums_before.txt"
"$HERE/EventDistributionReview" "$input" "$out"
sha256sum "$HERE/EventDistributionReview.cc" "$HERE/event_distribution_review.sh" "$HERE/EventDistributionReview" "$PROJECT/analysis/channel_response/adc_energy/FastPng.hh" "${inputs[@]}" > "$out/inputs/checksums_after.txt"
cmp "$out/inputs/checksums_before.txt" "$out/inputs/checksums_after.txt"
printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
