#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
label=${1:-20261013}
mode=${2:-}
[[ -z $mode || $mode == --recompute ]] || { echo "Invalid mode: $mode" >&2; exit 1; }
[[ $label != *[^a-zA-Z0-9_-]* && -n $label ]] || { echo 'Invalid label' >&2; exit 1; }
base="$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$label"
out="$base/digitization/temperature_estimate"
[[ ! -e "$out/summary.tsv" || $mode == --recompute ]] || { echo "Output exists: $out" >&2; exit 1; }
inputs=()
for energy in 5 20 100; do
  inputs+=("$base/integrated/comparison/20mm/${energy}GeV/events.tsv" "$base/digitization/comparison/20mm/${energy}GeV/paired_events.tsv")
done
for file in "${inputs[@]}"; do [[ -s $file ]]; done
mkdir -p "$out/inputs/analysis_source"
exec > >(tee -a "$out/run.log") 2>&1
printf 'Command: bash %q %q %q\nInput: %s\nOutput: %s\n' "$HERE/temperature_resolution_estimate.sh" "$label" "$mode" "$base" "$out"
build="$PROJECT/.scratch/temperature_resolution_estimate"
mkdir -p "$build"
binary=$(mktemp "$build/TemperatureResolutionEstimate.XXXXXX")
cleanup() { rm -f -- "$binary"; rmdir "$build" 2>/dev/null || :; }
trap cleanup EXIT HUP INT TERM
c++ -std=c++17 -O2 $(root-config --cflags) "$HERE/TemperatureResolutionEstimate.cc" $(root-config --libs) -o "$binary"
cp -p -- "$HERE/TemperatureResolutionEstimate.cc" "$HERE/temperature_resolution_estimate.sh" "$out/inputs/analysis_source/"
sha256sum "${inputs[@]}" "$HERE/TemperatureResolutionEstimate.cc" "$HERE/temperature_resolution_estimate.sh" > "$out/inputs/checksums_before.txt"
"$binary" "$base" "$out"
sha256sum "${inputs[@]}" "$HERE/TemperatureResolutionEstimate.cc" "$HERE/temperature_resolution_estimate.sh" > "$out/inputs/checksums_after.txt"
cmp "$out/inputs/checksums_before.txt" "$out/inputs/checksums_after.txt"
printf 'Completed: %s\n' "$(date --iso-8601=seconds)"
