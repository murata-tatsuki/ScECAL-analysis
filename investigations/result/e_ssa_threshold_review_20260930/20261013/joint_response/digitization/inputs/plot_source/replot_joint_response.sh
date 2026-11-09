#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P);PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013};[[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||exit 1
OUT=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL/joint_response/digitization
[[ -f $OUT/COMPARISON_COMPLETE ]]||{ echo 'Comparison incomplete' >&2;exit 1; }
mkdir -p -- "$OUT/inputs/plot_source"
cp -- "$HERE/JointResponseReplot.cc" "$HERE/replot_joint_response.sh" "$OUT/inputs/plot_source/"
sha256sum "$OUT/inputs/plot_source/"* > "$OUT/inputs/plot_source_checksums.txt"
printf 'bash %q %q\n' "$HERE/replot_joint_response.sh" "$LABEL" > "$OUT/inputs/replot_command.sh"
INPUTS=("$OUT/"*.tsv);for energy in 5 20 100;do INPUTS+=("$OUT/${energy}GeV/joint_response.root");done
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/replot_science_before_sha256.txt"
TMP=$(mktemp -d "$PROJECT/.scratch/joint_response_replot.XXXXXXXX");trap 'rm -rf -- "$TMP"' EXIT
read -r -a FLAGS <<< "$(root-config --cflags)";read -r -a LIBS <<< "$(root-config --libs)"
g++ -std=c++17 -O2 "${FLAGS[@]}" "$OUT/inputs/plot_source/JointResponseReplot.cc" "${LIBS[@]}" -o "$TMP/JointResponseReplot" > "$OUT/logs/build_replot.log" 2>&1
"$TMP/JointResponseReplot" "$OUT" > "$OUT/logs/replot.log" 2>&1
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/replot_science_after_sha256.txt"
cmp -s "$OUT/inputs/replot_science_before_sha256.txt" "$OUT/inputs/replot_science_after_sha256.txt"
touch "$OUT/REPLOT_COMPLETE"
echo 'Completed read-only PNG replot'
