#!/bin/bash
# Replot saved tail rates only; does not read raw/calibrated event trees.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energy=100; mc_tag=threshold; label=default
result_dir="$HERE/../../result/channel_response/hg_lg_tail"
while (($#)); do
 case "$1" in
  -h|--help) echo 'Usage: bash plot_layers.sh [--energy 50] [--mc-tag threshold] [--label default] [--result-dir PATH]'; exit 0;;
  --energy|--mc-tag|--label|--result-dir)
   (($#>=2)) || { echo "Missing value: $1" >&2;exit 1; }
   case "$1" in --energy) energy=$2;; --mc-tag) mc_tag=$2;; --label) label=$2;; --result-dir) result_dir=$2;; esac;shift 2;;
  *) echo "Unknown option: $1" >&2;exit 1;;
 esac
done
[[ $energy =~ ^[0-9]+([.][0-9]+)?$ && $mc_tag =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $label =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo 'Invalid energy/tag/label' >&2;exit 1; }
out="$result_dir/$mc_tag/$label/${energy}GeV"
[[ -f $out/tail_study.root && -f $out/channels.tsv ]] || { echo "Completed result required: $out" >&2;exit 1; }
( flock 9; make -C "$HERE" PlotHGLGTailLayers ) 9>"$HERE/.build.lock"
exec 8>"$out/.run.lock"
flock -n 8 || { echo 'Another process is writing this result' >&2;exit 1; }
"$HERE/PlotHGLGTailLayers" "$out"
