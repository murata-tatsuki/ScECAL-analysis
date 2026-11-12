#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
energy=all label=default tag=threshold result=$HERE/../result/impact_studies
while (($#));do
 case $1 in
  -h|--help) printf '%s\n' 'Usage: bash replot.sh [--energy all|VALUE] [--label NAME] [--mc-tag NAME] [--result-dir PATH]' 'Redraw event energy/change PNGs in completed results. all selects completed energies.' 'Original study.root and channel figures are preserved; new histogram counts go to event_histograms.root.';exit 0;;
  --energy|--label|--mc-tag|--result-dir) (($#>=2))||{ echo "Missing value for $1" >&2;exit 1; };case $1 in --energy) energy=$2;;--label) label=$2;;--mc-tag) tag=$2;;--result-dir) result=$2;;esac;shift 2;;
  *) echo "Unknown option $1" >&2;exit 1;;
 esac
done
[[ $label =~ ^[A-Za-z0-9_-]+$ && $tag =~ ^[A-Za-z0-9_-]+$ ]]||exit 1
[[ $energy == all || $energy =~ ^[0-9]+(\.[0-9]+)?$ ]]||exit 1
base=$(realpath -m -- "$result")/hl_intercept/$tag/$label
dirs=()
if [[ $energy == all ]];then
 shopt -s nullglob
 for file in "$base"/*GeV/study.root;do dirs+=("${file%/study.root}");done
else dirs=("$base/${energy}GeV");fi
((${#dirs[@]}))||{ echo "No completed studies: $base" >&2;exit 1; }
for dir in "${dirs[@]}";do [[ -f $dir/study.root ]]||{ echo "Study incomplete: $dir" >&2;exit 1; };done
make -C "$HERE" ReplotEventEnergy
stage=''
trap '[[ -z $stage ]] || rm -rf -- "$stage"' EXIT
for dir in "${dirs[@]}";do
 exec 9>"$dir/.lock"
 flock -n 9||{ echo "Study in use: $dir" >&2;exit 1; }
 stage=$(mktemp -d "$dir/.replot.XXXXXX")
 "$HERE/ReplotEventEnergy" --input "$dir" --out "$stage"
 mkdir -p "$dir/figures"
 for file in "$stage/figures/"*.png;do mv -f -- "$file" "$dir/figures/";done
 for name in event_histograms.root event_histogram_binning.tsv summary.tsv;do mv -f -- "$stage/$name" "$dir/$name";done
 rm -rf -- "$stage";stage=''
 exec 9>&-
done
