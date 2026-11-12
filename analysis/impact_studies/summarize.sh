#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
result_dir=$HERE/../result/impact_studies
mc_tag=threshold label=default study=hl_intercept
energies_csv=''
while (($#));do
 case $1 in
  -h|--help) echo 'Usage: bash summarize.sh [--result-dir PATH] [--mc-tag NAME] [--label NAME] [--study hl_intercept] [--energies LIST]';exit 0;;
  --energies) (($#>=2))||exit 1;energies_csv=$2;shift 2;;
  --result-dir|--mc-tag|--label|--study) (($#>=2))||exit 1;case $1 in --result-dir) result_dir=$2;; --mc-tag) mc_tag=$2;; --label) label=$2;; --study) study=$2;; esac;shift 2;;
  *) echo "Unknown option: $1" >&2;exit 1;;
 esac
done
[[ $study == hl_intercept && $mc_tag =~ ^[A-Za-z0-9_-]+$ && $label =~ ^[A-Za-z0-9_-]+$ ]]||exit 1
base=$(realpath -m -- "$result_dir")/$study/$mc_tag/$label
shopt -s nullglob
files=("$base"/*GeV/study.root)
if [[ -n $energies_csv ]];then
 [[ $energies_csv =~ ^[0-9]+(\.[0-9]+)?(,[0-9]+(\.[0-9]+)?)*$ ]]||{ echo 'Invalid energy list' >&2;exit 1; }
 IFS=, read -r -a energies <<< "$energies_csv"
 files=()
 for energy in "${energies[@]}";do
  file=$base/${energy}GeV/study.root
  [[ -f $file ]]||{ echo "Requested energy is incomplete: $energy GeV ($file)" >&2;exit 1; }
  files+=("$file")
 done
fi
((${#files[@]}))||{ echo "No completed results in $base" >&2;exit 1; }
make -C "$HERE" SummarizeImpact
"$HERE/SummarizeImpact" "$base/energy_scan" "${files[@]}"
