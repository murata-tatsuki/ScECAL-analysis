#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
project=$(cd -- "$here/../../.." && pwd -P)
base="$project/analysis/result/channel_response/hl_intercept"
reference=/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/share/all_hl_electron2023_hlratio_v3.root
energy=100
fit_options=()
while (($#)); do
 case "$1" in
  --base) base=$2; shift 2;;
  --reference) reference=$2; shift 2;;
  --energy) energy=$2; shift 2;;
  --low|--high|--min-bin|--min-bins|--min-span|--min-hits|--pivot|--switch-margin|--gain-switch|--top|--channels) fit_options+=("$1" "$2");shift 2;;
  --help) cat <<'HELP'
Usage: bash run.sh [--energy 100] [--base DIRECTORY] [--reference ROOT] [fit options]
Reads BASE/{before_correction,after_correction}/{data,simulation/threshold}/ENERGYGeV.root.
Creates BASE/InterCalib_{data,before,after}_ENERGYGeV.root and BASE/comparison/ENERGYGeV/.
Existing outputs are never overwritten. Use another --base containing the input directories
(or symlinks to them) for alternative fit ranges or repeat studies.
Fit options: --low 800 --high 2200 --min-bin 20 --min-bins 8 --min-span 400
 --min-hits 200 --pivot 1500 --switch-margin 600 --gain-switch 2600 --top 20 [--channels 240008,420008]
--channels is intended for development/small studies; default covers all 6300 channels.
HELP
   exit 0;;
  *) echo "Unknown option: $1" >&2; exit 2;;
 esac
done
[[ $energy =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo 'Invalid energy' >&2;exit 2; }
base=$(cd -- "$base" && pwd -P)
outputs=("InterCalib_data_${energy}GeV.root" "InterCalib_before_${energy}GeV.root" "InterCalib_after_${energy}GeV.root" "comparison/${energy}GeV")
for name in "${outputs[@]}"; do [[ ! -e "$base/$name" ]] || { echo "Refusing to overwrite: $base/$name" >&2;exit 1; }; done
make -C "$here" CompareHLIntercept
lock="$base/.comparison_${energy}GeV.lock"
mkdir "$lock" || { echo "Another comparison may be running: $lock" >&2;exit 1; }
stage=$(mktemp -d "$base/.comparison_${energy}GeV.XXXXXX")
cleanup() { local status=$?; rmdir "$lock"; if ((status));then echo "Comparison failed; staging files retained at $stage" >&2;else rm -rf -- "$stage";fi; }
trap cleanup EXIT
args=(--data "$base/after_correction/data/${energy}GeV.root"
 --data-before "$base/before_correction/data/${energy}GeV.root"
 --before "$base/before_correction/simulation/threshold/${energy}GeV.root"
 --after "$base/after_correction/simulation/threshold/${energy}GeV.root"
 --reference "$reference" --energy "$energy" --output "$stage/result" "${fit_options[@]}")
printf '%q ' "$here/CompareHLIntercept" "${args[@]}" > "$stage/command.txt"
printf '\n' >> "$stage/command.txt"
"$here/CompareHLIntercept" "${args[@]}" 2>&1 | tee "$stage/run.log"
cp "$stage/command.txt" "$stage/run.log" "$stage/result/comparison/${energy}GeV/"
stat -c '%n	%s bytes	%y' "$base/after_correction/data/${energy}GeV.root" "$base/before_correction/data/${energy}GeV.root" "$base/before_correction/simulation/threshold/${energy}GeV.root" "$base/after_correction/simulation/threshold/${energy}GeV.root" > "$stage/result/comparison/${energy}GeV/input_file_stats.txt"
sha256sum "$here/CompareHLIntercept.cc" "$here/CompareHLIntercept" "$here/run.sh" "$reference" > "$stage/result/comparison/${energy}GeV/checksums.txt"
# Recheck after computation. Only complete results are published.
for name in "${outputs[@]}"; do [[ ! -e "$base/$name" ]] || { echo "Output appeared during analysis: $base/$name" >&2;exit 1; }; done
mkdir -p -- "$base/comparison"
for name in "${outputs[@]}"; do mv -T -- "$stage/result/$name" "$base/$name";done
cat "$base/comparison/${energy}GeV/summary.txt"
