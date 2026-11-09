#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
energy=${1:-100}
[[ $energy =~ ^[0-9]+([.][0-9]+)?$ ]] || exit 1
out=${2:-$PROJECT/investigations/result/raw_calibration/${energy}GeV}
make -C "$HERE" InspectTails
"$HERE/InspectTails" "$out" "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$PROJECT/analysis/result/impact_studies/hl_intercept/threshold/production_v1/${energy}GeV/study.root" 2>&1 | tee "$out/tails.log"
