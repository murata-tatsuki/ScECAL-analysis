#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
base=${1:?Usage: inspect.sh INPUT_LABEL_DIRECTORY OUTPUT_DIRECTORY [ENERGIES_CSV]}
out=${2:?Output directory required}
energies=${3:-0.5,1,5,30,100,200,250}
params=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
args=(--hl "$params/all_hl_electron2023_hlratio_v3.root" --mip "$params/all_auto_muon_v4_trackfit.root" --pedestal "$params/pedestal2023_SPS.root")
IFS=, read -ra values <<< "$energies"
for energy in "${values[@]}";do
 count=$(wc -l < "$base/${energy}GeV/inputs/mc.tsv")
 mc_per_file=$(( (20000+count-1)/count ))
 data_per_file=2000
 case $energy in 0.5|1|2|3|4|5) data_per_file=100000;; esac
 "$HERE/InspectResponse" --input "$base/${energy}GeV" --out "$out/${energy}GeV" "${args[@]}" --data-per-file "$data_per_file" --mc-per-file "$mc_per_file" --comparison-source "$HERE/../../analysis/cog_channel_scan/SingleEnergyAnalysis.cc"
done
