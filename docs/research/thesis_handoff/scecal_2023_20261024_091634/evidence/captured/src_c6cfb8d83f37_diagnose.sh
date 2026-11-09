#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
if (($#!=3));then echo 'Usage: bash diagnose.sh SOURCE_RESULT LABEL CellID,...' >&2;exit 1;fi
source_dir=$(realpath -e -- "$1");label=$2;channels=$3
[[ $label =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && $channels =~ ^[0-9]+(,[0-9]+)*$ ]]||exit 1
out=$PROJECT/analysis/result/mip_calibration/$label
[[ -f $source_dir/COLLECTED && ! -e $out ]]||{ echo 'Completed collection and fresh output label required' >&2;exit 1; }
(flock 9;make -C "$HERE" Diagnose) 9>"$HERE/.build.lock"
exec 7<"$source_dir/.run.lock";flock -sn 7||exit 1
mkdir -p "$out/inputs"
exec 8>"$out/.run.lock";flock -n 8||exit 1
exec > >(tee "$out/run.log") 2>&1
printf '%q ' "$HERE/diagnose.sh" "$@" > "$out/inputs/command.sh";printf '\n' >> "$out/inputs/command.sh"
awk -F '\t' 'BEGIN{OFS="\t"} $1=="max-relative-error"{$2="0.10"} {print}' "$source_dir/inputs/configuration.tsv" > "$out/inputs/configuration.tsv"
awk -F '\t' 'BEGIN{OFS="\t"} $1=="min-fit-hits"{$2="50"} {print}' "$out/inputs/configuration.tsv" > "$out/inputs/run_configuration.tsv"
ln -s "$source_dir/inputs" "$out/inputs/collection_origin"
cp "$HERE/Diagnose" "$out/inputs/Diagnose"
tar -czf "$out/inputs/analysis_sources.tar.gz" -C "$PROJECT" analysis/mip_calibration/{Diagnose.cc,diagnose.sh,FitModel.hh,Common.hh,Makefile} analysis/channel_response/muon_mip_response/{MuonMIPCommon.hh,CurrentCalibration.hh,TrackGeometry.hh} analysis/channel_response/adc_energy/FastPng.hh
sha256sum "$source_dir/histograms.root" "$out/inputs/Diagnose" "$out/inputs/analysis_sources.tar.gz" > "$out/inputs/checksums.sha256"
root-config --version > "$out/inputs/root_version.txt"
"$out/inputs/Diagnose" "$source_dir/histograms.root" "$out/inputs/configuration.tsv" "$out" "$channels"
sha256sum "$out/diagnostics.root" "$out/models.tsv" "$out/profiles.tsv" "$out/runs.tsv" > "$out/inputs/result_checksums.sha256"
touch "$out/COMPLETE"
echo "Completed diagnosis: $out"
