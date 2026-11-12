#!/bin/bash
set -euo pipefail
export LC_ALL=C
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd "$HERE/../.." && pwd -P)
energy=5 max_events=6000 max_files=4 edge=10 seed=20261005 label=initial plots_only=0 intervals_only=0 maps_only=0 style_only=0
base=$PROJECT/analysis/result/shower_leak
tail_base=$PROJECT/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa
mask=$PROJECT/analysis/resolution/SingleEnergyAnalysis.cc
data_dir='' mc_dir=''
track_cut=off
source_label=''
track_args=()
usage(){ cat <<HELP
Usage: bash run.sh [--energy VALUE] [--max-events N] [--max-files N]
 [--edge-mm 10] [--seed 20261005] [--label initial] [--plots-only] [--intervals-only] [--maps-only]
 [--result-dir PATH] [--tail-results-base PATH] [--data-dir PATH] [--mc-dir PATH]
 [--style-only] [--track-cut off|on|both|v2] [--track-max-energy-gev 0.03] [--track-max-hits 60]
 [--track-min-layers 15] [--track-min-end-layers 2]
 [--track-max-hits-per-layer 2] [--track-max-rms-mm 5.3]
 [--low-signal-max-energy-gev 0.1] [--source-label LABEL --plots-only]
Defaults: 6000 events per sample from up to 4 evenly spaced files; 0 = all.
Supported energies: 0.5,1,2,3,4,5,10,20,30,40,50,60,70,80,100,120,150,200,250.
Events: one seeded random entry per equal-size stratum of the concatenated trees.
Both data and MC exclude existing low-LG tail-rich flags. v2 adds an explicit low-energy veto.
--plots-only reads saved events and mask; it does not access calibrated inputs.
--intervals-only redraws interval PNG/ROOT from saved canvases without refitting.
--maps-only redraws strip maps from saved geometry settings and mask.
--style-only redraws saved histograms and interval points without event scans or refitting.
--track-cut on vetoes sparse, low-signal through-going straight-track candidates.
This is NOT muon/pion identification. Identical thresholds apply to data and MC.
off preserves the original selection; on writes separate exclude_track results.
both aggregates once and generates off then on in the same job.
Use --plots-only --track-cut on to reuse saved events without scanning inputs.
--track-cut v2 also rejects total calibrated energy <= 0.1 GeV at beam >= 10 GeV.
At beam < 10 GeV, v2 keeps the old track cut. Output condition: exclude_track_v2.
--source-label reads events from another label and writes only to --label.
This requires --plots-only; max-events/max-files do not resample saved events.
HELP
}
args=("$@")
while (($#)); do case $1 in
 --help|-h) usage;exit 0;; --plots-only) plots_only=1;shift;; --intervals-only) plots_only=1;intervals_only=1;shift;; --maps-only) plots_only=1;maps_only=1;shift;; --style-only) plots_only=1;style_only=1;shift;;
 --*) (($#>=2))||exit 2;case $1 in
 --energy) energy=$2;; --max-events) max_events=$2;; --max-files) max_files=$2;;
 --edge-mm) edge=$2;; --seed) seed=$2;; --label) label=$2;; --result-dir) base=$2;;
 --tail-results-base) tail_base=$2;; --data-dir) data_dir=$2;; --mc-dir) mc_dir=$2;;
 --source-label) source_label=$2;; --track-cut) track_cut=$2;;
 --low-signal-max-energy-gev|--track-max-energy-gev|--track-max-hits|--track-min-layers|--track-min-end-layers|--track-max-hits-per-layer|--track-max-rms-mm) track_args+=("$1" "$2");;
 *) echo "Unknown option: $1" >&2;exit 2;;esac;shift 2;;
 *) echo "Unexpected argument $1" >&2;exit 2;; esac;done
[[ $energy =~ ^(0\.5|1|2|3|4|5|10|20|30|40|50|60|70|80|100|120|150|200|250)$ && $max_events =~ ^[0-9]+$ && $max_files =~ ^[0-9]+$ && $seed =~ ^[0-9]+$ && $label =~ ^[A-Za-z0-9_-]+$ ]]||{ echo 'Invalid options' >&2;exit 2; }
[[ $track_cut == off || $track_cut == on ||  $track_cut == both || $track_cut == v2 ]]||{ echo '--track-cut must be off, on, both or v2' >&2;exit 2; }
((style_only+intervals_only+maps_only<=1))||{ echo 'Choose only one display-only mode' >&2;exit 2; }
track_modes=("$track_cut");[[ $track_cut != both ]]||track_modes=(off on)
awk -v w="$edge" 'BEGIN{exit !(w>0&&w<110)}'||exit 2
result_root=$(realpath -m "$base")
base=$result_root/$label
if [[ -n $source_label ]];then
 [[ $source_label =~ ^[A-Za-z0-9_-]+$ && $source_label != "$label" && $plots_only == 1 ]]||{ echo '--source-label requires --plots-only and a different valid --label' >&2;exit 2; }
fi
logs=$PROJECT/analysis/log/shower_leak/$label
mkdir -p "$logs" "$base/inputs/${energy}GeV"
exec 9>"$base/inputs/${energy}GeV/.lock"
flock -n 9||{ echo 'Output is in use' >&2;exit 1; }
# Serialise builds shared by jobs for different energies.
(exec 8>"$HERE/.build.lock";flock 8;make -C "$HERE" all) > "$logs/build_${energy}GeV.log" 2>&1
record=$base/inputs/${energy}GeV
input_base=$base
if [[ -n $source_label ]];then
 input_base=$result_root/$source_label
 [[ ! -e $record/input_base.txt && ! -e $record/mask_source.cc ]]||{ echo 'Destination energy already initialized; omit --source-label to reuse it or choose a new label' >&2;exit 1; }
 [[ -f $input_base/inputs/${energy}GeV/mask_source.cc ]]||{ echo 'Missing source mask' >&2;exit 1; }
 for sample in data mc;do [[ -f $input_base/$sample/exclude/200mm/${energy}GeV/events.root ]]||{ echo 'Missing source events' >&2;exit 1; };done
elif [[ -f $record/input_base.txt ]];then
 input_base=$(cat "$record/input_base.txt")
 ((plots_only))||{ echo 'Saved-event destination requires --plots-only' >&2;exit 1; }
fi
if [[ $input_base != "$base" ]];then
 exec 7<"$input_base/inputs/${energy}GeV/.lock"
 flock -s -n 7||{ echo 'Source energy is in use' >&2;exit 1; }
fi
if [[ -n $source_label ]];then
 cp -a "$input_base/inputs/${energy}GeV" "$record/source_inputs"
 cp "$input_base/inputs/${energy}GeV/mask_source.cc" "$record/mask_source.cc"
 for sample in data mc;do cp "$input_base/$sample/exclude/200mm/${energy}GeV/inputs/files.tsv" "$record/source_${sample}_files.tsv";done
 printf '%s\n' "$input_base" > "$record/input_base.txt"
fi
plot_args=(--input-base "$input_base" --scratch-dir "$PROJECT/.scratch/shower_leak_style" --base "$base" --energy "$energy" --mask-source "$record/mask_source.cc" "${track_args[@]}")
if ((plots_only));then
 [[ -f $record/mask_source.cc ]]||{ echo 'Missing saved mask' >&2;exit 1; }
 for mode in "${track_modes[@]}";do
 condition=exclude;[[ $mode != on ]]||condition=exclude_track;[[ $mode != v2 ]]||condition=exclude_track_v2
 plot_record=$record/track_cut_$mode
 mkdir -p "$plot_record"
 for source in Plot.cc Common.hh IntervalPlots.hh TrackSelection.hh PlotStyle.hh DisplayPlots.hh run.sh submit.sh Makefile README.md;do cp "$HERE/$source" "$plot_record/$source";done
 { printf '%q ' bash "$HERE/run.sh" "${args[@]}";printf '\n';date --iso-8601=seconds;sha256sum "$HERE/Plot.cc" "$HERE/Common.hh" "$HERE/IntervalPlots.hh" "$HERE/TrackSelection.hh" "$HERE/PlotStyle.hh" "$HERE/DisplayPlots.hh" "$record/mask_source.cc" "$input_base"/{data,mc}/exclude/200mm/${energy}GeV/events.root; } >> "$plot_record/replot.log"
 "$HERE/Plot" "${plot_args[@]}" --track-cut "$mode" --intervals-only "$intervals_only" --maps-only "$maps_only" --style-only "$style_only" > "$logs/replot_${energy}GeV_track_${mode}.log" 2>&1
 printf 'Comparison: %s/comparison/%s/{200,20}mm/%sGeV/\n' "$base" "$condition" "$energy"
 done
 exit 0
fi
for sample in data mc;do [[ ! -e $base/$sample/exclude/200mm/${energy}GeV/events.root ]]||{ echo 'Result exists; choose another label or --plots-only' >&2;exit 1; };done
beam=sps;case $energy in 0.5|1|2|3|4|5) beam=ps;;esac
[[ -n $data_dir ]]||data_dir=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/$beam/calib/e-/${energy}GeV
[[ -n $mc_dir ]]||mc_dir=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold/${energy}GeV
tail=$(realpath -e "$tail_base/${energy}GeV")
# Verify the calibration parameter files recorded by the existing tail study.
# Only hashes/metadata are checked here; analysis consumes calibrated hits only.
awk '$2 ~ /\.root$/ {print}' "$tail/inputs/checksums.txt" | sha256sum -c - > "$record/tail_parameter_validation.txt"
cp "$mask" "$record/mask_source.cc"
cp "$tail/inputs/configuration.txt" "$record/tail_configuration.txt"
cp "$tail/inputs/data.tsv" "$record/tail_data.tsv"
cp "$tail/inputs/mc.tsv" "$record/tail_mc.tsv"
cp "$tail/inputs/checksums.txt" "$record/tail_checksums.txt"
for source in Common.hh IntervalPlots.hh TrackSelection.hh PlotStyle.hh DisplayPlots.hh Aggregate.cc Plot.cc run.sh submit.sh Makefile README.md;do cp "$HERE/$source" "$record/$source";done
cp "$PROJECT/analysis/resolution/HLTailSelection.hh" "$record/HLTailSelection.hh"
cp "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$record/StudyCommon.hh"
cp "$PROJECT/analysis/calibration_parameters/EBUdecode.cxx" "$record/EBUdecode.cxx"
cp "$PROJECT/analysis/calibration_parameters/EBUdecode.h" "$record/EBUdecode.h"
{ printf '%q ' bash "$HERE/run.sh" "${args[@]}";printf '\n';printf 'energy=%s\nedge_mm=%s\nmax_events=%s\nmax_files=%s\nseed=%s\ntail_results=%s\n' "$energy" "$edge" "$max_events" "$max_files" "$seed" "$tail";date --iso-8601=seconds;root-config --version;g++ --version | head -1;sha256sum "$tail/tail_study.root" "$mask" "$HERE"/*.cc "$HERE/Common.hh"; } > "$record/run.log"
for sample in data mc;do
 dir=$data_dir;[[ $sample != mc ]]||dir=$mc_dir
 [[ -d $dir ]]||{ echo "Missing input directory $dir" >&2;exit 1; }
 find "$(realpath -e "$dir")" -maxdepth 1 -type f -name '*.root' -print | sort > "$record/${sample}_available.txt"
 [[ -s $record/${sample}_available.txt ]]||{ echo 'No input ROOTs' >&2;exit 1; }
 out=$base/$sample/exclude/200mm/${energy}GeV
 "$HERE/Aggregate" --manifest "$record/${sample}_available.txt" --sample "$sample" --tail-results "$tail" --mask-source "$record/mask_source.cc" --out "$out" --energy "$energy" --edge-mm "$edge" --max-events "$max_events" --max-files "$max_files" --seed "$seed" > "$logs/${sample}_${energy}GeV.log" 2>&1
 printf 'Completed %s %s GeV: %s\n' "$sample" "$energy" "$out"
done
for mode in "${track_modes[@]}";do
condition=exclude;[[ $mode != on ]]||condition=exclude_track;[[ $mode != v2 ]]||condition=exclude_track_v2
"$HERE/Plot" "${plot_args[@]}" --track-cut "$mode" > "$logs/plot_${energy}GeV_track_${mode}.log" 2>&1
printf 'Comparison: %s/comparison/%s/{200,20}mm/%sGeV/\n' "$base" "$condition" "$energy"
done
