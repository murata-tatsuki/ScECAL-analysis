#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
PARAMS=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
DATA_ANALYSIS=${DATA_ANALYSIS:-/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO}
original=("$@")
fail(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
usage(){ cat <<HELP
Usage: bash run.sh [options]
  --label NAME                  New result label (default: default; no overwrites)
  --runs all|40,44               Default: all data muon runs
  --max-events N                 Evenly spaced events PER FILE; 0=all (default)
  --channels all|CELLID,...      All physical channels by default; tracking uses all layers
  --input-manifest FILE          run<TAB>absolute_decode_path; optional header run/path
  --decode-dir DIR               Default: 2023/sps/decode/mu-/100GeV
  --pedestal ROOT --current-mip ROOT --hl ROOT --threshold ROOT
  --temperature-mode current|none   Production convention (default) / diagnostic
  --min-layers N --max-angle RAD --max-track-rms MM (7 / 0.05 / 3)
  --long-margin MM --trans-margin MM (1 / 0.25; full strip chord)
  --factor-step VALUE            Correction stratum width (0.005; max 0.02)
  --min-fit-hits N               Counts in fit interval (200)
  --max-relative-error VALUE --max-stability VALUE (0.10 / 0.05)
  --max-correlation VALUE        MPV-width correlation limit (0.98)
  --fit-min-adc ADC             Raw pedestal-subtracted ADC fit floor (40); checked at +/-20
  --fit-max-10 ADC --fit-max-15 ADC (1000 / 2000)
  --workers N                   Fit processes (4; separate shard ROOTs)
  --refit-from RESULT_DIR        Reuse completed histograms; new label required
  --plots-only                  Redraw existing label from saved fits, no fitting
  --chip-plots-only             Add/redraw only 6x6 chip PNGs; implies --plots-only
  --collect-only                Stop after histogram collection
  --no-plots                    Retain numerical outputs without PNG generation
  --dry-run                     Validate configuration and paths; no writes/build/job
  --help
Output: analysis/result/mip_calibration/<label>/
Logs: analysis/log/mip_calibration/
No electron analysis, MC matching, or production parameter replacement.
HELP
}
declare -A opt=( [label]=default [runs]=all [max-events]=0 [channels]=all [input-manifest]=''
 [decode-dir]=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps/decode/mu-/100GeV
 [pedestal]="$PARAMS/pedestal2023_SPS.root" [current-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root" [threshold]="$PROJECT/analysis/result/threshold/threshold.root"
 [temperature-mode]=current [min-layers]=7 [max-angle]=0.05 [max-track-rms]=3 [long-margin]=1 [trans-margin]=0.25 [factor-step]=0.005
 [min-fit-hits]=200 [max-relative-error]=0.10 [max-stability]=0.05 [max-correlation]=0.98 [fit-min-adc]=40 [fit-max-10]=1000 [fit-max-15]=2000 [workers]=4 [refit-from]='' [plots-only]=0 [chip-plots-only]=0 [collect-only]=0 [no-plots]=0 )
declare -A given=()
dry=0
while (($#));do case $1 in
 --help|-h)usage;exit 0;; --dry-run)dry=1;shift;;
 --chip-plots-only)opt[plots-only]=1;given[plots-only]=1;opt[chip-plots-only]=1;given[chip-plots-only]=1;shift;;
 --plots-only|--collect-only|--no-plots)opt[${1#--}]=1;given[${1#--}]=1;shift;;
 --*)key=${1#--};[[ -v opt[$key] && $# -ge 2 ]]||fail "Unknown option or missing value: $1";opt[$key]=$2;given[$key]=1;shift 2;;
 *)fail "Unexpected argument $1";;esac;done
if [[ -n ${opt[refit-from]} ]];then
 [[ ${opt[plots-only]} == 0 && ${opt[collect-only]} == 0 ]]||fail 'refit mode cannot collect or be plots-only'
 source_dir=$(realpath -e -- "${opt[refit-from]}")
 [[ -f $source_dir/COLLECTED && -f $source_dir/histograms.root && -f $source_dir/inputs/configuration.tsv ]]||fail 'refit requires completed collection'
 for key in "${!given[@]}";do case $key in label|channels|temperature-mode|min-fit-hits|max-relative-error|max-stability|max-correlation|fit-min-adc|fit-max-10|fit-max-15|workers|refit-from|no-plots);;*)fail "--$key changes collection; recollect under a new label";;esac;done
 while IFS=$'\t' read -r key value;do
  [[ -v opt[$key] ]]||continue
  case $key in label|workers|refit-from|plots-only|chip-plots-only|collect-only|no-plots|max-relative-error)continue;;esac
  [[ -v given[$key] ]]||opt[$key]=$value
 done < "$source_dir/inputs/configuration.tsv"
fi
[[ ${opt[label]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]||fail 'Invalid label'
[[ ${opt[runs]} == all || ${opt[runs]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid run list'
[[ ${opt[channels]} == all || ${opt[channels]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid channel list'
[[ ${opt[temperature-mode]} == current || ${opt[temperature-mode]} == none ]]||fail 'Unknown temperature mode'
for key in max-events min-layers min-fit-hits workers;do [[ ${opt[$key]} =~ ^(0|[1-9][0-9]*)$ ]]||fail "Invalid integer $key";done
((opt[workers]>0&&opt[workers]<=64&&opt[min-layers]>=3&&opt[min-layers]<=14&&opt[min-fit-hits]>=50))||fail 'Invalid worker/statistics/track count'
for key in max-angle max-track-rms long-margin trans-margin factor-step max-relative-error max-stability max-correlation fit-min-adc fit-max-10 fit-max-15;do [[ ${opt[$key]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $key";done
awk -v a="${opt[max-angle]}" -v r="${opt[max-track-rms]}" -v l="${opt[long-margin]}" -v t="${opt[trans-margin]}" -v f="${opt[factor-step]}" -v e="${opt[max-relative-error]}" -v s="${opt[max-stability]}" -v fmin="${opt[fit-min-adc]}" -v c="${opt[max-correlation]}" 'BEGIN{exit !(a>0&&a<.5&&r>0&&l<22.5&&t<2.5&&f>0&&f<=.02&&e>0&&e<1&&s>0&&s<1&&c>0&&c<1&&fmin>=0&&fmin<500)}'||fail 'Invalid numerical limit'
[[ -z ${SLURM_CPUS_PER_TASK:-} || ${opt[workers]} -le $SLURM_CPUS_PER_TASK ]]||fail 'workers exceeds allocated Slurm CPUs'
out=$PROJECT/analysis/result/mip_calibration/${opt[label]}
if ((opt[plots-only]));then
 for key in "${!given[@]}";do case $key in plots-only|chip-plots-only|label|channels|workers);;*)fail "plots-only cannot change --$key";;esac;done
 [[ -f $out/COMPLETE && -f $out/fits.root ]]||fail 'plots-only requires completed fits'
 printf 'Redraw only: %s\n' "$out";((!dry))||exit 0
 (flock 9;make -C "$HERE" -j "${opt[workers]}" all) 9>"$HERE/.build.lock"
 exec 8>"$out/.run.lock";flock -n 8||fail 'Output is in use'
 sha256sum "$out/histograms.root" "$out/fits.root" "$out/mip.root" > "$out/inputs/plot_checksums_before.txt"
 tar -czf "$out/inputs/plot_sources.tar.gz" -C "$PROJECT" analysis/mip_calibration/{Draw.cc,Common.hh,run.sh,Makefile} analysis/channel_response/adc_energy/FastPng.hh
 sha256sum "$HERE/Draw" "$out/inputs/plot_sources.tar.gz" > "$out/inputs/plot_code_checksums.sha256"
 printf '%q ' "$HERE/run.sh" "${original[@]}" >> "$out/inputs/plot_commands.sh";printf '\n' >> "$out/inputs/plot_commands.sh"
 draw_options=();((!opt[chip-plots-only]))||draw_options+=(--chip-plots-only)
 "$HERE/Draw" "$out/histograms.root" "$out/fits.root" "$out" "${opt[channels]}" "${draw_options[@]}" 2>&1 | tee -a "$out/run.log"
 sha256sum "$out/histograms.root" "$out/fits.root" "$out/mip.root" > "$out/inputs/plot_checksums_after.txt"
 cmp "$out/inputs/plot_checksums_before.txt" "$out/inputs/plot_checksums_after.txt"
 exit 0
fi
[[ ! -e $out ]]||fail "Output already exists: $out; use a new label or --plots-only"
for key in pedestal current-mip hl threshold;do [[ -f ${opt[$key]} ]]||fail "Missing $key: ${opt[$key]}";done
manifest_lines=()
include_run(){ [[ ${opt[runs]} == all || ,${opt[runs]}, == *,"$1",* ]]; }
if [[ -n ${opt[refit-from]} ]];then
 while IFS= read -r line;do manifest_lines+=("$line");done < "$source_dir/inputs/manifest.tsv"
elif [[ -n ${opt[input-manifest]} ]];then
 [[ -f ${opt[input-manifest]} ]]||fail 'Missing manifest'
 while IFS=$'\t' read -r run path extra;do
  [[ -z $run || $run == \#* || $run == run ]]&&continue
  [[ $run =~ ^[0-9]+$ && $path == /* && -f $path && -z $extra ]]||fail 'Malformed manifest row'
  include_run "$run"&&manifest_lines+=("$run"$'\t'"$path")
 done < "${opt[input-manifest]}"
else
 shopt -s nullglob
 for path in "${opt[decode-dir]}"/ECAL_Run*.root;do name=${path##*/};[[ $name =~ ^ECAL_Run([0-9]+)_ ]]||fail "Cannot infer run: $path";run=${BASH_REMATCH[1]};include_run "$run"&&manifest_lines+=("$run"$'\t'"$(realpath -e -- "$path")");done
fi
((${#manifest_lines[@]}>0))||fail 'No selected input files'
printf 'Output: %s\nInput files: %s; max events/file: %s; workers: %s; temperature: %s\n' "$out" "${#manifest_lines[@]}" "${opt[max-events]}" "${opt[workers]}" "${opt[temperature-mode]}"
((!dry))||exit 0
(flock 9;make -C "$HERE" -j "${opt[workers]}" all) 9>"$HERE/.build.lock"
mkdir -p -- "$out/inputs/parameters" "$out/inputs/bin" "$PROJECT/analysis/log/mip_calibration/${opt[label]}"
exec 8>"$out/.run.lock";flock -n 8||fail 'Output is in use'
exec > >(tee -a "$out/run.log") 2>&1
printf 'Started: %s\n' "$(date -Is)"
printf '%q ' "$HERE/run.sh" "${original[@]}" > "$out/inputs/command.sh";printf '\n' >> "$out/inputs/command.sh"
printf '%s\n' "${manifest_lines[@]}" > "$out/inputs/manifest.tsv"
for key in pedestal current-mip hl threshold;do
 printf '%s\t%s\n' "$key" "${opt[$key]}" >> "$out/inputs/parameter_sources.tsv"
 cp -- "${opt[$key]}" "$out/inputs/parameters/$key.root";opt[$key]=$out/inputs/parameters/$key.root
done
for key in "${!opt[@]}";do printf '%s\t%s\n' "$key" "${opt[$key]}";done | LC_ALL=C sort > "$out/inputs/configuration.tsv"
cp -- "$HERE/Collect" "$HERE/Fit" "$HERE/Export" "$HERE/Draw" "$out/inputs/bin/"
tar -czf "$out/inputs/analysis_sources.tar.gz" -C "$PROJECT" analysis/mip_calibration/{Common.hh,Collect.cc,FitModel.hh,Fit.cc,Export.cc,Draw.cc,Diagnose.cc,diagnose.sh,Makefile,run.sh,submit.sh,README.md} analysis/channel_response/muon_mip_response/{MuonMIPCommon.hh,TrackGeometry.hh,CurrentCalibration.hh} analysis/channel_response/{root.mk,adc_energy/FastPng.hh}
mkdir "$out/inputs/current_reader_sources"
for source in "$DATA_ANALYSIS"/{src/Calibration.cxx,src/EBUdecode.cxx,include/Calibration.h,include/EBUdecode.h,run/global_config} "$PARAMS/../src/Extract.cxx" "$PARAMS/../include/Extract.h" "$PARAMS/../src/MCDigi.cxx";do cp -- "$source" "$out/inputs/current_reader_sources/";done
sha256sum "$out"/inputs/parameters/* "$out"/inputs/bin/* "$out"/inputs/current_reader_sources/* "$out/inputs/analysis_sources.tar.gz" > "$out/inputs/checksums.sha256"
root-config --version > "$out/inputs/root_version.txt"
if [[ -n ${opt[refit-from]} ]];then
 exec 7<"$source_dir/.run.lock";flock -sn 7||fail 'Source result is in use'
 ln -s -- "$source_dir/histograms.root" "$out/histograms.root"
 ln -s -- "$source_dir/dependencies.tsv" "$out/dependencies.tsv"
 ln -s -- "$source_dir/event_audit.tsv" "$out/event_audit.tsv"
 ln -s -- "$source_dir/inputs" "$out/inputs/collection_origin"
 sha256sum "$source_dir/histograms.root" > "$out/inputs/refit_input.sha256"
else
 while IFS=$'\t' read -r run path;do stat --printf '%n\t%s\t%y\n' -- "$path";done < "$out/inputs/manifest.tsv" > "$out/inputs/input_stats_before.tsv"
 "$out/inputs/bin/Collect" "$out/inputs/manifest.tsv" "$out/inputs/configuration.tsv" "$out"
 while IFS=$'\t' read -r run path;do stat --printf '%n\t%s\t%y\n' -- "$path";done < "$out/inputs/manifest.tsv" > "$out/inputs/input_stats_after.tsv"
 cmp "$out/inputs/input_stats_before.tsv" "$out/inputs/input_stats_after.tsv"
fi
touch "$out/COLLECTED"
((!opt[collect-only]))||exit 0
mkdir -p -- "$PROJECT/.scratch/mip_calibration_runs"
work=$(mktemp -d "$PROJECT/.scratch/mip_calibration_runs/${opt[label]}.XXXXXX")
printf 'Worker workspace: %s\n' "$work"
trap 'status=$?; if ((status));then printf "Incomplete run; retained worker files: %s\n" "$work";fi' EXIT
pids=();shards=()
for ((i=0;i<opt[workers];++i));do
 shard=$work/worker$i.root;shards+=("$shard")
 "$out/inputs/bin/Fit" "$out/histograms.root" "$out/inputs/configuration.tsv" "$shard" "$i" "${opt[workers]}" > "$PROJECT/analysis/log/mip_calibration/${opt[label]}/fit-$i.log" 2>&1 & pids+=("$!")
done
failed=0;for pid in "${pids[@]}";do wait "$pid"||failed=1;done
((!failed))||fail 'A fit worker failed; see analysis/log/mip_calibration'
hadd -v 0 "$out/fits.root" "${shards[@]}"
for type in tsv variants.tsv;do
 destination=$out/fits.tsv;[[ $type == tsv ]]||destination=$out/fit_variants.tsv
 awk 'FNR==1&&NR!=1{next}{print}' "${shards[@]/%/.$type}" > "$destination"
done
"$out/inputs/bin/Export" "$out/histograms.root" "$out/fits.root" "$out"
if ((!opt[no-plots]));then "$out/inputs/bin/Draw" "$out/histograms.root" "$out/fits.root" "$out" "${opt[channels]}";fi
sha256sum "$out/histograms.root" "$out/fits.root" "$out/mip.root" > "$out/inputs/result_checksums.sha256"
touch "$out/COMPLETE"
for shard in "${shards[@]}";do rm -f -- "$shard" "$shard.tsv" "$shard.variants.tsv";done
rmdir -- "$work";rmdir "$PROJECT/.scratch/mip_calibration_runs" 2>/dev/null || :
trap - EXIT
printf 'Completed: %s\nReader status: %s/reader_status.tsv\n' "$(date -Is)" "$out"
