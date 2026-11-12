#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../../.." && pwd -P)
MC=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
PARAMS=$MC/Analysis_edit/share
original=("$@")
fail(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
usage(){ cat <<HELP
Usage: bash run.sh [options]
 --label NAME                 New result label (default default); never overwritten
 --reuse-collection LABEL    Reuse completed collection under a NEW label; never overwrite failed results
 --reuse-fits LABEL           Reuse collected input and completed fit rows under a NEW label
 --histogram-storage MODE     compact (default; exact sparse ROOT table) or legacy TH1 objects
 --fit-workers N              Fit processes (1; Slurm CPUs-1 when CPUs>1); parent writes all ROOT/TSV
 --trend-axes standard|all     standard omits timestamp_ticks response trends; timestamp audit retained
 --memory-mb N                Aggregate/sort memory budget (512; range 16..8192), plus ROOT overhead
 --input-manifest TSV         7 columns: sample(0/1), particle(muon/electron/pedestal), run,
                              absolute decode, calib(- for pedestal), truth(- for data), evidence(- if unavailable)
 --energy GEV                 Electron energy (100); does not change muon energy
 --muon-energy GEV            Muon energy (100), actual files required
 --electron-runs LIST         Data electron runs (all)
 --muon-runs LIST             Muon runs (all)
 --max-files N                Per sample/particle (1; 0 all)
 --max-events N               Contiguous raw events per file (1000; 0 all)
 --start-entry N              First raw entry (0)
 --event-bin N                Fixed raw-entry interval width (250)
 --time-bin SEC               Fixed elapsed-time interval width (60)
 --tick-bin N                 Saved timestamp tick interval width (60)
 --seconds-per-tick S         Default 0: physical time unconfirmed, ticks/event order only
 --clock-evidence TEXT        Required for nonzero seconds-per-tick; no epoch is inferred
 --temperature-bin C          Common temperature bin width (0.25)
 --block-events N             Correlated error / train-validation block width (250)
 --cog MM                     Electron layer9/10 square CoG cut (200=nocut)
 --data-tail 0|1 --mc-tail 0|1 Default electron: data exclude(1), MC keep(0); muon unaffected
 --tail-results DIR           Existing hg_lg_tail/threshold/e_ssa/<energy>GeV
 --gain-table TSV             Existing sipm_gain gain_estimates.tsv (optional, preserve fit statuses)
 --mask FILE                  Physical CellIDs, one per line (default none)
 --data-base DIR --mc-base DIR
 --data-pedestal FILE --mc-pedestal FILE --data-mip FILE --mc-mip FILE
 --data-hl FILE --mc-hl FILE --threshold FILE
 --min-layers N --max-angle RAD --max-track-rms MM --long-margin MM --trans-margin MM
 --hg-max ADC                 Raw HG cap (2200)
 --min-fit-hits N             Plateau entries (100)
 --fit-chi2max X              Landau-Gaussian fit quality (5)
 --min-slope-hits N --min-temperature-span C --min-residual-hits N --min-plot-hits N
 --channels all|CELLID,...    Plot selection only; all channels collected and fit
 --diagnose-correction 0|1    Optional held-out electron time model (default 0)
 --split blocks|runs          Whole alternating blocks / even-vs-odd runs
 --bootstrap N                Paired validation block bootstrap replicates (200)
 --collect-only               Collect once without statistics/plots
 --analyze-only               Analyze completed collection at same label, immutable saved settings
 --plots-only                 Saved tables only, no collect/refit, new plot directory
 --plot-label NAME            Render label (initial; plots-only requires a new name)
 --no-plots                   Skip PNG
 --dry-run                    Validate paths/settings, no build/write/job
HELP
}
declare -A o=( [label]=default [reuse-collection]='' [reuse-fits]='' [histogram-storage]=compact [trend-axes]=standard [fit-workers]=1 [memory-mb]=512 [input-manifest]='' [energy]=100 [muon-energy]=100 [electron-runs]=all [muon-runs]=all [max-files]=1 [max-events]=1000 [start-entry]=0 [event-bin]=250 [time-bin]=60 [tick-bin]=60 [seconds-per-tick]=0 [clock-evidence]=- [temperature-bin]=0.25 [block-events]=250 [cog]=200 [data-tail]=1 [mc-tail]=0 [tail-results]='' [mask]=- [gain-table]=- [data-base]=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps [mc-base]="$MC/Result_MC" [data-pedestal]="$PARAMS/pedestal2023_SPS.root" [mc-pedestal]="$PARAMS/pedestal2023_SPS.root" [data-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [mc-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [data-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root" [mc-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root" [threshold]="$PROJECT/analysis/result/threshold/threshold.root" [min-layers]=7 [max-angle]=0.05 [max-track-rms]=3 [long-margin]=1 [trans-margin]=0.25 [hg-max]=2200 [min-fit-hits]=100 [fit-chi2max]=5 [min-slope-hits]=40 [min-temperature-span]=0.5 [min-residual-hits]=10 [min-plot-hits]=40 [channels]=all [diagnose-correction]=0 [split]=blocks [bootstrap]=200 [plot-label]=initial )
declare -A given=()
dry=0; mode=all; plots=1
while (($#));do case $1 in
 --help|-h)usage;exit 0;; --dry-run)dry=1;shift;; --no-plots)plots=0;shift;;
 --collect-only|--analyze-only|--plots-only)[[ $mode == all ]]||fail 'Conflicting modes';mode=${1#--};shift;;
 --*)k=${1#--};[[ -v o[$k] && $# -ge 2 ]]||fail "Unknown/missing option $1";o[$k]=$2;given[$k]=1;shift 2;;
 *)fail "Unexpected argument $1";;esac;done
for k in label plot-label;do [[ ${o[$k]} =~ ^[a-zA-Z0-9][a-zA-Z0-9_.-]*$ && ${o[$k]} != .. ]]||fail "Invalid $k";done
[[ ${o[memory-mb]} =~ ^[0-9]+$ ]] && ((o[memory-mb]>=16&&o[memory-mb]<=8192)) || fail 'memory-mb must be 16..8192'
[[ ${o[histogram-storage]} == compact || ${o[histogram-storage]} == legacy ]]||fail 'histogram-storage must be compact or legacy'
[[ ${o[trend-axes]} == standard || ${o[trend-axes]} == all ]]||fail 'trend-axes must be standard or all'
if [[ -v given[tick-bin] && ! -v given[trend-axes] ]];then o[trend-axes]=all;fi
if [[ ! -v given[fit-workers] && ${SLURM_CPUS_PER_TASK:-1} -gt 1 ]];then o[fit-workers]=$((SLURM_CPUS_PER_TASK-1));fi
[[ ${o[fit-workers]} =~ ^[0-9]+$ ]]&&((o[fit-workers]>=1&&o[fit-workers]<=32))||fail 'fit-workers must be 1..32'
if ((o[fit-workers]>1));then available_cpus=${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN)};((o[fit-workers]+1<=available_cpus))||fail 'Fit workers plus parent exceed allocated CPUs';fi
FIT_SOURCE=-
if [[ -n ${o[reuse-fits]} ]];then
 [[ -z ${o[reuse-collection]} || ${o[reuse-collection]} == "${o[reuse-fits]}" ]]||fail 'reuse-fits and reuse-collection must use the same label'
 o[reuse-collection]=${o[reuse-fits]}
fi
fit_definition(){
 sed -n '/^struct Summary /,/^};/p; /^void count(/p; /^template<int N>struct Histogram/,/^};/p; /^double langau(/,/^struct Tables /p; /^ for(Long64_t i=0;i<nHits/,/^ }ht->ResetBranchAddresses();/p' "$1" | sed '/^#include /d; /^bool includeTickTrends=/d; /^std::vector<Axis> selectedAxes(/d; s/selectedAxes(/axes(/g'
}
OUT=$PROJECT/analysis/result/channel_response/time_stability/${o[label]}
if [[ -n ${o[reuse-collection]} ]];then
 [[ $mode == all ]]||fail 'reuse-collection cannot be combined with another mode'
 [[ ${o[reuse-collection]} =~ ^[a-zA-Z0-9][a-zA-Z0-9_.-]*$ && ${o[reuse-collection]} != .. ]]||fail 'Invalid source label'
 for k in "${!given[@]}";do case $k in label|reuse-collection|reuse-fits|histogram-storage|trend-axes|fit-workers|memory-mb|channels|plot-label);;*)fail "reuse-collection uses saved collection settings; --$k not permitted";;esac;done
 SOURCE=$PROJECT/analysis/result/channel_response/time_stability/${o[reuse-collection]}
 [[ -f $SOURCE/COLLECTED && -f $SOURCE/inputs/configuration.tsv && -f $SOURCE/inputs/manifest.tsv && -f $SOURCE/inputs/collection.sha256 && -f $SOURCE/flat.root ]]||fail 'Completed source collection missing'
 if [[ -n ${o[reuse-fits]} ]];then
  [[ -s $SOURCE/fits.tsv && -f $SOURCE/inputs/statistics_sources/Analyze.cc && -f $SOURCE/inputs/statistics_sources/Common.hh && -f $SOURCE/inputs/statistics_sources/StreamingStats.hh ]]||fail 'Saved fit rows and numerical definitions missing'
  cmp -s "$SOURCE/inputs/statistics_sources/Common.hh" "$HERE/Common.hh"||fail 'Cached fit axis/input definitions changed'
  cmp -s "$SOURCE/inputs/statistics_sources/StreamingStats.hh" "$HERE/StreamingStats.hh"||fail 'Cached fit numerical definitions changed'
  cmp -s <(fit_definition "$SOURCE/inputs/statistics_sources/Analyze.cc") <(fit_definition "$HERE/Analyze.cc")||fail 'Cached histogram or fit definitions changed'
  FIT_SOURCE=$SOURCE/fits.tsv
  printf 'Reuse completed attempted fits: %s\n' "$FIT_SOURCE"
 fi
 [[ ! -e $OUT ]]||fail "Result label already exists: $OUT"
 printf 'Reuse collection: %s\nNew output: %s\nMemory budget: %s MiB plus ROOT overhead\n' "$SOURCE" "$OUT" "${o[memory-mb]}"
 (( !dry ))||exit 0
 exec 7<"$SOURCE/.run.lock";flock -sn 7||fail 'Source result is being written'
 sha256sum -c "$SOURCE/inputs/collection.sha256" > /dev/null
 mkdir -p -- "$(dirname -- "$OUT")";mkdir -- "$OUT";mkdir -- "$OUT/inputs"
 cp -- "$SOURCE/inputs/configuration.tsv" "$SOURCE/inputs/manifest.tsv" "$OUT/inputs/"
 cp -- "$SOURCE/inputs/collection.sha256" "$OUT/inputs/source_collection.sha256"
 printf 'source_result\t%s\nflat_file\t%s\n' "$SOURCE" "$(realpath -e -- "$SOURCE/flat.root")" > "$OUT/inputs/reused_collection.tsv"
 ln -s -- "$(realpath -e -- "$SOURCE/flat.root")" "$OUT/flat.root"
 for f in calibration.tsv counters.tsv input_audit.tsv production_audit.tsv root_metadata.tsv run_boundaries.tsv;do [[ ! -f $SOURCE/$f ]]||cp -- "$SOURCE/$f" "$OUT/$f";done
 sha256sum "$OUT/flat.root" "$OUT/inputs/configuration.tsv" "$OUT/inputs/manifest.tsv" > "$OUT/inputs/collection.sha256"
 mkdir -- "$OUT/inputs/execution_sources"
 cp -- "$HERE/run.sh" "$HERE/submit.sh" "$HERE/Makefile" "$HERE/README.md" "$OUT/inputs/execution_sources/"
 if [[ $FIT_SOURCE != - ]];then
  printf 'source_fits\t%s\nsource_collection_checksum\t%s\n' "$FIT_SOURCE" "$SOURCE/inputs/collection.sha256" > "$OUT/inputs/reused_fits.tsv"
  stat -Lc '%n\t%s\t%y\t%i' "$FIT_SOURCE" > "$OUT/inputs/reused_fits_stat_before.tsv"
  fit_definition "$HERE/Analyze.cc" | sha256sum > "$OUT/inputs/fit_definition.sha256"
 fi
 touch "$OUT/COLLECTED";mode=analyze-only
elif [[ $mode == plots-only || $mode == analyze-only ]];then
 [[ -f $OUT/COLLECTED && -f $OUT/inputs/configuration.tsv ]]||fail 'Completed collection missing'
 for k in "${!given[@]}";do case $k in label|channels|plot-label|memory-mb);;histogram-storage|trend-axes|fit-workers)[[ $mode == analyze-only ]]||fail "plots-only uses saved axes and statistics; --$k not permitted";;*)fail "$mode uses saved settings; --$k not permitted";;esac;done
 if [[ $mode == plots-only ]];then [[ -f $OUT/ANALYZED ]]||fail 'Completed statistics missing';[[ ! -e $OUT/plots/${o[plot-label]} ]]||fail 'Choose a new --plot-label';else [[ ! -e $OUT/statistics.root ]]||fail 'Statistics already exist';fi
 printf '%s: %s; render=%s\n' "$mode" "$OUT" "${o[plot-label]}"
 (( !dry ))||exit 0
else
 for k in max-files max-events start-entry event-bin block-events min-layers min-fit-hits min-slope-hits min-residual-hits min-plot-hits bootstrap data-tail mc-tail diagnose-correction;do [[ ${o[$k]} =~ ^[0-9]+$ ]]||fail "Invalid integer $k";done
 for k in energy muon-energy time-bin tick-bin seconds-per-tick temperature-bin cog hg-max fit-chi2max min-temperature-span max-angle max-track-rms long-margin trans-margin;do [[ ${o[$k]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k";done
 awk -v e="${o[energy]}" -v m="${o[muon-energy]}" -v t="${o[time-bin]}" -v z="${o[tick-bin]}" -v b="${o[temperature-bin]}" -v hg="${o[hg-max]}" -v angle="${o[max-angle]}" -v rms="${o[max-track-rms]}" -v lm="${o[long-margin]}" -v tm="${o[trans-margin]}" 'BEGIN{exit !(e>0&&m>0&&t>0&&z>0&&b>0&&hg>0&&hg<=2600&&angle>0&&angle<.5&&rms>0&&lm<22.5&&tm<2.5)}'||fail 'Invalid ranges'
 ((o[event-bin]>0&&o[block-events]>0&&o[min-layers]>=3&&o[min-layers]<=14&&o[min-fit-hits]>=20&&o[bootstrap]>=20&&o[data-tail]<=1&&o[mc-tail]<=1&&o[diagnose-correction]<=1))||fail 'Invalid counts'
 [[ ${o[split]} == blocks || ${o[split]} == runs ]]||fail 'Invalid split'
 [[ ${o[channels]} == all || ${o[channels]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid channels'
 [[ ${o[seconds-per-tick]} == 0 || ${o[clock-evidence]} != - ]]||fail 'Supply --clock-evidence for time unit'
 [[ ! -e $OUT ]]||fail "Result label already exists: $OUT"
 [[ -n ${o[tail-results]} ]]||o[tail-results]=$PROJECT/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/${o[energy]}GeV
 for k in data-pedestal mc-pedestal data-mip mc-mip data-hl mc-hl threshold;do [[ -f ${o[$k]} ]]||fail "Missing $k ${o[$k]}";done
 [[ ${o[mask]} == - || -f ${o[mask]} ]]||fail 'Missing mask'
 [[ ${o[gain-table]} == - || -f ${o[gain-table]} ]]||fail 'Missing gain table'
 rows=()
 if [[ -n ${o[input-manifest]} ]];then
 [[ -f ${o[input-manifest]} ]]||fail 'Missing manifest';while IFS= read -r line;do [[ -z $line || $line == \#* || $line == sample$'\t'* ]]||rows+=("$line");done < "${o[input-manifest]}"
 else
 shopt -s nullglob
 for particle in muon electron;do
 if [[ $particle == muon ]];then species=mu-;energy=${o[muon-energy]};tag=mu_track;wanted=${o[muon-runs]};else species=e-;energy=${o[energy]};tag=e_ssa;wanted=${o[electron-runs]};fi
 count=0
 for raw in "${o[data-base]}/decode/$species/${energy}GeV/"*.root;do
 name=${raw##*/};[[ $name =~ ^ECAL_Run([0-9]+) ]]||continue;run=${BASH_REMATCH[1]};[[ $wanted == all || ,$wanted, == *,$run,* ]]||continue
 rows+=("0"$'\t'"$particle"$'\t'"$run"$'\t'"$raw"$'\t'"${o[data-base]}/calib/$species/${energy}GeV/$name"$'\t-\t-');((++count));if ((o[max-files]>0&&count>=o[max-files]));then break;fi
 done
 ((count))||fail "No actual $particle data at ${energy}GeV"
 count=0
 for raw in "${o[mc-base]}/decode/$species/$tag/threshold/${energy}GeV/"*.root;do
 name=${raw##*/};run=$count
 if [[ $particle == muon ]];then [[ $name =~ ^ECAL_Run([0-9]+) ]]||continue;run=${BASH_REMATCH[1]};[[ $wanted == all || ,$wanted, == *,$run,* ]]||continue;else run=${name#${energy}GeV_};run=${run%.root};fi
 rows+=("1"$'\t'"$particle"$'\t'"$run"$'\t'"$raw"$'\t'"${o[mc-base]}/calib/$species/$tag/threshold/${energy}GeV/$name"$'\t'"${o[mc-base]}/generate/$species/$tag/${energy}GeV/$name"$'\t-');((++count));if ((o[max-files]>0&&count>=o[max-files]));then break;fi
 done
 ((count))||fail "No actual $particle MC at ${energy}GeV; use manifest, never electron MC for muons"
 done
 fi
 ((${#rows[@]}))||fail 'Empty input selection'
 declare -A seen=()
 for row in "${rows[@]}";do IFS=$'\t' read -r sample particle run raw cal truth evidence extra <<< "$row";[[ -z ${extra:-} && -n ${evidence:-} ]]||fail 'Manifest requires 7 columns';[[ $sample == 0 || $sample == 1 ]]||fail 'Invalid sample';[[ $particle == muon || $particle == electron || $particle == pedestal ]]||fail 'Invalid particle';[[ $run =~ ^[0-9]+$ ]]||fail 'Invalid run';[[ $particle == pedestal || $cal != - ]]||fail 'Calibration required';[[ $sample == 0 || $truth != - ]]||fail 'MC truth required'
 for path in "$raw" "$cal" "$truth" "$evidence";do [[ $path == - || ( $path == /* && -f $path ) ]]||fail "Input must be an existing absolute path: $path";done
 canonical=$(realpath -e -- "$raw");[[ ! -v seen[$canonical] ]]||fail 'Duplicate raw file';seen[$canonical]=1
 if [[ $particle == electron && ( ( $sample == 0 && ${o[data-tail]} == 1 ) || ( $sample == 1 && ${o[mc-tail]} == 1 ) ) ]];then [[ -f ${o[tail-results]}/tail_study.root ]]||fail 'Missing existing tail results';fi
 done
 printf 'Output: %s\nContiguous events/file: %s; event-bin=%s; seconds-per-tick=%s\n' "$OUT" "${o[max-events]}" "${o[event-bin]}" "${o[seconds-per-tick]}"
 printf '%s\n' "${rows[@]}"
 (( !dry ))||exit 0
 mkdir -p -- "$(dirname -- "$OUT")";mkdir -- "$OUT";mkdir -- "$OUT/inputs"
 printf '%s\n' "${rows[@]}" > "$OUT/inputs/manifest.tsv"
 for k in "${!o[@]}";do printf '%s\t%s\n' "$k" "${o[$k]}";done | sort > "$OUT/inputs/configuration.tsv"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0
exec 8>"$OUT/.run.lock";flock -n 8||fail 'Output is in use'
(flock 9;make -C "$HERE" all) 9>"$HERE/.build.lock"
{
 printf 'cd %q\n' "$PWD";printf '%q ' "$HERE/run.sh" "${original[@]}";printf '\n'
} >> "$OUT/run.log"
if [[ $mode == all || $mode == collect-only ]];then
 mkdir -- "$OUT/inputs/parameters" "$OUT/inputs/sources"
 printf 'parameter\tsource\n' > "$OUT/inputs/parameter_sources.tsv"
 for k in data-pedestal mc-pedestal data-mip mc-mip data-hl mc-hl threshold;do printf '%s\t%s\n' "$k" "${o[$k]}" >> "$OUT/inputs/parameter_sources.tsv";cp -- "${o[$k]}" "$OUT/inputs/parameters/$k.root";o[$k]=$OUT/inputs/parameters/$k.root;done
 if [[ ${o[gain-table]} != - ]];then printf '%s\n' "${o[gain-table]}" > "$OUT/inputs/gain_source.txt";cp -- "${o[gain-table]}" "$OUT/inputs/gain_estimates.tsv";o[gain-table]=$OUT/inputs/gain_estimates.tsv;fi
 for k in "${!o[@]}";do printf '%s\t%s\n' "$k" "${o[$k]}";done | sort > "$OUT/inputs/configuration.tsv"
 cp -- "$HERE"/*.{cc,hh,sh} "$HERE/Makefile" "$OUT/inputs/sources/"
 for source in "$HERE/../muon_mip_response/"{MuonMIPCommon,TrackGeometry,CurrentCalibration}.hh "$HERE/../../resolution/HLTailSelection.hh" /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/{Decode,Calibration,EBUdecode}.cxx "$MC/Analysis_edit/src/MCDigi.cxx";do cp -- "$source" "$OUT/inputs/sources/";done
 [[ ${o[mask]} == - ]]||cp -- "${o[mask]}" "$OUT/inputs/mask.txt"
 if [[ -f ${o[tail-results]}/tail_study.root ]];then sha256sum "${o[tail-results]}/tail_study.root" "${o[tail-results]}/inputs/"{data,mc}.tsv > "$OUT/inputs/tail_checksums.sha256";fi
 sha256sum "$OUT/inputs/parameters/"* "$OUT/inputs/sources/"* > "$OUT/inputs/code_constants.sha256"
 stats(){ while IFS=$'\t' read -r sample particle run raw cal truth evidence;do for path in "$raw" "$cal" "$truth" "$evidence";do [[ $path == - ]]||stat -Lc '%n\t%s\t%y\t%i' -- "$path";done;done < "$OUT/inputs/manifest.tsv"; }
 stats > "$OUT/inputs/input_stats_before.tsv"
 "$HERE/Collect" "$OUT/inputs/manifest.tsv" "$OUT/inputs/configuration.tsv" "$OUT" 2>&1 | tee -a "$OUT/run.log"
 stats > "$OUT/inputs/input_stats_after.tsv";cmp "$OUT/inputs/input_stats_before.tsv" "$OUT/inputs/input_stats_after.tsv"||fail 'Inputs changed during collection'
 sha256sum "$OUT/flat.root" "$OUT/inputs/configuration.tsv" "$OUT/inputs/manifest.tsv" > "$OUT/inputs/collection.sha256"
 touch "$OUT/COLLECTED"
fi
[[ $mode != collect-only ]]||exit 0
if [[ $mode != plots-only ]];then
 mkdir -- "$OUT/inputs/statistics_sources"
 cp -- "$HERE/Analyze.cc" "$HERE/Common.hh" "$HERE/DiskReduce.hh" "$HERE/StreamingStats.hh" "$HERE/HistogramStore.hh" "$HERE/FitReuse.hh" "$HERE/FitPool.hh" "$OUT/inputs/statistics_sources/"
 sha256sum "$OUT/inputs/statistics_sources/"* > "$OUT/inputs/statistics_code.sha256"
 sha256sum -c "$OUT/inputs/collection.sha256" > "$OUT/inputs/collection_check.log"
 printf 'histogram_storage\t%s\nfit_reuse\t%s\ntrend_axes\t%s\nfit_workers\t%s\n' "${o[histogram-storage]}" "$FIT_SOURCE" "${o[trend-axes]}" "${o[fit-workers]}" > "$OUT/inputs/statistics_execution.tsv"
 "$HERE/Analyze" "$OUT/inputs/configuration.tsv" "$OUT" "${o[memory-mb]}" "${o[histogram-storage]}" "$FIT_SOURCE" "${o[trend-axes]}" "${o[fit-workers]}" 2>&1 | tee -a "$OUT/run.log"
 if [[ $FIT_SOURCE != - ]];then stat -Lc '%n\t%s\t%y\t%i' "$FIT_SOURCE" > "$OUT/inputs/reused_fits_stat_after.tsv";cmp "$OUT/inputs/reused_fits_stat_before.tsv" "$OUT/inputs/reused_fits_stat_after.tsv"||fail 'Source fit cache changed';fi
 sha256sum "$OUT/statistics.root" "$OUT/"*.tsv > "$OUT/inputs/statistics.sha256"
 touch "$OUT/ANALYZED"
fi
if ((plots));then
 sha256sum -c "$OUT/inputs/statistics.sha256" > /dev/null
 "$HERE/Draw" "$OUT/inputs/configuration.tsv" "$OUT" "$OUT/plots/${o[plot-label]}" "${o[channels]}" "${o[memory-mb]}" 2>&1 | tee -a "$OUT/run.log"
 sha256sum -c "$OUT/inputs/statistics.sha256" > "$OUT/plots/${o[plot-label]}/numerical_checksum_check.log"
 cp -- "$HERE/Draw.cc" "$HERE/Common.hh" "$HERE/DiskReduce.hh" "$HERE/../adc_energy/FastPng.hh" "$OUT/plots/${o[plot-label]}/"
 printf '%q ' "$HERE/run.sh" "${original[@]}" > "$OUT/plots/${o[plot-label]}/command.sh";printf '\n' >> "$OUT/plots/${o[plot-label]}/command.sh"
 sha256sum "$HERE/Draw.cc" "$HERE/Common.hh" "$HERE/DiskReduce.hh" "$HERE/../adc_energy/FastPng.hh" > "$OUT/plots/${o[plot-label]}/code.sha256"
fi
[[ $mode == plots-only ]]||touch "$OUT/COMPLETE"
