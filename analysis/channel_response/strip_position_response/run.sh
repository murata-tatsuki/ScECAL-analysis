#!/bin/bash
set -euo pipefail
original=("$@")
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../../.." && pwd -P)
MC=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
PARAMS=$MC/Analysis_edit/share
fail(){ echo "ERROR: $*" >&2; exit 1; }
usage(){ cat <<HELP
Usage: bash run.sh [options]
Muon ADC strip-position analysis; no production constants updated.
  --label NAME                  Fresh result label (default: pilot)
  --runs 44|40,44|all            Real muon runs (default: 44)
  --energy GeV                  Actual muon input energy (default: 100)
  --sample both|data|mc          Default: both
  --mc-tag threshold            Explicit MC production directory
  --input-manifest PATH          Six TSV columns: 0/1,run,decode,calib,-,truth
  --data-base PATH --mc-base PATH
  --data-pedestal PATH --mc-pedestal PATH
  --data-mip PATH --mc-mip PATH --data-hl PATH --mc-hl PATH
  --data-threshold PATH --mc-threshold PATH
  --max-events N                Evenly spaced per file (default: 1000; 0=all)
  --channels all|ID,...          Collected target channels; track uses all unmasked
  --mask none|ID,...             Same excluded channels in data and MC/track
  --min-layers 7 --max-angle 0.05 --max-track-rms 3
  --long-margin 1 --trans-margin 0.25
  --u-bins 3 --v-bins 1         Odd bin counts; central disjoint reference
  --position-resolution 1.5     Assumed resolution floor mm; bin width >= 2x
  --temperature-bin 2 --slope-bin 0.02 --min-fit-hits 100
  --fit-chi2max 3 --max-relative-error 0.15 --max-threshold-shift 0.05
  --split-mode event|run --validation-runs 44
  --map-axis u|v --map-condition position_matched|unadjusted|data_mc_matched
  --workers 2                   Separate fit processes; one ROOT per worker
  --plot-channels all|ID,... --plot-limit 12 (0=all collected channels)
  --plot-label NAME              PNG generation directory (default: default)
  --refit-from PATH              Refit saved histograms in a fresh label
  --reuse-collection PATH        Completed label from THIS schema; no raw reread
  --collect-only --no-plots --plots-only --dry-run
  --result-dir PATH              Analysis result base, before label
  --help
plots-only reads saved fits/histograms/TSV. Use a NEW --plot-label to avoid overwrite.
HELP
}
declare -A opt=( [label]=pilot [runs]=44 [energy]=100 [sample]=both [mc-tag]=threshold
 [input-manifest]='' [data-base]=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps [mc-base]="$MC/Result_MC"
 [result-dir]="$PROJECT/analysis/result/channel_response/strip_position_response" [max-events]=1000 [channels]=all [mask]=none
 [min-layers]=7 [max-angle]=0.05 [max-track-rms]=3 [long-margin]=1 [trans-margin]=0.25
 [u-bins]=3 [v-bins]=1 [position-resolution]=1.5 [temperature-bin]=2 [slope-bin]=0.02
 [min-fit-hits]=100 [fit-chi2max]=3 [max-relative-error]=0.15 [max-threshold-shift]=0.05
 [split-mode]=event [validation-runs]=44 [map-axis]=u [map-condition]=position_matched
 [workers]=2 [plot-channels]=all [plot-limit]=12 [plot-label]=default [reuse-collection]='' [refit-from]=''
 [data-pedestal]="$PARAMS/pedestal2023_SPS.root" [mc-pedestal]="$PARAMS/pedestal2023_SPS.root"
 [data-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [mc-mip]="$PARAMS/all_auto_muon_v4_trackfit.root"
 [data-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root" [mc-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root"
 [data-threshold]="$PROJECT/analysis/result/threshold/threshold.root" [mc-threshold]="$PROJECT/analysis/result/threshold/threshold.root" )
declare -A given=()
dry=0; plots=0; collect=0; noplot=0
while (($#));do case $1 in
 --help|-h)usage;exit 0;; --dry-run)dry=1;shift;; --plots-only)plots=1;shift;; --collect-only)collect=1;shift;; --no-plots)noplot=1;shift;;
 --*)k=${1#--};[[ -v opt[$k] && $# -ge 2 ]]||fail "Unknown/missing $1";opt[$k]=$2;given[$k]=1;shift 2;; *)fail "Unexpected argument $1";;esac;done
for k in label mc-tag plot-label;do [[ ${opt[$k]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && ${opt[$k]} != .. ]]||fail "Invalid $k";done
out=$(realpath -m -- "${opt[result-dir]}")/${opt[label]}
if ((plots));then
 ((collect==0 && noplot==0))||fail 'Incompatible plotting flags'
 for k in "${!given[@]}";do case $k in label|result-dir|plot-label|plot-channels|plot-limit|workers);;*)fail "plots-only cannot change --$k";;esac;done
 [[ -f $out/COMPLETE && -f $out/fits.root && -f $out/histograms.root && -f $out/validation.root ]]||fail 'Missing completed fit outputs'
 [[ ! -e $out/figures/${opt[plot-label]} ]]||fail 'Choose a new --plot-label; existing figures are preserved'
 printf 'Plots-only: %s -> figures/%s (no original input reads)\n' "$out" "${opt[plot-label]}"
 ((dry==0))||exit 0
 exec 8>"$out/.run.lock";flock -n 8||fail 'Output in use'
 mkdir -p "$out/figures/${opt[plot-label]}"
 plotdir=$out/figures/${opt[plot-label]}
 tar -czf "$plotdir/sources.tar.gz" -C "$PROJECT" analysis/channel_response/strip_position_response/{Draw.cc,Common.hh,run.sh,Makefile} analysis/channel_response/adc_energy/FastPng.hh analysis/channel_response/root.mk
 cp "$out/inputs/configuration.tsv" "$plotdir/configuration.tsv"
 for k in plot-channels plot-limit;do if [[ -v given[$k] ]];then sed -i "/^$k"$'\t'"/d" "$plotdir/configuration.tsv";printf '%s\t%s\n' "$k" "${opt[$k]}" >> "$plotdir/configuration.tsv";fi;done
 (flock 9;make -C "$HERE" Draw) 9>"$HERE/.build.lock" > "$plotdir/build.log" 2>&1
 sha256sum "$out"/{histograms.root,fits.root,validation.root,bins.tsv,fits.tsv,comparison.tsv,groups.tsv,candidate_map.tsv,validation.tsv} > "$plotdir/checksums_before.txt"
 printf 'cd %q\n' "$PWD" > "$plotdir/command.sh"
 printf '%q ' bash "$HERE/run.sh" "${original[@]}" >> "$plotdir/command.sh";printf '\n' >> "$plotdir/command.sh"
 "$HERE/Draw" "$out" "$plotdir/configuration.tsv" "$plotdir" > "$plotdir/run.log" 2>&1
 sha256sum "$out"/{histograms.root,fits.root,validation.root,bins.tsv,fits.tsv,comparison.tsv,groups.tsv,candidate_map.tsv,validation.tsv} > "$plotdir/checksums_after.txt"
 cmp "$plotdir/checksums_before.txt" "$plotdir/checksums_after.txt"||fail 'Numerical outputs changed during plotting'
 printf 'Completed redraw: %s\n' "$plotdir";exit 0
fi
source_path=''; refit=0
[[ -z ${opt[reuse-collection]} || -z ${opt[refit-from]} ]]||fail 'Choose reuse-collection or refit-from'
if [[ -n ${opt[refit-from]} ]];then opt[reuse-collection]=${opt[refit-from]};refit=1;fi
if [[ -n ${opt[reuse-collection]} ]];then
 source_path=$(realpath -e -- "${opt[reuse-collection]}")
 [[ -f $source_path/COLLECTED && -f $source_path/inputs/configuration.tsv ]]||fail 'Reuse requires completed collection and saved configuration'
 while IFS=$'\t' read -r k v;do case $k in
 runs|sample|energy|mc-tag|max-events|channels|mask|min-layers|max-angle|max-track-rms|long-margin|trans-margin|input-manifest|data-base|mc-base|data-pedestal|mc-pedestal|data-mip|mc-mip|data-hl|mc-hl|data-threshold|mc-threshold)
 [[ ! -v given[$k] || ${opt[$k]} == "$v" ]]||fail "Incompatible collection option --$k";opt[$k]=$v;;esac;done < "$source_path/inputs/configuration.tsv"
 [[ -f $source_path/inputs/collection_code.sha256 && -f $source_path/inputs/collection_products.sha256 ]]||fail 'Missing collection provenance'
 sha256sum -c "$source_path/inputs/collection_products.sha256" >/dev/null||fail 'Reused collection/configuration changed'
 if ((refit));then
  [[ -f $source_path/COMPLETE && -f $source_path/histograms.root ]]||fail 'Refit requires completed histograms'
  [[ ! -f $source_path/inputs/analysis_products.sha256 ]]||sha256sum -c "$source_path/inputs/analysis_products.sha256" >/dev/null||fail 'Saved histogram/selection provenance changed'
  while IFS=$'\t' read -r k v;do case $k in u-bins|v-bins|position-resolution|temperature-bin|slope-bin|split-mode|validation-runs)
   [[ ! -v given[$k] || ${opt[$k]} == "$v" ]]||fail "Refit cannot change --$k";opt[$k]=$v;;esac;done < "$source_path/inputs/configuration.tsv"
 else sha256sum -c "$source_path/inputs/collection_code.sha256" >/dev/null||fail 'Collection definition changed; recollect inputs';fi
fi
for k in max-events min-layers u-bins v-bins min-fit-hits workers plot-limit;do [[ ${opt[$k]} =~ ^(0|[1-9][0-9]*)$ ]]||fail "Invalid integer $k";done
((opt[workers]>=1 && opt[workers]<=32 && opt[min-fit-hits]>=20 && opt[min-layers]>=3 && opt[min-layers]<=14 && opt[u-bins]>0 && opt[u-bins]%2==1 && opt[v-bins]>0 && opt[v-bins]%2==1))||fail 'Invalid counts'
[[ -z ${SLURM_CPUS_PER_TASK:-} || ${opt[workers]} -le $SLURM_CPUS_PER_TASK ]]||fail 'Workers exceed allocated CPUs'
for k in energy max-angle max-track-rms long-margin trans-margin position-resolution temperature-bin slope-bin fit-chi2max max-relative-error max-threshold-shift;do [[ ${opt[$k]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k (use 0.05 form)";done
awk -v e="${opt[energy]}" -v a="${opt[max-angle]}" -v r="${opt[max-track-rms]}" -v l="${opt[long-margin]}" -v t="${opt[trans-margin]}" -v pr="${opt[position-resolution]}" -v u="${opt[u-bins]}" -v v="${opt[v-bins]}" -v tb="${opt[temperature-bin]}" -v sb="${opt[slope-bin]}" 'BEGIN{exit !(e>0&&a>0&&a<0.5&&r>0&&l<22.5&&t<2.5&&pr>0&&45/u>=2*pr&&5/v>=2*pr&&tb>0&&sb>0)}'||fail 'Invalid geometry/bin settings'
[[ ${opt[sample]} =~ ^(both|data|mc)$ && ${opt[split-mode]} =~ ^(event|run)$ && ${opt[map-axis]} =~ ^(u|v)$ && ${opt[map-condition]} =~ ^(unadjusted|position_matched|data_mc_matched)$ ]]||fail 'Invalid sample/split/map setting'
for k in runs channels mask plot-channels validation-runs;do [[ ${opt[$k]} == all || ${opt[$k]} == none || ${opt[$k]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail "Invalid $k";done
for s in data mc;do for p in pedestal mip hl threshold;do [[ -f ${opt[$s-$p]} ]]||fail "Missing $s-$p: ${opt[$s-$p]}";done;done
rows=()
if [[ -z $source_path ]];then
 if [[ -n ${opt[input-manifest]} ]];then
  while IFS= read -r line;do [[ -z $line || $line == \#* || $line == sample$'\t'* ]]||rows+=("$line");done < "${opt[input-manifest]}"
 else
  shopt -s nullglob
  if [[ ${opt[sample]} == mc ]];then files=("${opt[mc-base]}/decode/mu-/mu_track/${opt[mc-tag]}/${opt[energy]}GeV/"ECAL_Run*.root);else files=("${opt[data-base]}/decode/mu-/${opt[energy]}GeV/"ECAL_Run*.root);fi
  shopt -u nullglob
  for path in "${files[@]}";do name=${path##*/};[[ $name =~ ^ECAL_Run([0-9]+) ]]||continue;run=${BASH_REMATCH[1]};[[ ${opt[runs]} == all || ,${opt[runs]}, == *,$run,* ]]||continue
   if [[ ${opt[sample]} != mc ]];then rows+=("0"$'\t'"$run"$'\t'"$path"$'\t'"${opt[data-base]}/calib/mu-/${opt[energy]}GeV/$name"$'\t-\t-');fi
   if [[ ${opt[sample]} != data ]];then stub=mu-/mu_track/${opt[mc-tag]}/${opt[energy]}GeV/ECAL_Run${run}.root;rows+=("1"$'\t'"$run"$'\t'"${opt[mc-base]}/decode/$stub"$'\t'"${opt[mc-base]}/calib/$stub"$'\t-'$'\t'"${opt[mc-base]}/generate/mu-/mu_track/${opt[energy]}GeV/ECAL_Run${run}.root");fi
  done
 fi
 declare -A seen=() run_seen=() paths=();selected_rows=();data_n=0;mc_n=0
 for row in "${rows[@]}";do IFS=$'\t' read -r s r raw cal unused truth extra <<< "$row";[[ -z ${extra:-} && $s =~ ^[01]$ && $r =~ ^[0-9]+$ && -n $truth ]]||fail 'Malformed manifest'
  [[ ${opt[sample]} != data || $s == 0 ]]||continue;[[ ${opt[sample]} != mc || $s == 1 ]]||continue;[[ ${opt[runs]} == all || ,${opt[runs]}, == *,$r,* ]]||continue
  [[ ! -v seen[$s,$r] ]]||fail 'Repeated sample/run';seen[$s,$r]=1;run_seen[$r]=1
  for path in "$raw" "$cal" "$truth";do [[ $path == - ]]&&continue;[[ $path == /* && -f $path ]]||fail "Missing/relative input: $path";canonical=$(realpath -e -- "$path");[[ ! -v paths[$canonical] ]]||fail 'Duplicate input file';paths[$canonical]=1;done
  [[ $s == 0 || $truth != - ]]||fail 'MC requires actual muon truth';selected_rows+=("$row");if [[ $s == 0 ]];then data_n=$((data_n+1));else mc_n=$((mc_n+1));fi
 done
 rows=("${selected_rows[@]}");((${#rows[@]}))||fail 'No real muon input at requested energy/runs'
 [[ ${opt[sample]} == mc || $data_n -gt 0 ]]||fail 'No data inputs';[[ ${opt[sample]} == data || $mc_n -gt 0 ]]||fail 'No MC inputs'
 if [[ ${opt[runs]} != all ]];then IFS=, read -r -a want <<< "${opt[runs]}";for r in "${want[@]}";do [[ -v run_seen[$r] ]]||fail "Run $r absent";done;fi
else printf 'Reuse collection: %s\n' "$source_path";fi
[[ ! -e $out ]]||fail "Result directory exists; choose a new label: $out"
printf 'Output: %s\nenergy=%s GeV; max-events/file=%s; workers=%s; u/v bins=%s/%s; no CoG, tail or energy cuts\n' "$out" "${opt[energy]}" "${opt[max-events]}" "${opt[workers]}" "${opt[u-bins]}" "${opt[v-bins]}"
if [[ -z $source_path ]];then printf '%s\n' "${rows[@]}";fi
((dry==0))||exit 0
mkdir -p "$out/inputs"
exec 8>"$out/.run.lock";flock -n 8||fail 'Output in use'
printf 'cd %q\n' "$PWD" > "$out/inputs/command.sh"
printf '%q ' bash "$HERE/run.sh" "${original[@]}" >> "$out/inputs/command.sh";printf '\n' >> "$out/inputs/command.sh"
{
 date -Is
 (flock 9;make -C "$HERE" -j"${opt[workers]}" all) 9>"$HERE/.build.lock"
 mkdir -p "$out/inputs/constants" "$out/inputs/current_sources"
 declare -A copied=()
 for s in data mc;do for p in pedestal mip hl threshold;do k=$s-$p;src=$(realpath -e -- "${opt[$k]}");if [[ ! -v copied[$src] ]];then dest="$out/inputs/constants/${src##*/}";[[ ! -e $dest ]]||dest="$out/inputs/constants/$k.root";cp -- "$src" "$dest";copied[$src]=$dest;fi;opt[$k]=${copied[$src]};done;done
 for k in "${!opt[@]}";do printf '%s\t%s\n' "$k" "${opt[$k]}";done | LC_ALL=C sort > "$out/inputs/configuration.tsv"
 cfg=$out/inputs/configuration.tsv
 sha256sum "$out"/inputs/constants/* > "$out/inputs/constant_checksums.txt"
 sha256sum "$HERE/Collect.cc" "$HERE/Common.hh" "$HERE/../muon_mip_response/"{MuonMIPCommon.hh,TrackGeometry.hh,CurrentCalibration.hh} /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/EBUdecode.cxx > "$out/inputs/collection_code.sha256"
 tar -czf "$out/inputs/analysis_sources.tar.gz" -C "$PROJECT" analysis/channel_response/strip_position_response/{Common.hh,Collect.cc,Aggregate.cc,Fit.cc,Draw.cc,Makefile,run.sh,submit.sh,README.md} analysis/channel_response/muon_mip_response/{MuonMIPCommon.hh,CurrentCalibration.hh,TrackGeometry.hh} analysis/channel_response/{root.mk,execute.sl} analysis/channel_response/adc_energy/FastPng.hh analysis/mip_calibration/{Common.hh,FitModel.hh}
 for file in "$MC/Analysis_edit/src/"{MCDigi.cxx,Extract.cxx} "$MC/Analysis_edit/"{digitization_sps.sh,digitization_e.sh} /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/{Calibration.cxx,EBUdecode.cxx};do [[ ! -f $file ]]||cp "$file" "$out/inputs/current_sources/";done
 sha256sum "$out"/inputs/current_sources/* > "$out/inputs/current_source_checksums.txt"
 if [[ -z $source_path ]];then
  printf '%s\n' "${rows[@]}" > "$out/inputs/manifest.tsv"
  input_stats(){ while IFS=$'\t' read -r s r raw cal unused truth;do for f in "$raw" "$cal" "$truth";do [[ $f == - ]]||stat -L --printf='%n\t%s\t%Y\n' -- "$f";done;done < "$out/inputs/manifest.tsv"; }
  input_stats > "$out/inputs/input_stats_before.tsv"
  "$HERE/Collect" "$out/inputs/manifest.tsv" "$cfg" "$out"
  input_stats > "$out/inputs/input_stats_after.tsv"
  cmp "$out/inputs/input_stats_before.tsv" "$out/inputs/input_stats_after.tsv"||fail 'Inputs changed during collection'
 else
  for f in crossings.root parameters.tsv event_audit.tsv production_audit.tsv;do ln -s "$source_path/$f" "$out/$f";done
  cp "$source_path/inputs/manifest.tsv" "$out/inputs/manifest.tsv"
  cp "$source_path/inputs/collection_code.sha256" "$out/inputs/source_collection_code.sha256"
  sha256sum "$source_path/crossings.root" > "$out/inputs/reused_collection.sha256"
 fi
 sha256sum "$out"/{crossings.root,parameters.tsv,event_audit.tsv,production_audit.tsv} "$out"/inputs/{configuration.tsv,manifest.tsv} > "$out/inputs/collection_products.sha256"
 touch "$out/COLLECTED"
 if ((collect==0));then
  if ((refit));then
   for f in histograms.root bins.tsv matching.tsv channel_coverage.tsv event_summary.tsv;do ln -s "$source_path/$f" "$out/$f";done
   sha256sum "$source_path/histograms.root" "$source_path/bins.tsv" "$source_path/matching.tsv" > "$out/inputs/refit_histograms.sha256"
  else "$HERE/Aggregate" "$out/crossings.root" "$cfg" "$out";fi
  mkdir -p "$PROJECT/.scratch/strip_position_response"
  scratch=$(mktemp -d "$PROJECT/.scratch/strip_position_response/fit.XXXXXX")
  printf 'Fit shards: %s\n' "$scratch"
  pids=();shards=();for ((i=0;i<opt[workers];++i));do shard=$scratch/worker$i.root;shards+=("$shard");"$HERE/Fit" worker "$out" "$cfg" "$i" "${opt[workers]}" "$shard" > "$out/fit_worker$i.log" 2>&1 & pids+=("$!");done
  failed=0;for pid in "${pids[@]}";do wait "$pid"||failed=1;done;((failed==0))||fail "Fit failed; retained shards: $scratch"
  hadd -f "$out/fits.root" "${shards[@]}"
  "$HERE/Fit" summarize "$out" "$cfg"
  rm -- "${shards[@]}";rmdir "$scratch"
  if ((noplot==0));then "$HERE/Draw" "$out" "$cfg" "$out/figures/${opt[plot-label]}";fi
  if ((refit));then sha256sum -c "$out/inputs/refit_histograms.sha256";fi
  sha256sum "$out"/{histograms.root,bins.tsv,matching.tsv} "$out/inputs/configuration.tsv" > "$out/inputs/analysis_products.sha256"
  touch "$out/COMPLETE"
 fi
 [[ -z $source_path ]]||sha256sum -c "$out/inputs/reused_collection.sha256"
 date -Is
} 2>&1 | tee "$out/run.log"
printf 'Completed: %s\n' "$out"
