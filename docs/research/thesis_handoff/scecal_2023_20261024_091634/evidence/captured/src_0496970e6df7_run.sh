#!/bin/bash
set -euo pipefail
original_args=("$@")
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../../.." && pwd -P)
MC=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
PARAMS=$MC/Analysis_edit/share
fail(){ echo "ERROR: $*" >&2; exit 1; }
usage(){ cat <<HELP
Usage: bash run.sh [options]
  --runs all|44|29,44       Muon run numbers (default: all)
  --sample both|data|mc     Same analysis for either/both samples (default: both)
  --energy 100              Muon beam energy, GeV (default: 100)
  --mc-tag threshold        MC digitization directory
  --label NAME              Output label (default: default)
  --result-dir PATH         Default: analysis/result/channel_response/muon_mip_response
  --data-base PATH          Default: ECAL_data/analysed/2023/sps
  --mc-base PATH            Default: Result_MC
  --input-manifest PATH     Optional six columns: sample,run,decode,calib,-,truth
  --data-mip PATH --mc-mip PATH --digi-mip PATH
  --data-pedestal PATH --mc-pedestal PATH
  --data-hl PATH --mc-hl PATH
  --data-threshold PATH --mc-threshold PATH
  --max-events N            Evenly spaced events per file (0: all; default: 0)
  --min-layers N            Other layers required per projection (default: 7)
  --max-angle RAD           Projected angle limit (default: 0.05)
  --max-track-rms MM        Track residual limit (default: 3)
  --long-margin MM --trans-margin MM   Strip edge margins (1, 0.25)
  --long-bin MM --transverse-bin MM --slope-bin VALUE (5, 1.5, 0.01)
  --min-bin N               Minimum data/MC hits in a common geometry bin (3)
  --min-fit-hits N          Effective entries in fit interval (100)
  --plateau-sigmas N        Fit above threshold + N sigma (3)
  --fit-chi2max N           Maximum fit chi2/NDF (3)
  --channels all|CELLID,... Channels to fit/draw (default: all)
  --dependency-fits 0|1     Fit MPV vs position/temperature/run (default: 1)
  --workers N              Independent channel-fit processes (default: 4)
  --fit-cache 0|1          Reuse exactly identical fits (default: 1)
  --fit-cache-mb N         Cache budget per worker (default: 32 MiB)
  --reuse-flat PATH        Analyze completed collected ROOT without recollection
  --allow-recovered         Include recovered ROOT files, diagnostic only
  --collect-only            Input audit + collection; omit fits/figures
  --plots-only              Redraw all existing muon PNGs; no fits/collection
  --dry-run                 Check paths and show configuration; no writes/build/job
  --help

No saved trackFit3D is required. Inputs are read-only. PNG only; ROOT/TSV retained.
HELP
}
declare -A opt=(
 [runs]=all [sample]=both [energy]=100 [mc-tag]=threshold [label]=default
 [result-dir]="$HERE/../../result/channel_response/muon_mip_response"
 [data-base]=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps [mc-base]="$MC/Result_MC" [input-manifest]=''
 [data-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [mc-mip]="$PARAMS/all_auto_muon_v4_trackfit.root" [digi-mip]="$PARAMS/all_auto_muon_v4_trackfit.root"
 [data-pedestal]="$PARAMS/pedestal2023_SPS.root" [mc-pedestal]="$PARAMS/pedestal2023_SPS.root"
 [data-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root" [mc-hl]="$PARAMS/all_hl_electron2023_hlratio_v3.root"
 [data-threshold]="$PROJECT/analysis/result/threshold/threshold.root" [mc-threshold]="$PROJECT/analysis/result/threshold/threshold.root"
 [max-events]=0 [min-layers]=7 [max-angle]=0.05 [max-track-rms]=3 [long-margin]=1 [trans-margin]=0.25
 [long-bin]=5 [transverse-bin]=1.5 [slope-bin]=0.01 [min-bin]=3 [min-fit-hits]=100 [plateau-sigmas]=3 [fit-chi2max]=3
 [channels]=all [dependency-fits]=1 [workers]=4 [fit-cache]=1 [fit-cache-mb]=32 [reuse-flat]='' [allow-recovered]=0 [collect-only]=0 [plots-only]=0 )
dry=0
declare -A given=()
while (($#));do case $1 in
 --help|-h)usage;exit 0;; --dry-run)dry=1;shift;;
 --plots-only)opt[plots-only]=1;shift;;
 --allow-recovered)opt[allow-recovered]=1;given[allow-recovered]=1;shift;; --collect-only)opt[collect-only]=1;shift;;
 --*) key=${1#--};[[ -v opt[$key] && $# -ge 2 ]]||fail "Unknown/missing option $1";opt[$key]=$2;given[$key]=1;shift 2;;
 *)fail "Unexpected argument $1";; esac;done
for k in label mc-tag;do [[ ${opt[$k]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && ${opt[$k]} != .. ]]||fail "Invalid $k";done
case ${opt[sample]} in both|data|mc);;*)fail 'sample must be both, data or mc';;esac
[[ ${opt[runs]} == all || ${opt[runs]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid runs'
[[ ${opt[channels]} == all || ${opt[channels]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid channels'
for k in max-events min-layers min-bin min-fit-hits dependency-fits workers fit-cache fit-cache-mb allow-recovered collect-only plots-only;do [[ ${opt[$k]} =~ ^(0|[1-9][0-9]*)$ ]]||fail "Invalid integer $k";done
((opt[workers]>=1 && opt[workers]<=64 && opt[fit-cache]<=1 && opt[fit-cache-mb]<=4096 && opt[min-layers]>=3 && opt[min-layers]<=14 && opt[min-bin]>=1 && opt[min-fit-hits]>=20 && opt[dependency-fits]<=1 && opt[allow-recovered]<=1 && opt[collect-only]<=1 && opt[plots-only]<=1))||fail 'Invalid count/cut value'
for k in energy max-angle max-track-rms long-margin trans-margin long-bin transverse-bin slope-bin plateau-sigmas fit-chi2max;do
 [[ ${opt[$k]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k"
done
awk -v e="${opt[energy]}" -v a="${opt[max-angle]}" -v r="${opt[max-track-rms]}" -v l="${opt[long-margin]}" -v t="${opt[trans-margin]}" -v b="${opt[long-bin]}" -v c="${opt[transverse-bin]}" -v s="${opt[slope-bin]}" -v q="${opt[fit-chi2max]}" 'BEGIN{exit !(e>0&&a>0&&a<0.5&&r>0&&l<22.5&&t<2.5&&b>0&&c>0&&s>0&&q>0)}'||fail 'Invalid cut/bin value'
if [[ -n ${SLURM_CPUS_PER_TASK:-} && ${opt[workers]} -gt ${SLURM_CPUS_PER_TASK} ]];then fail "workers=${opt[workers]} exceeds allocated CPUs=${SLURM_CPUS_PER_TASK}";fi
if ((opt[plots-only]));then
 [[ -z ${opt[reuse-flat]} && ${opt[collect-only]} == 0 ]]||fail 'plots-only cannot be combined with collection/reuse-flat'
 for key in "${!given[@]}";do
  case $key in workers|channels|label|result-dir);;*)fail "--$key would change analysis settings; plots-only uses saved fits";;esac
 done
 out=$(realpath -m -- "${opt[result-dir]}")/${opt[label]}
 for input_name in spectra.root fits.tsv inputs/configuration.tsv;do [[ -f $out/$input_name ]]||fail "Missing saved plot input: $out/$input_name";done
 printf 'Plots only: %s\nworkers=%s; channels=%s; existing response, MPV/mean profiles, and channel map PNGs\n' "$out" "${opt[workers]}" "${opt[channels]}"
 (( !dry ))||exit 0
 (flock 9;make -C "$HERE" all) 9>"$HERE/.build.lock"
 exec 8>"$out/.run.lock";flock -n 8||fail 'Another process is writing this output'
 plot_inputs=("$out/spectra.root" "$out/fits.tsv")
 [[ ! -f $out/dependency_fits.tsv ]]||plot_inputs+=("$out/dependency_fits.tsv")
 sha256sum "${plot_inputs[@]}" > "$out/inputs/plot_input_checksums_before.txt"
 sha256sum "$HERE/RedrawMuonMIP.cc" "$HERE/MuonMIPPlotStyle.hh" "$HERE/MuonMIPCommon.hh" "$HERE/run.sh" "$HERE/Makefile" "$HERE/../adc_energy/FastPng.hh" > "$out/inputs/plot_code_checksums.txt"
 printf 'cd %q\n' "$PWD" > "$out/inputs/plot_command.sh"
 printf '%q ' "$HERE/run.sh" "${original_args[@]}" >> "$out/inputs/plot_command.sh"
 printf '\n' >> "$out/inputs/plot_command.sh"
 {
  date -Is
  "$HERE/RedrawMuonMIP" "$out" --workers "${opt[workers]}" --channels "${opt[channels]}"
  sha256sum "${plot_inputs[@]}" > "$out/inputs/plot_input_checksums_after.txt"
  cmp -s "$out/inputs/plot_input_checksums_before.txt" "$out/inputs/plot_input_checksums_after.txt"||fail 'Saved numerical plot input changed during redraw'
  date -Is
  printf 'Plots updated from saved histograms and fits; input checksums unchanged.\n'
 } 2>&1 | tee -a "$out/run.log"
 exit 0
fi
reuse_source=''
if [[ -n ${opt[reuse-flat]} ]];then
 [[ -f ${opt[reuse-flat]} ]]||fail 'Missing --reuse-flat input'
 opt[reuse-flat]=$(realpath -e -- "${opt[reuse-flat]}")
 reuse_source=$(dirname -- "${opt[reuse-flat]}")
 [[ -f $reuse_source/inputs/configuration.tsv && ( -f $reuse_source/inputs/input_stats_after.tsv || -f $reuse_source/COMPLETE ) ]]||fail 'reuse-flat requires completed collection and its saved inputs/configuration.tsv'
 if [[ -f $reuse_source/inputs/input_stats_after.tsv ]];then
  [[ -f $reuse_source/inputs/input_stats_before.tsv ]]&&cmp -s -- "$reuse_source/inputs/input_stats_before.tsv" "$reuse_source/inputs/input_stats_after.tsv"||fail 'Reused collection failed the input-stability check'
 fi
 ((opt[collect-only]==0))||fail 'reuse-flat and collect-only cannot be combined'
 while IFS=$'\t' read -r key stored;do
  case $key in sample|runs|energy|mc-tag|max-events|min-layers|max-angle|max-track-rms|long-margin|trans-margin|allow-recovered|data-base|mc-base|input-manifest)
   if [[ -v given[$key] && ${opt[$key]} != "$stored" ]];then fail "Collection option --$key differs from reused flat ROOT";fi
   opt[$key]=$stored;;
  data-mip|mc-mip|digi-mip|data-pedestal|mc-pedestal|data-hl|mc-hl|data-threshold|mc-threshold)
   [[ -v given[$key] ]]||opt[$key]=$stored;;
  esac
 done < "$reuse_source/inputs/configuration.tsv"
fi
for k in data-mip mc-mip digi-mip data-pedestal mc-pedestal data-hl mc-hl data-threshold mc-threshold;do
 [[ -f ${opt[$k]} ]]||fail "Missing $k: ${opt[$k]}";opt[$k]=$(realpath -e -- "${opt[$k]}")
done
rows=()
if [[ -z ${opt[reuse-flat]} ]];then
if [[ -n ${opt[input-manifest]} ]];then
 [[ -f ${opt[input-manifest]} ]]||fail 'Missing input manifest'
 while IFS= read -r row;do [[ -z $row || $row == \#* ]]||rows+=("$row");done < "${opt[input-manifest]}"
else
 shopt -s nullglob
 if [[ ${opt[sample]} == mc ]];then
  files=("${opt[mc-base]}/decode/mu-/mu_track/${opt[mc-tag]}/${opt[energy]}GeV/"ECAL_Run*.root)
 else
  files=("${opt[data-base]}/decode/mu-/${opt[energy]}GeV/"ECAL_Run*.root)
 fi
 shopt -u nullglob
 ((${#files[@]}))||fail 'No muon runs; use --data-base, --mc-base or --input-manifest'
 declare -A found=()
 for raw in "${files[@]}";do
  name=${raw##*/};[[ $name =~ ^ECAL_Run([0-9]+) ]]||continue;run=${BASH_REMATCH[1]}
  [[ ${opt[runs]} == all || ,${opt[runs]}, == *,$run,* ]]||continue
  [[ ! -v found[$run] ]]||fail "Multiple files for run $run; use --input-manifest";found[$run]=1
  if [[ ${opt[sample]} != mc ]];then rows+=("0"$'\t'"$run"$'\t'"$raw"$'\t'"${opt[data-base]}/calib/mu-/${opt[energy]}GeV/$name"$'\t-'$'\t-');fi
  if [[ ${opt[sample]} != data ]];then
   mpath=mu-/mu_track/${opt[mc-tag]}/${opt[energy]}GeV/ECAL_Run${run}.root
   rows+=("1"$'\t'"$run"$'\t'"${opt[mc-base]}/decode/$mpath"$'\t'"${opt[mc-base]}/calib/$mpath"$'\t-'$'\t'"${opt[mc-base]}/generate/mu-/mu_track/${opt[energy]}GeV/ECAL_Run${run}.root")
  fi
 done
 if [[ ${opt[runs]} != all ]];then IFS=, read -r -a wanted <<< "${opt[runs]}";for run in "${wanted[@]}";do [[ -v found[$run] ]]||fail "Muon run $run not found";done;fi
fi
((${#rows[@]}))||fail 'No requested inputs'
filtered=()
declare -A seen_runs=() path_ok=()
for row in "${rows[@]}";do
 IFS=$'\t' read -r sample run raw cal unused truth extra <<< "$row"
 [[ -z ${extra:-} && ( $sample == 0 || $sample == 1 ) && $run =~ ^[0-9]+$ && -n $raw && -n $cal && -n $truth ]]||fail 'Manifest requires six tab-separated columns'
 [[ ${opt[sample]} != data || $sample == 0 ]]||continue
 [[ ${opt[sample]} != mc || $sample == 1 ]]||continue
 [[ ${opt[runs]} == all || ,${opt[runs]}, == *,$run,* ]]||continue
 for path in "$raw" "$cal" "$truth";do [[ $path == - || $path == /* ]]||fail 'Manifest input paths must be absolute';done
 filtered+=("$row");seen_runs[$run]=1
 [[ -v path_ok[$run] ]]||path_ok[$run]=1
 for path in "$raw" "$cal" "$truth";do [[ $path == - ]]&&continue;[[ -f $path ]]||{ echo "Unavailable input (run $run): $path" >&2;path_ok[$run]=0; };done
 printf '%s run %s: %s\n' "$([[ $sample == 0 ]]&&echo data||echo mc)" "$run" "$raw"
done
rows=("${filtered[@]}")
((${#rows[@]}))||fail 'No manifest rows match sample/runs'
if [[ ${opt[runs]} != all ]];then IFS=, read -r -a wanted <<< "${opt[runs]}";for run in "${wanted[@]}";do [[ -v seen_runs[$run] ]]||fail "Requested run $run absent from selected manifest";done;fi
usable_data=0
usable_mc=0
for row in "${rows[@]}";do IFS=$'\t' read -r sample run unused <<< "$row";if ((path_ok[$run]));then if [[ $sample == 0 ]];then usable_data=$((usable_data+1));else usable_mc=$((usable_mc+1));fi;fi;done
[[ ${opt[sample]} == mc || $usable_data -gt 0 ]]||fail 'No complete data input run remains'
[[ ${opt[sample]} == data || $usable_mc -gt 0 ]]||fail 'No complete MC input run remains'
else
 printf 'Reusing collected input: %s\n' "${opt[reuse-flat]}"
fi
out=$(realpath -m -- "${opt[result-dir]}")/${opt[label]}
printf 'Output: %s\nmax-events/file=%s; leave-layer-out min-layers/projection=%s; PNG only; workers=%s; fit-cache=%s\n' "$out" "${opt[max-events]}" "${opt[min-layers]}" "${opt[workers]}" "${opt[fit-cache]}"
printf 'Run pairs with missing/recovered input are excluded by the ROOT metadata audit.\n'
((!dry))||exit 0
[[ ! -e $out/flat.root && ! -e $out/COMPLETE ]]||fail "Output exists; choose another --label: $out"
( flock 9; make -C "$HERE" all ) 9>"$HERE/.build.lock"
mkdir -p -- "$out/inputs"
exec 8>"$out/.run.lock";flock -n 8||fail 'Another process is writing this output'
[[ ! -e $out/flat.root && ! -e $out/COMPLETE ]]||fail 'Output appeared while waiting for lock'
if [[ -z ${opt[reuse-flat]} ]];then
 printf '%s\n' "${rows[@]}" > "$out/inputs/requested.tsv"
else
 mkdir -p -- "$out/inputs/collection"
 cp -- "$reuse_source/inputs/configuration.tsv" "$out/inputs/collection/configuration.tsv"
 for name in manifest.tsv requested.tsv file_metadata.tsv parameter_checksums.txt code_checksums.txt input_stats_before.tsv input_stats_after.tsv;do [[ ! -f $reuse_source/inputs/$name ]]||cp -- "$reuse_source/inputs/$name" "$out/inputs/collection/$name";done
 [[ ! -f $reuse_source/validation.tsv ]]||cp -- "$reuse_source/validation.tsv" "$out/validation.tsv"
 ln -s -- "${opt[reuse-flat]}" "$out/flat.root"
 stat -L --printf='%n\t%s\t%Y\n' -- "${opt[reuse-flat]}" > "$out/inputs/flat_stats_before.tsv"
fi
for k in "${!opt[@]}";do [[ ${opt[$k]} != *$'\t'* && ${opt[$k]} != *$'\n'* ]]||fail 'Tabs/newlines in configuration';printf '%s\t%s\n' "$k" "${opt[$k]}";done | sort > "$out/inputs/configuration.tsv"
sha256sum "$HERE/"*.cc "$HERE/"*.hh "$HERE/"*.sh "$HERE/Makefile" "$HERE/../root.mk" "$HERE/../execute.sl" "$PROJECT/analysis/calibration_parameters/EBUdecode.cxx" "$HERE/../adc_energy/FastPng.hh" > "$out/inputs/code_checksums.txt"
for k in data-mip mc-mip digi-mip data-pedestal mc-pedestal data-hl mc-hl data-threshold mc-threshold;do sha256sum "${opt[$k]}";done > "$out/inputs/parameter_checksums.txt"
for source in "$MC/Analysis_edit/src/MCDigi.cxx" "$MC/Analysis_edit/src/Extract.cxx" /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx /home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/MIPFit.cxx;do
 [[ -f $source ]]||continue;cp -- "$source" "$out/inputs/current_${source##*/}";sha256sum "$source"
done > "$out/inputs/current_source_checksums.txt"
printf 'cd %q\n' "$PWD" > "$out/inputs/command.sh"
printf '%q ' "$HERE/run.sh" "${original_args[@]}" >> "$out/inputs/command.sh"
printf '\n' >> "$out/inputs/command.sh"
{
 date -Is
 if [[ -z ${opt[reuse-flat]} ]];then
 "$HERE/InspectMuonInputs" "$out/inputs/requested.tsv" "$out/inputs/configuration.tsv" "$out/inputs/manifest.tsv" "$out/inputs/file_metadata.tsv"
 stat_inputs(){
  while IFS=$'\t' read -r sample_id run_id raw_path calib_path unused truth_path;do
   for input_path in "$raw_path" "$calib_path" "$truth_path";do [[ $input_path == - ]]&&continue;stat -L --printf='%n\t%s\t%Y\n' -- "$input_path";done
  done < "$out/inputs/manifest.tsv"
 }
 stat_inputs > "$out/inputs/input_stats_before.tsv"
 "$HERE/CollectMuonMIP" "$out/inputs/manifest.tsv" "$out/inputs/configuration.tsv" "$out/flat.root" "$out/validation.tsv"
 stat_inputs > "$out/inputs/input_stats_after.tsv"
 cmp -s "$out/inputs/input_stats_before.tsv" "$out/inputs/input_stats_after.tsv"||fail 'An input changed during collection; output is not accepted'

 fi
 if ((opt[collect-only]==0));then "$HERE/AnalyzeMuonMIP" "$out/flat.root" "$out/inputs/configuration.tsv" "$out";fi
 if [[ -n ${opt[reuse-flat]} ]];then
  stat -L --printf='%n\t%s\t%Y\n' -- "${opt[reuse-flat]}" > "$out/inputs/flat_stats_after.tsv"
  cmp -s "$out/inputs/flat_stats_before.tsv" "$out/inputs/flat_stats_after.tsv"||fail 'Reused flat input changed during analysis'
 fi
 date -Is
} 2>&1 | tee "$out/run.log"
touch "$out/COMPLETE"
echo "Completed: $out"
