#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../../.." && pwd -P)
args=("$@")
fail(){ echo "ERROR: $*" >&2; exit 1; }
usage(){ cat <<'HELP'
Usage: bash run.sh [options]
 --label NAME                  Fresh output label (default default)
 --energy GeV                  Actual muon energy; default 100, never electron fallback
 --sample both|data|mc          Default both
 --runs all|44|40,44            Default all; each sample/run must have one decode
 --input-manifest PATH          Six TSV columns: 0/1,run,decode,calib or -,-,truth or -
 --data-base PATH               Existing SPS analysed base
 --mc-base PATH                 Existing muon Result_MC base
 --mc-tag TAG                   Muon production configuration directory (threshold)
 --max-events N                Equally spaced events/file; 0=all (default)
 --min-layers N --max-angle X --max-track-rms X
 --long-margin X --trans-margin X --ambiguity-sigma X
 --stability 0|1               Nominal + loose/tight margins + tighter track (default 1)
 --groups '0,1,2,3;10,11,12,13;26,27,28,29'  Excluded layer groups, or none
 --mask PATH                   Optional TSV: sample(0/1),physical_CellID,reason
 --long-bin X --trans-bin X --slope-bin X --rms-bin X
 --min-common-bin N            Outcome-blind equal-count matching threshold (3)
 --min-count N --confidence X  Status threshold (20), interval confidence (0.95)
 --workers N                   Plot processes (default 2; bounded by allocated CPUs)
 --collect-only --no-plots     Stop after collection / omit drawing
 --reuse-crossings DIR         Reuse this analysis's completed collection; inherit settings
 --plots-only                  Read saved statistics only; PNG into a fresh redraw directory
 --result-dir PATH             Output base, label appended
 --dry-run                     Validate paths/options, no output/build/submission
 --build-only                  Compile only
HELP
}
declare -A o=( [label]=default [energy]=100 [sample]=both [runs]=all [mc-tag]=threshold
 [result-dir]="$PROJECT/analysis/result/channel_response/detection_efficiency"
 [data-base]=/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps
 [mc-base]=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC
 [input-manifest]='' [max-events]=0 [min-layers]=7 [max-angle]=0.05 [max-track-rms]=3
 [long-margin]=1 [trans-margin]=0.25 [ambiguity-sigma]=1 [stability]=1
 [groups]='0,1,2,3;10,11,12,13;26,27,28,29' [mask]=''
 [long-bin]=10 [trans-bin]=1 [slope-bin]=0.02 [rms-bin]=1 [min-common-bin]=3 [min-count]=20 [confidence]=0.95
 [workers]=2 [reuse-crossings]='' )
declare -A given=()
dry=0; plots=0; collect=0; noplot=0; build=0
while (($#));do case $1 in
 --help|-h)usage;exit 0;; --dry-run)dry=1;shift;; --plots-only)plots=1;shift;; --collect-only)collect=1;shift;; --no-plots)noplot=1;shift;; --build-only)build=1;shift;;
 --*)k=${1#--};[[ -v o[$k] && $# -ge 2 ]]||fail "Unknown/missing option $1";o[$k]=$2;given[$k]=1;shift 2;; *)fail "Unexpected $1";;esac;done
for k in label mc-tag;do [[ ${o[$k]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && ${o[$k]} != .. ]]||fail "Invalid $k";done
[[ ${o[workers]} =~ ^[1-9][0-9]*$ ]]||fail 'Invalid workers'
cpus=${SLURM_CPUS_PER_TASK:-$(nproc)};((o[workers]<=cpus))||fail "workers exceed allocated CPUs ($cpus)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
compile(){ (flock 9;make -C "$HERE" -j "$((${o[workers]}>1?2:1))") 9>"$HERE/.build.lock"; }
if ((build));then ((dry))||compile;exit 0;fi
out=$(realpath -m -- "${o[result-dir]}")/${o[label]}
if ((plots));then
 ((collect==0&&noplot==0))&&[[ -z ${o[reuse-crossings]} ]]||fail 'Incompatible plots-only mode'
 for k in "${!given[@]}";do case $k in label|result-dir|workers);;*)fail "plots-only cannot alter --$k";;esac;done
 [[ -f $out/COMPLETE && -f $out/statistics.root && -f $out/inputs/configuration.tsv ]]||fail 'Need completed saved statistics'
 printf 'Plots-only, read-only statistics: %s\nFresh redraw directory; workers=%s\n' "$out" "${o[workers]}"
 ((dry==0))||exit 0
 exec 8>"$out/.run.lock";flock -n 8||fail 'Output is in use'
 compile
 plotout=$(mktemp -d "$out/redraw_XXXXXXXX")
 sha256sum "$out/statistics.root" "$out"/*.tsv "$out/inputs/configuration.tsv" > "$plotout/numerical_before.sha256"
 printf '%q ' "$HERE/run.sh" "${args[@]}" > "$plotout/command.sh";printf '\n' >> "$plotout/command.sh"
 sha256sum "$HERE/Draw.cc" "$HERE/Common.hh" "$HERE/Draw" > "$plotout/plot_code.sha256"
 tar -cf "$plotout/plot_source.tar" -C "$PROJECT" analysis/channel_response/detection_efficiency/{Draw.cc,Common.hh,Makefile,run.sh} analysis/channel_response/muon_mip_response/{TrackGeometry.hh,MuonMIPCommon.hh} analysis/channel_response/adc_energy/FastPng.hh analysis/calibration_parameters/{EBUdecode.cxx,EBUdecode.h} analysis/channel_response/root.mk
 pids=();for ((i=0;i<o[workers];i++));do "$HERE/Draw" "$out/statistics.root" "$out/inputs/configuration.tsv" "$plotout" "$i" "${o[workers]}" > "$plotout/worker$i.log" 2>&1 & pids+=("$!");done
 failed=0;for pid in "${pids[@]}";do wait "$pid"||failed=1;done;((failed==0))||fail "Drawing failed; see $plotout"
 sha256sum "$out/statistics.root" "$out"/*.tsv "$out/inputs/configuration.tsv" > "$plotout/numerical_after.sha256"
 cmp "$plotout/numerical_before.sha256" "$plotout/numerical_after.sha256"||fail 'Numerical results changed during plots-only'
 printf 'Plots-only passed; numerical inputs unchanged. PNG: %s\n' "$plotout" | tee "$plotout/run.log"
 exit 0
fi
if [[ -e $out ]];then
 printf 'ERROR: Existing output is preserved: %s\n' "$out" >&2
 if [[ -f $out/COMPLETE && -f $out/statistics.root ]];then
  printf 'This label already has completed statistics; repeating collection is unnecessary.\nTo redraw the saved result:\n  bash ' >&2
  printf '%q ' "$HERE/run.sh" --result-dir "${o[result-dir]}" --label "${o[label]}" --plots-only >&2
  printf '\n' >&2
 elif [[ -f $out/COLLECTED ]];then
  printf 'This label already has a completed collection. Analyze it with an unused --label and --reuse-crossings %q.\n' "$out" >&2
 else
  printf 'This directory exists but has no completion marker; inspect its run.log before choosing another label.\n' >&2
 fi
 printf 'For a separate new run, choose an unused --label. No existing files were changed.\n' >&2
 exit 1
fi
reuse=${o[reuse-crossings]}
if [[ -n $reuse ]];then
 reuse=$(realpath -e -- "$reuse");[[ -f $reuse/COLLECTED && -f $reuse/inputs/collection.sha256 ]]||fail 'Reuse requires completed detection_efficiency collection'
 (cd -- "$reuse";sha256sum -c inputs/collection.sha256 >/dev/null)||fail 'Reused collection checksum mismatch'
 while IFS=$'\t' read -r k value;do
  [[ -v o[$k] ]]||continue
  if [[ -v given[$k] && ${o[$k]} != "$value" ]];then fail "--$k differs from saved collection";fi
  o[$k]=$value
 done < "$reuse/inputs/configuration.tsv"
fi
[[ ${o[sample]} == both || ${o[sample]} == data || ${o[sample]} == mc ]]||fail 'Invalid sample'
[[ ${o[runs]} == all || ${o[runs]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid runs'
for k in max-events min-layers stability min-common-bin min-count;do [[ ${o[$k]} =~ ^[0-9]+$ ]]||fail "Invalid integer $k";done
for k in energy max-angle max-track-rms long-margin trans-margin ambiguity-sigma long-bin trans-bin slope-bin rms-bin confidence;do [[ ${o[$k]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k";done
awk -v e="${o[energy]}" -v a="${o[max-angle]}" -v r="${o[max-track-rms]}" -v l="${o[long-margin]}" -v t="${o[trans-margin]}" -v u="${o[ambiguity-sigma]}" -v b="${o[long-bin]}" -v d="${o[trans-bin]}" -v s="${o[slope-bin]}" -v q="${o[rms-bin]}" -v cf="${o[confidence]}" 'BEGIN{exit !(e>0&&a>0&&a<.5&&r>0&&l<22.5&&t<2.5&&u>=0&&b>0&&d>0&&s>0&&q>0&&cf>0&&cf<1)}'||fail 'Invalid cuts/bin/confidence'
((o[min-layers]>=3&&o[min-layers]<=14&&o[stability]<=1&&o[min-common-bin]>=1&&o[min-count]>=2))||fail 'Invalid count/layer limit'
[[ ${o[groups]} == none || ${o[groups]} =~ ^[0-9]+(,[0-9]+)+(\;[0-9]+(,[0-9]+)+)*$ ]]||fail 'Invalid excluded groups'
[[ -z ${o[mask]} || -f ${o[mask]} ]]||fail 'Missing mask'
rows=()
if [[ -z $reuse ]];then
 if [[ -n ${o[input-manifest]} ]];then
  [[ -f ${o[input-manifest]} ]]||fail 'Missing manifest'
  while IFS= read -r row||[[ -n $row ]];do [[ -z $row || $row == \#* || $row == sample$'\t'* ]]||rows+=("$row");done < "${o[input-manifest]}"
 else
  shopt -s nullglob
  if [[ ${o[sample]} == mc ]];then files=("${o[mc-base]}/decode/mu-/mu_track/${o[mc-tag]}/${o[energy]}GeV/"ECAL_Run*.root);else files=("${o[data-base]}/decode/mu-/${o[energy]}GeV/"ECAL_Run*.root);fi
  for raw in "${files[@]}";do name=${raw##*/};[[ $name =~ ^ECAL_Run([0-9]+) ]]||continue;run=${BASH_REMATCH[1]};[[ ${o[runs]} == all || ,${o[runs]}, == *,$run,* ]]||continue
   if [[ ${o[sample]} != mc ]];then rows+=("0"$'\t'"$run"$'\t'"$raw"$'\t'"${o[data-base]}/calib/mu-/${o[energy]}GeV/$name"$'\t-'$'\t-');fi
   if [[ ${o[sample]} != data ]];then mp="mu-/mu_track/${o[mc-tag]}/${o[energy]}GeV/ECAL_Run$run.root";rows+=("1"$'\t'"$run"$'\t'"${o[mc-base]}/decode/$mp"$'\t'"${o[mc-base]}/calib/$mp"$'\t-'$'\t'"${o[mc-base]}/generate/mu-/mu_track/${o[energy]}GeV/ECAL_Run$run.root");fi
  done
 fi
 filtered=();declare -A seen=() seen_paths=();input_paths=()
 for row in "${rows[@]}";do IFS=$'\t' read -r s run raw cal unused truth extra <<< "$row"
  [[ ( $s == 0 || $s == 1 ) && $run =~ ^[0-9]+$ && $unused == - && -n $truth && -z ${extra:-} ]]||fail 'Invalid six-column manifest row'
  [[ ${o[sample]} != data || $s == 0 ]]||continue;[[ ${o[sample]} != mc || $s == 1 ]]||continue;[[ ${o[runs]} == all || ,${o[runs]}, == *,$run,* ]]||continue
  [[ ! -v seen[$s,$run] ]]||fail "Repeated sample/run $s/$run";seen[$s,$run]=1
  [[ $s == 0 || $truth != - ]]||fail 'MC needs explicit muon truth'
  for p in "$raw" "$cal" "$truth";do [[ $p == - ]]&&continue;[[ $p == /* && -f $p ]]||fail "Missing/relative input $p";input_paths+=("$p");done
  canon=$(realpath -e -- "$raw");[[ ! -v seen_paths[$canon] ]]||fail 'Repeated decode path';seen_paths[$canon]=1
  filtered+=("$row")
 done
 rows=("${filtered[@]}");((${#rows[@]}))||fail 'No actual requested muon inputs; no electron energy fallback'
 if [[ ${o[runs]} != all ]];then IFS=, read -r -a wanted <<< "${o[runs]}";for run in "${wanted[@]}";do for s in 0 1;do [[ ${o[sample]} != data || $s == 0 ]]||continue;[[ ${o[sample]} != mc || $s == 1 ]]||continue;[[ -v seen[$s,$run] ]]||fail "Requested sample/run missing: $s/$run";done;done;fi
 printf '%s\n' "${rows[@]}"
fi
printf 'Output: %s\nworkers=%s; muon energy=%s GeV; max-events/file=%s\n' "$out" "${o[workers]}" "${o[energy]}" "${o[max-events]}"
((dry==0))||exit 0
mkdir -p -- "$(dirname -- "$out")";mkdir -- "$out"||fail 'Output race: directory exists';mkdir -- "$out/inputs"
exec 8>"$out/.run.lock";flock -n 8||fail 'Output in use'
printf 'cd %q\n' "$PWD" > "$out/inputs/command.sh";printf '%q ' "$HERE/run.sh" "${args[@]}" >> "$out/inputs/command.sh";printf '\n' >> "$out/inputs/command.sh"
stat_inputs(){ for p in "${input_paths[@]}";do stat -L --printf='%n\t%s\t%y\t%i\n' -- "$p";done; }
collection_sources=("$HERE/Collect.cc" "$HERE/Common.hh" "$HERE/../muon_mip_response/TrackGeometry.hh" "$HERE/../muon_mip_response/MuonMIPCommon.hh" "$PROJECT/analysis/calibration_parameters/EBUdecode.cxx" "$PROJECT/analysis/calibration_parameters/EBUdecode.h")
{
 date -Is
 compile
 root-config --version > "$out/inputs/root_version.txt"
 tar -cf "$out/inputs/source.tar" -C "$PROJECT" analysis/channel_response/detection_efficiency/{Common.hh,Collect.cc,Analyze.cc,Draw.cc,Makefile,run.sh,submit.sh,README.md} analysis/channel_response/muon_mip_response/{TrackGeometry.hh,MuonMIPCommon.hh} analysis/channel_response/adc_energy/FastPng.hh analysis/calibration_parameters/{EBUdecode.cxx,EBUdecode.h} analysis/channel_response/{root.mk,execute.sl}
 sha256sum "${collection_sources[@]}" > "$out/inputs/collection_code.sha256"
 sha256sum "$HERE/Collect" "$HERE/Analyze" "$HERE/Draw" "$out/inputs/source.tar" > "$out/inputs/build.sha256"
 if [[ -n $reuse ]];then
  cmp "$reuse/inputs/collection_code.sha256" "$out/inputs/collection_code.sha256"||fail 'Reused collector/geometry definitions differ from current source'
  cp -- "$reuse/inputs/configuration.tsv" "$out/inputs/configuration.tsv";cp -- "$reuse/inputs/manifest.tsv" "$out/inputs/manifest.tsv"
  for name in input_audit.tsv schema.tsv sample_selection.tsv variants.tsv groups.tsv chip_readout_audit.tsv;do cp -- "$reuse/$name" "$out/$name";done
  ln -s -- "$reuse/crossings.root" "$out/crossings.root";printf '%s\n' "$reuse" > "$out/inputs/reused_from.txt"
 else
  printf '%s\n' "${rows[@]}" > "$out/inputs/manifest.tsv"
  if [[ -n ${o[mask]} ]];then cp -- "${o[mask]}" "$out/inputs/mask.tsv";o[mask]="$out/inputs/mask.tsv";fi
  for k in "${!o[@]}";do case $k in label|result-dir|workers|reuse-crossings);;*)printf '%s\t%s\n' "$k" "${o[$k]}";;esac;done | LC_ALL=C sort > "$out/inputs/configuration.tsv"
  stat_inputs > "$out/inputs/input_stats_before.tsv"
  "$HERE/Collect" "$out/inputs/manifest.tsv" "$out/inputs/configuration.tsv" "$out"
  stat_inputs > "$out/inputs/input_stats_after.tsv";cmp "$out/inputs/input_stats_before.tsv" "$out/inputs/input_stats_after.tsv"||fail 'Inputs changed while collecting'
  producer="${o[mc-base]}/../Analysis_edit/src/MCDigi.cxx"
  if [[ -f $producer ]];then cp -- "$producer" "$out/inputs/current_MCDigi.cxx";sha256sum "$producer" > "$out/inputs/current_producer.sha256";fi
 fi
 cat > "$out/inputs/provenance.tsv" <<'PROVENANCE'
item	status	evidence_limit
calibration	ID_audited_if_supplied	Run_Num,Event_Time,TriggerID/Event_Num; ambiguous reduced IDs not paired; acquisition denominator independent
truth	muon_and_energy_audited	No saved event ID; no entry-order join; truth not used as an efficiency denominator
production	historical_pairing_unverified	Paths/current_sources_and_equal_IDs_do_not_prove_a_shared_generation_receipt
CT_correction	unknown_for_input_production	No historical settings/constants receipt established; current producer copy is context only
legacy_crossings	not_reused	Existing flat.root lacks neighboring-strip outcomes, excluded-group refits, and complete identity/uncertainty fields
SSA	not_measured	Decode record count and deduplicated physical strip count are saved; virtual SSA hit count is unavailable (-1)
PROVENANCE
 (cd -- "$out";sha256sum crossings.root inputs/configuration.tsv inputs/manifest.tsv inputs/collection_code.sha256 > inputs/collection.sha256; if [[ -f inputs/mask.tsv ]];then sha256sum inputs/mask.tsv >> inputs/collection.sha256;fi)
 touch "$out/COLLECTED"
 if ((collect==0));then
  "$HERE/Analyze" "$out/crossings.root" "$out/inputs/configuration.tsv" "$out"
  if ((noplot==0));then
   pids=();for ((i=0;i<o[workers];i++));do "$HERE/Draw" "$out/statistics.root" "$out/inputs/configuration.tsv" "$out" "$i" "${o[workers]}" > "$out/plot_worker$i.log" 2>&1 & pids+=("$!");done
   failed=0;for pid in "${pids[@]}";do wait "$pid"||failed=1;done;((failed==0))||fail 'Plot worker failed; inspect plot_worker logs'
  fi
  touch "$out/COMPLETE"
 fi
 date -Is
 printf 'Completed: %s\n' "$out"
} 2>&1 | tee "$out/run.log"
