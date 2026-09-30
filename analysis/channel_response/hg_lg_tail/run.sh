#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
MC=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
DATA=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
params=$MC/Analysis_edit/share
usage(){ cat <<HELP
Usage: bash run.sh [options]
  --energy VALUE            One beam energy (default: 100); use submit.sh for all
  --beam auto|ps|sps         Data beam (default: auto, PS for 0.5-5 GeV)
  --mc-tag NAME             MC subdirectory (default: threshold)
  --label NAME              Output label (default: default)
  --result-dir PATH         Default: analysis/result/channel_response/hg_lg_tail
  --data-decode PATH         Raw data ROOT file or directory
  --data-calib PATH          Corresponding calibrated file or directory
  --mc-decode PATH           Raw MC ROOT file or directory
  --mc-calib PATH            Corresponding calibrated file or directory
  --pedestal PATH --hl PATH --mip PATH    Calibration parameter files
  --bad-channel-source PATH  comparison source file or none
  --max-events N            Maximum events per sample, across files (0: all)
  --max-files N             Maximum files per sample (0: all)
  --events-per-file N       Deterministic 10-block sample per file (0: all)
  --hg-min ADC --hg-max ADC Overlap range in raw HG (default: 800, 2200)
  --nsigma N               Tail threshold after per-channel/sample median centering (5)
  --min-hits N             Required overlap hits in EACH sample/channel (100)
  --min-tail-channels N    Distinct low-tail channels for event flag (3)
  --selection MODE         all | cog20 | shower-cog20 (default)
  --reference-energy GEV   Shower reference; 0 uses sampled MC good-energy median
  --scratch-root PATH      Local cache location, removed on exit (default: /tmp)
  --dry-run                Validate inputs and show configuration; no writes/submission

All connected channels are inspected. A single common channel set and tail rule
are used for both hit and event counts. No analysis output from investigations/
or impact_studies is required. Figures are PNG only; trees/histograms are ROOT.
HELP
}
fail(){ echo "ERROR: $*" >&2; exit 1; }
declare -A opt=( [energy]=100 [beam]=auto [mc-tag]=threshold [label]=default
 [result-dir]="$HERE/../../result/channel_response/hg_lg_tail"
 [data-decode]='' [data-calib]='' [mc-decode]='' [mc-calib]=''
 [pedestal]="$params/pedestal2023_SPS.root" [hl]="$params/all_hl_electron2023_hlratio_v3.root" [mip]="$params/all_auto_muon_v4_trackfit.root"
 [bad-channel-source]="$HERE/../../comparison/SingleEnergyAnalysis.cc"
 [max-events]=0 [max-files]=0 [events-per-file]=0 [hg-min]=800 [hg-max]=2200 [nsigma]=5 [min-hits]=100 [min-tail-channels]=3
 [selection]=shower-cog20 [reference-energy]=0 [scratch-root]=/tmp )
dry=0
while (($#));do
 case $1 in -h|--help) usage;exit 0;; --dry-run) dry=1;shift;; --*) key=${1#--};[[ -v opt[$key] && $# -ge 2 ]]||fail "Unknown/missing option $1";opt[$key]=$2;shift 2;; *) fail "Unexpected argument $1";; esac
done
[[ ${opt[energy]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail 'energy must be one positive number'
for k in mc-tag label;do [[ ${opt[$k]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]||fail "Invalid $k";done
for k in max-events max-files events-per-file min-hits min-tail-channels;do [[ ${opt[$k]} =~ ^(0|[1-9][0-9]*)$ ]]||fail "Invalid integer $k";done
((opt[min-hits]>=2&&opt[min-tail-channels]>=1))||fail 'min-hits >=2 and min-tail-channels >=1 required'
for k in hg-min hg-max nsigma reference-energy;do [[ ${opt[$k]} =~ ^[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k";done
awk -v e="${opt[energy]}" -v lo="${opt[hg-min]}" -v hi="${opt[hg-max]}" -v sig="${opt[nsigma]}" 'BEGIN{exit !(e>0&&lo<hi&&hi<=2600&&sig>0)}' || fail 'Invalid beam/HG range/nsigma'
case ${opt[selection]} in all|cog20|shower-cog20);;*) fail 'Invalid selection';;esac
if [[ ${opt[beam]} == auto ]];then case ${opt[energy]} in 0.5|1|2|3|4|5) opt[beam]=ps;;*) opt[beam]=sps;;esac;fi
[[ ${opt[beam]} == ps || ${opt[beam]} == sps ]]||fail 'Invalid beam'
[[ -n ${opt[data-decode]} ]]||opt[data-decode]=$DATA/${opt[beam]}/decode/e-/${opt[energy]}GeV
[[ -n ${opt[data-calib]} ]]||opt[data-calib]=$DATA/${opt[beam]}/calib/e-/${opt[energy]}GeV
[[ -n ${opt[mc-decode]} ]]||opt[mc-decode]=$MC/Result_MC/decode/e-/sps/${opt[mc-tag]}/${opt[energy]}GeV
[[ -n ${opt[mc-calib]} ]]||opt[mc-calib]=$MC/Result_MC/calib/e-/sps/${opt[mc-tag]}/${opt[energy]}GeV
for k in pedestal hl mip bad-channel-source;do
 [[ $k != bad-channel-source || ${opt[$k]} != none ]]||continue
 [[ -f ${opt[$k]} ]]||fail "Missing $k: ${opt[$k]}";opt[$k]=$(realpath -e -- "${opt[$k]}")
done
declare -a data_pairs=() mc_pairs=()
for sample in data mc;do
 raw=${opt[$sample-decode]};cal=${opt[$sample-calib]};files=()
 if [[ -f $raw ]];then files=("$raw");elif [[ -d $raw && -d $cal ]];then shopt -s nullglob;files=("$raw/"*.root);shopt -u nullglob;else fail "Invalid input pair: $raw / $cal";fi
 ((${#files[@]}))||fail "No ROOT files: $raw"
 if ((${opt[max-files]}>0 && ${#files[@]}>${opt[max-files]}));then files=("${files[@]:0:${opt[max-files]}}");fi
 declare -n pairs=${sample}_pairs
 for a in "${files[@]}";do b=$cal;[[ ! -d $cal ]]||b=$cal/${a##*/};[[ -f $b ]]||fail "Missing paired file: $b";a=$(realpath -e -- "$a");b=$(realpath -e -- "$b");[[ $a$b != *$'\t'* && $a$b != *$'\n'* ]]||fail 'Tabs/newlines in input paths';pairs+=("$a"$'\t'"$b");done
 printf '%s: %s file pair(s), beam=%s, energy=%s GeV\n' "$sample" "${#pairs[@]}" "${opt[beam]}" "${opt[energy]}"
 unset -n pairs
done
out=$(realpath -m -- "${opt[result-dir]}")/${opt[mc-tag]}/${opt[label]}/${opt[energy]}GeV
[[ ! -e $out/tail_study.root ]]||fail "Completed output exists; choose another --label: $out"
printf 'Output: %s\nSelection: %s; all channels; max-events=%s; common min-hits=%s; low-tail z < -%s; event flag >=%s channels\n' "$out" "${opt[selection]}" "${opt[max-events]}" "${opt[min-hits]}" "${opt[nsigma]}" "${opt[min-tail-channels]}"
((!dry))||exit 0
( flock 9; make -C "$HERE" all ) 9>"$HERE/.build.lock"
mkdir -p -- "$out/inputs" "${opt[scratch-root]}"
exec 8>"$out/.run.lock";flock -n 8||fail 'Another job is writing this output'
[[ ! -e $out/tail_study.root ]]||fail 'Completed output appeared while waiting'
scratch=$(mktemp -d "${opt[scratch-root]%/}/scecal_hg_lg_tail.XXXXXX")
trap 'rm -rf -- "$scratch"' EXIT
printf '%s\n' "${data_pairs[@]}" > "$out/inputs/data.tsv"
printf '%s\n' "${mc_pairs[@]}" > "$out/inputs/mc.tsv"
for k in "${!opt[@]}";do printf '%s=%s\n' "$k" "${opt[$k]}";done | sort > "$out/inputs/configuration.txt"
sha256sum "$HERE/HGLGTailStudy.cc" "$HERE/HGLGTailLayerMaps.hh" "$HERE/../../calibration_parameters/EBUdecode.cxx" "$HERE/run.sh" "$HERE/../../impact_studies/StudyCommon.hh" "${opt[pedestal]}" "${opt[hl]}" "${opt[mip]}" > "$out/inputs/checksums.txt"
if [[ ${opt[bad-channel-source]} != none ]];then cp -- "${opt[bad-channel-source]}" "$out/inputs/comparison_source.cc";opt[bad-channel-source]=$out/inputs/comparison_source.cc;fi
command=("$HERE/HGLGTailStudy" --data-manifest "$out/inputs/data.tsv" --mc-manifest "$out/inputs/mc.tsv" --out "$out" --scratch "$scratch")
for k in energy pedestal hl mip bad-channel-source max-events max-files events-per-file hg-min hg-max nsigma min-hits min-tail-channels selection reference-energy;do command+=("--$k" "${opt[$k]}");done
"${command[@]}" 2>&1 | tee "$out/run.log"
