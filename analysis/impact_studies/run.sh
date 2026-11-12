#!/bin/bash
set -euo pipefail
export LC_ALL=C
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
MC=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
DATA=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
declare -A opt=([study]=hl_intercept [energy]=100 [beam]=auto [mc-tag]=threshold [label]=default
 [max-events]=0 [plot-channels]=all [result-dir]="$HERE/../result/impact_studies"
 [hl]="$MC/Analysis_edit/share/all_hl_electron2023_hlratio_v3.root"
 [mip]="$MC/Analysis_edit/share/all_auto_muon_v4_trackfit.root" [pedestal]="$MC/Analysis_edit/share/pedestal2023_SPS.root"
 [data-decode]='' [data-calib]='' [mc-decode]='' [mc-calib]='')
plots=1 dry=0
fail(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
usage(){ cat <<EOF
Usage: bash run.sh [options]
  --study hl_intercept     Study name (currently only hl_intercept)
  --energy VALUE|all       Beam energy in GeV (default: 100); all runs the full
                         beam-energy list sequentially, then summarizes it
  --beam auto|ps|sps       Data beam; auto uses ps at 0.5,1,2,3,4,5 GeV
  --mc-tag NAME           MC condition (default: threshold)
  --label NAME            Named run/configuration (default: default)
  --max-events N          Per sample, across files; 0 = all (default)
  --plot-channels IDS     all (default), or physical CellIDs separated by commas
                         Changes detailed PNG selection ONLY; all hits are analysed
  --no-plots              Save ROOT and TSV only
  --data-decode PATH      Override raw data directory or ROOT file
  --data-calib PATH       Override matching calibrated data directory or file
  --mc-decode PATH        Override raw MC directory or ROOT file
  --mc-calib PATH         Override matching calibrated MC directory or file
  --hl PATH              MC InterCalib parameter ROOT
  --mip PATH             MC MIP_Fit parameter ROOT
  --pedestal PATH        MC ChnLevel pedestal ROOT (double vectors)
  --result-dir PATH      Default: analysis/result/impact_studies
  --dry-run              Validate inputs and print command; no outputs/build
  -h, --help             Show this help

Output: <result-dir>/hl_intercept/<mc-tag>/<label>/<energy>GeV/
Never overwrites a completed result. Use a different --label for another run.
MC low-energy paths are under sps, including the PS energies.
EOF
}
original_args=("$@")
while (($#)); do case $1 in
 -h|--help) usage; exit 0;; --no-plots) plots=0;shift;; --dry-run) dry=1;shift;;
 --*) key=${1#--};[[ -v opt[$key] ]]||fail "Unknown option $1";(($#>=2))||fail "Missing value for $1";opt[$key]=$2;shift 2;;
 *) fail "Unexpected argument $1";; esac;done
if [[ ${opt[energy]} == all ]];then
 exec bash "$HERE/run_all.sh" "${original_args[@]}"
fi
[[ ${opt[study]} == hl_intercept ]]||fail 'Only hl_intercept is currently implemented'
[[ ${opt[energy]} =~ ^[0-9]+(\.[0-9]+)?$ ]]||fail 'Invalid energy'
awk -v e="${opt[energy]}" 'BEGIN{exit !(e>0 && e<1e6)}'||fail 'Energy must be positive'
[[ ${opt[max-events]} =~ ^[0-9]{1,18}$ ]]||fail 'Invalid max-events'
for key in mc-tag label;do [[ ${opt[$key]} =~ ^[A-Za-z0-9_-]+$ ]]||fail "Invalid $key";done
[[ ${opt[plot-channels]} == all || ${opt[plot-channels]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid plot-channels'
case ${opt[beam]} in auto) case ${opt[energy]} in 0.5|1|2|3|4|5) opt[beam]=ps;; *) opt[beam]=sps;; esac;; ps|sps) :;; *) fail 'Invalid beam';; esac
for key in hl mip pedestal;do [[ -f ${opt[$key]} ]]||fail "Missing ${opt[$key]}";opt[$key]=$(realpath -e -- "${opt[$key]}");done
out=$(realpath -m -- "${opt[result-dir]}")/${opt[study]}/${opt[mc-tag]}/${opt[label]}/${opt[energy]}GeV
[[ ! -e $out/study.root ]]||fail "Completed output exists: $out; use another --label"
declare -a data_pairs=() mc_pairs=()
for sample in data mc;do
 for stage in decode calib;do key=$sample-$stage;if [[ -z ${opt[$key]} ]];then
  if [[ $sample == data ]];then opt[$key]=$DATA/${opt[beam]}/$stage/e-/${opt[energy]}GeV;else opt[$key]=$MC/Result_MC/$stage/e-/sps/${opt[mc-tag]}/${opt[energy]}GeV;fi
 fi;done
 raw=${opt[$sample-decode]} cal=${opt[$sample-calib]}
 [[ -e $raw && -e $cal ]]||fail "Missing input: $raw or $cal"
 files=();if [[ -f $raw ]];then files=("$raw");else [[ -d $cal ]]||fail 'Decode directory requires calib directory';shopt -s nullglob;files=("$raw"/*.root);shopt -u nullglob;fi
 ((${#files[@]}))||fail "No ROOT files: $raw"
 declare -n pairs=${sample}_pairs
 for a in "${files[@]}";do b=$cal;[[ ! -d $cal ]]||b=$cal/${a##*/};[[ -f $b ]]||fail "Missing calibration $b";a=$(realpath -e -- "$a");b=$(realpath -e -- "$b");[[ $a$b != *$'\t'* && $a$b != *$'\n'* ]]||fail 'Tabs/newlines in paths are unsupported';pairs+=("$a"$'\t'"$b");done
 printf '%s: %s input pairs\n' "$sample" "${#pairs[@]}";unset -n pairs
done
command=("$HERE/HLInterceptStudy" --mc-manifest "$out/inputs/mc.tsv" --data-manifest "$out/inputs/data.tsv" --out "$out" --plots "$plots")
for key in energy max-events plot-channels hl mip pedestal;do command+=("--$key" "${opt[$key]}");done
printf '%q ' "${command[@]}";printf '\n'
((!dry))||exit 0
make -C "$HERE" HLInterceptStudy
mkdir -p -- "$out"
exec 9>"$out/.lock"
flock -n 9||fail "Another process is using $out"
[[ ! -e $out/study.root ]]||fail 'Output completed by another process'
mkdir -p -- "$out/inputs"
printf '%s\n' "${data_pairs[@]}" > "$out/inputs/data.tsv"
printf '%s\n' "${mc_pairs[@]}" > "$out/inputs/mc.tsv"
{ for key in "${!opt[@]}";do printf '%s=%q\n' "$key" "${opt[$key]}";done;printf 'plots=%s\n' "$plots";printf '%q ' "${command[@]}";printf '\n';sha256sum -- "${opt[hl]}" "${opt[mip]}" "${opt[pedestal]}" "$HERE/HLInterceptStudy.cc" "$HERE/StudyCommon.hh" "$HERE/EventHistograms.hh"; } > "$out/inputs/configuration.txt"
"${command[@]}"
