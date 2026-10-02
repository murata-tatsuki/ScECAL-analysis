#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
MC=/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
DATA=/megraid01/users/data_beamtest/ECAL_data/analysed/2023
params=$MC/Analysis_edit/share
usage(){ cat <<HELP
Usage: bash run.sh [options]
  --energy VALUE             Default: 100 (submit.sh accepts all or a list)
  --beam auto|ps|sps          Auto: PS for 0.5,1,2,3,4,5; otherwise SPS
  --sample data|mc|both       Default: both
  --mc-tag NAME              Default: threshold (MC campaign: e_ssa_noHLintercept)
  --label NAME               Default: default; completed outputs are not overwritten
  --result-dir PATH          Default: analysis/result/channel_response/calibration_residual
  --data-decode PATH --data-calib PATH --mc-decode PATH --mc-calib PATH
                             Input file or directory; files paired by basename
  --data-pedestal PATH --data-mip PATH --data-hl PATH
  --mc-pedestal PATH --mc-mip PATH --mc-hl PATH
                             Defaults: current Analysis_edit/share constants
  --data-period TEXT --mc-period TEXT
                             Documented constants validity; default: UNCONFIRMED
  --data-provenance PATH --mc-provenance PATH
                             Optional production config/log; copied and hashed
  --channels CSV             Physical CellIDs, default: all
  --layers LIST              e.g. 4,9-10; default: 0-29
  --temp-width C             Fixed temperature bins [k*w,(k+1)*w), default: 0.5
  --adc-min ADC --adc-max ADC --adc-width ADC   Ped-sub HG range: -100,4200,50
  --switch-adc ADC           Raw HG threshold, default: 2600
  --switch-window ADC        Half-width for fits, default: 200
  --switch-bin-width ADC     Default: 25; must divide switch-window
  --overlap-min ADC --overlap-max ADC    Raw HG linear-fit range: 800,2200
  --min-hits N               Plot/bin-fit minimum count, default: 30
  --events-per-file N        Evenly spaced midpoint sample (0: all, default)
  --max-files N              Per sample (0: all, default); prefer all data runs
  --plots 0|1                Default: 1; 0 saves ROOT/TSV only
  --dry-run                  Validate input paths; no build/output/submission

Each file, Run_Num and temperature bin is kept separate. Energy is read from
saved branches. Exact quantiles include tails; no MC intercept is removed.
Different constants for different periods: run separately with explicit inputs,
constants and --label. A provenance attachment is evidence, not automatic proof.
HELP
}
fail(){ echo "ERROR: $*" >&2;exit 1; }
declare -A opt=( [energy]=100 [beam]=auto [sample]=both [mc-tag]=threshold [label]=default
 [result-dir]="$HERE/../../result/channel_response/calibration_residual"
 [data-decode]='' [data-calib]='' [mc-decode]='' [mc-calib]=''
 [data-pedestal]="$params/pedestal2023_SPS.root" [mc-pedestal]="$params/pedestal2023_SPS.root"
 [data-mip]="$params/all_auto_muon_v4_trackfit.root" [mc-mip]="$params/all_auto_muon_v4_trackfit.root"
 [data-hl]="$params/all_hl_electron2023_hlratio_v3.root" [mc-hl]="$params/all_hl_electron2023_hlratio_v3.root"
 [data-period]=UNCONFIRMED [mc-period]=UNCONFIRMED [data-provenance]='' [mc-provenance]=''
 [channels]=all [layers]=0-29 [temp-width]=0.5 [adc-min]=-100 [adc-max]=4200 [adc-width]=50
 [switch-adc]=2600 [switch-window]=200 [switch-bin-width]=25 [overlap-min]=800 [overlap-max]=2200
 [min-hits]=30 [events-per-file]=0 [max-files]=0 [plots]=1 )
dry=0
while (($#));do case $1 in
 -h|--help)usage;exit 0;;--dry-run)dry=1;shift;;
 --*)key=${1#--};[[ -v opt[$key] && $# -ge 2 ]]||fail "Unknown/missing option $1";opt[$key]=$2;shift 2;;
 *)fail "Unexpected argument $1";;esac;done
for k in energy temp-width adc-min adc-max adc-width switch-adc switch-window switch-bin-width overlap-min overlap-max;do
 [[ ${opt[$k]} =~ ^-?[0-9]+([.][0-9]+)?$ ]]||fail "Invalid number $k";done
for k in min-hits events-per-file max-files;do [[ ${opt[$k]} =~ ^(0|[1-9][0-9]*)$ ]]||fail "Invalid integer $k";done
[[ ${opt[plots]} == 0 || ${opt[plots]} == 1 ]]||fail 'plots must be 0 or 1'
for k in label mc-tag;do [[ ${opt[$k]} =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]||fail "Invalid $k";done
[[ ${opt[layers]} =~ ^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$ ]]||fail 'Invalid layers'
[[ ${opt[channels]} == all || ${opt[channels]} =~ ^[0-9]+(,[0-9]+)*$ ]]||fail 'Invalid channels'
awk -v e="${opt[energy]}" -v t="${opt[temp-width]}" -v lo="${opt[adc-min]}" -v hi="${opt[adc-max]}" -v w="${opt[adc-width]}" -v cut="${opt[switch-adc]}" -v sw="${opt[switch-window]}" -v bw="${opt[switch-bin-width]}" -v ol="${opt[overlap-min]}" -v oh="${opt[overlap-max]}" -v n="${opt[min-hits]}" 'BEGIN { if (!(e>0&&t>0&&hi>lo&&w>0&&bw>0&&sw>=3*bw&&oh>ol&&oh<=cut&&n>=3)) exit 1; r=sw/bw; exit (r-int(r+0.5)>1e-8 || int(r+0.5)-r>1e-8) }' || fail 'Invalid numeric ranges'
case ${opt[sample]} in both)samples=(data mc);;data|mc)samples=("${opt[sample]}");;*)fail 'Invalid sample';;esac
if [[ ${opt[beam]} == auto ]];then case ${opt[energy]} in 0.5|1|2|3|4|5)opt[beam]=ps;;*)opt[beam]=sps;;esac;fi
[[ ${opt[beam]} == ps || ${opt[beam]} == sps ]]||fail 'Invalid beam'
[[ -n ${opt[data-decode]} ]]||opt[data-decode]=$DATA/${opt[beam]}/decode/e-/${opt[energy]}GeV
[[ -n ${opt[data-calib]} ]]||opt[data-calib]=$DATA/${opt[beam]}/calib/e-/${opt[energy]}GeV
[[ -n ${opt[mc-decode]} ]]||opt[mc-decode]=$MC/Result_MC/decode/e-/e_ssa_noHLintercept/${opt[mc-tag]}/${opt[energy]}GeV
[[ -n ${opt[mc-calib]} ]]||opt[mc-calib]=$MC/Result_MC/calib/e-/e_ssa_noHLintercept/${opt[mc-tag]}/${opt[energy]}GeV
declare -a data_pairs=() mc_pairs=()
for sample in "${samples[@]}";do
 for k in pedestal mip hl provenance;do key=$sample-$k;[[ $k != provenance || -n ${opt[$key]} ]]||continue;[[ -f ${opt[$key]} ]]||fail "Missing $key: ${opt[$key]}";opt[$key]=$(realpath -e -- "${opt[$key]}");done
 raw=${opt[$sample-decode]};cal=${opt[$sample-calib]};files=()
 if [[ -f $raw ]];then files=("$raw");elif [[ -d $raw && -d $cal ]];then shopt -s nullglob;files=("$raw/"*.root);shopt -u nullglob;else fail "Invalid pair $raw / $cal";fi
 ((${#files[@]}))||fail "No ROOT inputs: $raw"
 if ((${opt[max-files]}>0 && ${#files[@]}>${opt[max-files]}));then files=("${files[@]:0:${opt[max-files]}}");fi
 declare -n pairs=${sample}_pairs
 for a in "${files[@]}";do b=$cal;[[ ! -d $cal ]]||b=$cal/${a##*/};[[ -f $b ]]||fail "Missing paired file: $b";
  a=$(realpath -e -- "$a");b=$(realpath -e -- "$b");[[ $a$b != *$'\t'* && $a$b != *$'\n'* ]]||fail 'Tabs/newlines in paths';pairs+=("$a"$'\t'"$b");done
 printf '%s: %s file pairs\n  decode: %s\n  calib: %s\n' "$sample" "${#pairs[@]}" "$raw" "$cal"
 unset -n pairs
done
out=$(realpath -m -- "${opt[result-dir]}")/${opt[mc-tag]}/${opt[label]}/${opt[energy]}GeV
[[ ! -e $out ]]||fail "Output directory exists; choose another label: $out"
printf 'Output: %s\nTemperature width: %s C; switch: %s +/- %s raw ADC; plots: %s\n' "$out" "${opt[temp-width]}" "${opt[switch-adc]}" "${opt[switch-window]}" "${opt[plots]}"
((!dry))||exit 0
(flock 9;make -C "$HERE" all) 9>"$HERE/.build.lock"
mkdir -p -- "$(dirname -- "$out")"
mkdir -- "$out" || fail 'Output was created by another process'
mkdir -- "$out/inputs"
for k in "${!opt[@]}";do printf '%s=%q\n' "$k" "${opt[$k]}";done | sort > "$out/inputs/configuration.txt"
sha256sum "$HERE/CalibrationResidual" "$HERE/CalibrationResidual.cc" "$HERE/RootInput.hh" "$HERE/CalibrationAudit.hh" "$HERE/ResidualStats.hh" "$HERE/../adc_energy/FastPng.hh" "$HERE/run.sh" > "$out/inputs/code_checksums.txt"
root-config --version > "$out/inputs/root_version.txt"
# These snapshots identify the current extraction policy, not the saved files' history.
reference=/megraid01/users/data_beamtest/analysis/ECAL_Analysis_LCIO
for relative in src/Calibration.cxx include/Calibration.h run/global_config run_simulation/global_config;do
 source_file=$reference/$relative
 if [[ -f $source_file ]];then
  target_name=${relative//\//_}
  cp -- "$source_file" "$out/inputs/current_$target_name"
  sha256sum "$source_file" >> "$out/inputs/current_reference_checksums.txt"
 fi
done
printf 'sample\tparameter\tpath\tsha256\tdeclared_validity\tproduction_evidence\tstatus\n' > "$out/inputs/calibration_provenance.tsv"
for sample in "${samples[@]}";do
 declare -n pairs=${sample}_pairs
 printf '%s\n' "${pairs[@]}" > "$out/inputs/$sample.tsv"
 unset -n pairs
 evidence=${opt[$sample-provenance]}
 if [[ -n $evidence ]];then cp -- "$evidence" "$out/inputs/${sample}_production_evidence.txt";sha256sum "$evidence" >> "$out/inputs/evidence_checksums.txt";else evidence=UNCONFIRMED;fi
 for k in pedestal mip hl;do path=${opt[$sample-$k]};digest=$(sha256sum -- "$path");digest=${digest%% *};period=${opt[$sample-period]};
  [[ $period$path$evidence != *$'\t'* && $period$path$evidence != *$'\n'* ]]||fail 'Tabs/newlines in provenance fields'
  printf '%s\t%s\t%s\t%s\t%s\t%s\tcandidate_not_verified_production\n' "$sample" "$k" "$path" "$digest" "$period" "$evidence" >> "$out/inputs/calibration_provenance.tsv"
 done
 command=("$HERE/CalibrationResidual" --manifest "$out/inputs/$sample.tsv" --out "$out/$sample" --sample "$sample")
 for k in pedestal mip hl;do command+=("--$k" "${opt[$sample-$k]}");done
 for k in channels layers temp-width adc-min adc-max adc-width switch-adc switch-window switch-bin-width overlap-min overlap-max min-hits events-per-file plots;do command+=("--$k" "${opt[$k]}");done
 printf '%q ' "${command[@]}" >> "$out/inputs/commands.sh";printf '\n' >> "$out/inputs/commands.sh"
 "${command[@]}" 2>&1 | tee "$out/${sample}.log"
done
touch "$out/COMPLETE"
printf 'Completed: %s\n' "$out"
