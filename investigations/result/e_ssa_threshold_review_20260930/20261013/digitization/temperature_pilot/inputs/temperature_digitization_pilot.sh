#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
ACTION=${2:-prepare}
[[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||exit 1
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
PILOT=$BASE/digitization
OUT=$PILOT/temperature_pilot
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
if [[ $ACTION == prepare ]];then
 [[ -f $PILOT/COMPARISON_COMPLETE && ! -e $OUT/inputs/jobs.tsv ]]||{ echo 'MissingCTbaseline orpilotexists' >&2;exit 1; }
 mkdir -p -- "$OUT/inputs/source/src" "$OUT/inputs/source/include" "$OUT/inputs/binaries" "$OUT/logs"
 cp -- "$PILOT/inputs/digitizer_source/"{main.cxx,CMakeLists.txt} "$OUT/inputs/source/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.cxx "$OUT/inputs/source/src/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.h "$OUT/inputs/source/include/"
 cp -- "$HERE/TemperaturePhysical.hh" "$OUT/inputs/source/include/"
 cp -- "$HERE/temperature_digitization_pilot.sh" "$OUT/inputs/"
 perl -0pi -e 's/#include "Extract.h"/#include "Extract.h"\n#include "TemperaturePhysical.hh"/;s/(int SiPMType = \(\(layerNo < 4\) \|\| \(layerNo >27\)\) \? 1 : 0;)/$1\n Point thermalPosition = CellID2Position(cellID);\n double thermalT = TemperaturePhysical::local(layerNo,thermalPosition.posX,thermalPosition.posY,TemperaturePhysical::current);/;s/(photonNo = round\(photonNo\);)/$1\n pCrosstalk10 = TemperaturePhysical::pct(4,thermalT);\n pCrosstalk15 = TemperaturePhysical::pct(0,thermalT);\n pde10 *= TemperaturePhysical::pdeRatio(4,thermalT);\n pde15 *= TemperaturePhysical::pdeRatio(0,thermalT);/;s/(double ADCMean  = peNo \* ex.SiPMGain\[layerNo\]\[chipNo\]\[chanNo\])/$1 * TemperaturePhysical::gainRatio(layerNo,thermalT)/;s/InDigiBranchClear\(\);/InDigiBranchClear();\n TemperaturePhysical::set(i_entry);/;s/_digiTemperature.assign\(30, vector<double>\(16, 20\)\);/_digiTemperature = TemperaturePhysical::current;/' "$OUT/inputs/source/src/MCDigi.cxx"
 perl -0pi -e 's/#include "Extract.h"/#include "Extract.h"\n#include "TemperaturePhysical.hh"/;s/(EXTR.SiPMGainExtract\(argv\[5\]\);)/$1\n TemperaturePhysical::initialize(EXTR,argv[5]);/' "$OUT/inputs/source/main.cxx"
 grep -q 'TemperaturePhysical::set(i_entry)' "$OUT/inputs/source/src/MCDigi.cxx"
 grep -q 'TemperaturePhysical::gainRatio(layerNo,thermalT)' "$OUT/inputs/source/src/MCDigi.cxx"
 SNAP=$OUT/inputs/reviewer_source
 mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan"
 cp -- "$HERE/"{TemperatureDigitizationReview.cc,DigitizationImpactReview.cc,ResolutionImpact.cc,CurrentParameters.hh} "$SNAP/investigations/e_ssa_threshold_review_20260930/"
 cp -- "$PILOT/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
 cp -- "$PILOT/inputs/reviewer_source/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
 g++ -std=c++20 -O2 "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/TemperatureDigitizationReview.cc" "${ROOT_LIBS[@]}" -o "$OUT/inputs/binaries/TemperatureDigitizationReview" > "$OUT/logs/reviewer_build.log" 2>&1
 "$OUT/inputs/binaries/TemperatureDigitizationReview" "$BASE" "$OUT" bank > "$OUT/logs/bank.log" 2>&1
 printf 'task\tcase\tenergy_GeV\tfile_index\told_decode\told_calib\ttruth\tnew_decode\tnew_calib\n' > "$OUT/inputs/manifest.tsv"
 task=0
 for condition in ref20 real_pde_constant real_pde_linear;do
  real_temp=0;linear_pde=0
  [[ $condition != ref20 ]]&&real_temp=1
  [[ $condition == real_pde_linear ]]&&linear_pde=1
  g++ -std=c++20 -O2 "${ROOT_CFLAGS[@]}" -DPHYS_TEMP_REAL="$real_temp" -DPHYS_PDE_LINEAR="$linear_pde" -I"$OUT/inputs/source/include" "$OUT/inputs/source/main.cxx" "$OUT/inputs/source/src/"{MCDigi,Extract}.cxx "${ROOT_LIBS[@]}" -o "$OUT/inputs/binaries/MCAnalysis_$condition" > "$OUT/logs/build_$condition.log" 2>&1
  while IFS=$'\t' read -r index energy fi oldraw oldcal truth raw cal;do
   [[ $index == task || $energy != 100 ]]&&continue
   dest=$OUT/mc/$condition/100GeV
   mkdir -p -- "$dest/decode" "$dest/calib"
   printf '%s\t%s\t100\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$condition" "$fi" "$oldraw" "$oldcal" "$truth" "$dest/decode/${truth##*/}" "$dest/calib/${truth##*/}" >> "$OUT/inputs/manifest.tsv"
   task=$((task+1))
  done < "$PILOT/inputs/manifest.tsv"
 done
 [[ $task == 6 ]]||exit 1
 sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* > "$OUT/inputs/reviewer_code_checksums.txt"
 sha256sum "$OUT/inputs/source/"*.* "$OUT/inputs/source/src/"* "$OUT/inputs/source/include/"* "$OUT/inputs/binaries/"* "$OUT/inputs/temperature_digitization_pilot.sh" "$OUT/inputs/temperature_bank.root" "$PILOT/inputs/calib" "$PILOT/inputs/reviewer_parameters/"* > "$OUT/inputs/source_parameter_checksums.txt"
 printf 'Fixed-optical-yieldtoythermalpilot100GeVonly; temperaturesjointdata30x16bank; reference20C gain correctedusingLEDtemp; pCT quadraticovervoltage; PDEconstant/linearbracket; calibrationproduction unchanged; MIPinput unchanged;1-p0 normalisation20C only; independent RNG SetSeed(0).\n' > "$OUT/inputs/model_conditions.txt"
 job=$(sbatch --parsable --partition=all --array=0-5 --mem=4G --time=02:00:00 --cpus-per-task=1 --ntasks=1 --job-name=temperature_pilot --output="$OUT/logs/slurm-%A_%a.out" --error="$OUT/logs/slurm-%A_%a.err" --wrap="bash '$HERE/temperature_digitization_pilot.sh' '$LABEL' worker")
 printf '%s\t%s\n' "$(date -Is)" "$job" > "$OUT/inputs/jobs.tsv"
 echo "Submitted bounded6task thermal diagnostic $job"
elif [[ $ACTION == worker ]];then
 index=${SLURM_ARRAY_TASK_ID:?}
 line=$(sed -n "$((index+2))p" "$OUT/inputs/manifest.tsv")
 IFS=$'\t' read -r task condition energy fi oldraw oldcal truth raw cal <<< "$line"
 [[ $task == "$index" && ! -e $raw && ! -e $cal ]]||exit 1
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$oldraw" "$oldcal" > "$OUT/logs/input_${task}_before.tsv"
 export TEMPERATURE_BANK=$OUT/inputs/temperature_bank.root
 export TEMPERATURE_OFFSET=0
 [[ $fi == 5 ]]&&TEMPERATURE_OFFSET=1000
 PARAM=$PILOT/inputs/reviewer_parameters
 "$OUT/inputs/binaries/MCAnalysis_$condition" "$truth" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$PARAM/20230504_1425_Calib_LED_all.root" "$PARAM/threshold.root" "$raw" 0
 "$PILOT/inputs/calib" "$raw" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$cal"
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$oldraw" "$oldcal" > "$OUT/logs/input_${task}_after.tsv"
 cmp -s "$OUT/logs/input_${task}_before.tsv" "$OUT/logs/input_${task}_after.tsv"
 touch "$OUT/mc/$condition/100GeV/${truth##*/}.COMPLETE"
elif [[ $ACTION == review ]];then
 [[ -f $OUT/inputs/jobs.tsv && ! -e $OUT/COMPARISON_COMPLETE ]]||exit 1
 while IFS=$'\t' read -r task condition energy fi oldraw oldcal truth raw cal;do
  [[ $task == task ]]&&continue
  [[ -f $OUT/mc/$condition/100GeV/${truth##*/}.COMPLETE ]]||{ echo 'Pilotjobsnotcomplete' >&2;exit 1; }
 done < "$OUT/inputs/manifest.tsv"
 sha256sum -c "$OUT/inputs/source_parameter_checksums.txt" > "$OUT/inputs/source_verification.log"
 sha256sum -c "$OUT/inputs/reviewer_code_checksums.txt" > "$OUT/inputs/reviewer_code_verification.log"
 INPUTS=("$OUT/inputs/manifest.tsv" "$OUT/inputs/temperature_bank.root" "$BASE/integrated/comparison/20mm/100GeV/events.tsv" "$PILOT/comparison/20mm/100GeV/paired_events.tsv")
 while IFS=$'\t' read -r task condition energy fi oldraw oldcal truth raw cal;do [[ $task == task ]]&&continue;INPUTS+=("$raw" "$cal" "$truth" "$oldraw" "$oldcal");done < "$OUT/inputs/manifest.tsv"
 sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/review_inputs_before_sha256.txt"
 "$OUT/inputs/binaries/TemperatureDigitizationReview" "$BASE" "$OUT" review > "$OUT/logs/review.log" 2>&1
 sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/review_inputs_after_sha256.txt"
 cmp -s "$OUT/inputs/review_inputs_before_sha256.txt" "$OUT/inputs/review_inputs_after_sha256.txt"
 touch "$OUT/COMPARISON_COMPLETE"
 echo "Thermal diagnostic complete: $OUT"
elif [[ $ACTION == plots ]];then
 [[ -f $OUT/COMPARISON_COMPLETE ]]||exit 1
 sha256sum "$OUT/"*.tsv "$OUT/temperature_resolution.root" > "$OUT/inputs/plot_inputs_before_sha256.txt"
 "$OUT/inputs/binaries/TemperatureDigitizationReview" "$BASE" "$OUT" plots > "$OUT/logs/plots.log" 2>&1
 sha256sum "$OUT/"*.tsv "$OUT/temperature_resolution.root" > "$OUT/inputs/plot_inputs_after_sha256.txt"
 cmp -s "$OUT/inputs/plot_inputs_before_sha256.txt" "$OUT/inputs/plot_inputs_after_sha256.txt"
else
 echo 'Usage: temperature_digitization_pilot.sh LABEL prepare|worker|review|plots'  >&2;exit 1
fi
