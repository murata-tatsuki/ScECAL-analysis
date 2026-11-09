#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
ACTION=${2:-comparison}
[[ $# -le 2 && $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||{ echo 'Usage: bash run_digitization_impact_review.sh [LABEL] [comparison|modelscan|adccontrols]' >&2;exit 1; }
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
OUT=$BASE/digitization
[[ $ACTION == comparison || $ACTION == modelscan || $ACTION == adccontrols ]]||{ echo 'Unknown action; use comparison, modelscan or adccontrols' >&2;exit 1; }
if [[ $ACTION == adccontrols ]];then
 PILOT=$OUT
 MODEL=$PILOT/model_scan
 [[ -f $MODEL/MODEL_COMPARISON_COMPLETE && ! -e $MODEL/ADC_CONTROL_COMPLETE ]]||{ echo 'Model comparison incomplete or ADC supplement already completed' >&2;exit 1; }
 exec 8>"$MODEL/inputs/adc_control.lock"
 flock -n 8||{ echo 'ADC supplement already running' >&2;exit 1; }
 SNAP=$MODEL/inputs/adc_review_source
 mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan"
 cp -- "$HERE/DigitizationImpactReview.cc" "$HERE/run_digitization_impact_review.sh" "$HERE/ResolutionImpact.cc" "$HERE/CurrentParameters.hh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
 cp -- "$MODEL/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
 cp -- "$MODEL/inputs/reviewer_source/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
 sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* > "$MODEL/inputs/adc_code_checksums.txt"
 printf 'bash %q %q adccontrols\n' "$HERE/run_digitization_impact_review.sh" "$LABEL" > "$MODEL/inputs/adc_command.sh"
 INPUTS=("$PILOT/comparison/20mm/100GeV/paired_events.tsv" "$MODEL/inputs/manifest.tsv" "$PILOT/inputs/manifest.tsv" "$MODEL/inputs/adc_command.sh")
 while IFS=$'\t' read -r task case_name energy file_index old_raw old_cal truth new_raw new_cal;do [[ $task == task ]]&&continue;INPUTS+=("$old_raw" "$old_cal" "$truth" "$new_raw" "$new_cal");done < "$MODEL/inputs/manifest.tsv"
 while IFS=$'\t' read -r task energy file_index old_raw old_cal truth new_raw new_cal;do [[ $task == task || $energy != 100 ]]&&continue;INPUTS+=("$new_raw" "$new_cal");done < "$PILOT/inputs/manifest.tsv"
 INPUTS+=("$PILOT/inputs/reviewer_parameters/"*.root)
 sha256sum -- "${INPUTS[@]}" > "$MODEL/inputs/adc_inputs_before_sha256.txt"
 mkdir -p -- "$PROJECT/.scratch"
 TASK_TMP=$(mktemp -d "$PROJECT/.scratch/digitization_adc_controls.XXXXXXXX")
 trap 'rm -rf -- "$TASK_TMP"' EXIT
 read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
 read -r -a ROOT_LIBS <<< "$(root-config --libs)"
 g++ -std=c++20 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/DigitizationImpactReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/DigitizationImpactReview" > "$MODEL/logs/adc_review_build.log" 2>&1
 "$TASK_TMP/DigitizationImpactReview" "$BASE" "$MODEL" --adc-controls > "$MODEL/logs/adc_review.log" 2>&1
 sha256sum -- "${INPUTS[@]}" > "$MODEL/inputs/adc_inputs_after_sha256.txt"
 cmp -s "$MODEL/inputs/adc_inputs_before_sha256.txt" "$MODEL/inputs/adc_inputs_after_sha256.txt"
 touch "$MODEL/ADC_CONTROL_COMPLETE"
 printf 'Completed read-only ADC supplement: %s\n' "$MODEL"
 exit 0
fi
if [[ $ACTION == modelscan ]];then
 PILOT=$OUT
 MODEL=$PILOT/model_scan
 [[ -f $PILOT/COMPARISON_COMPLETE && -f $MODEL/inputs/manifest.tsv && ! -e $MODEL/MODEL_COMPARISON_COMPLETE ]]||{ echo 'CTbaseline/modelmanifest incomplete or model comparison already completed' >&2;exit 1; }
 exec 8>"$MODEL/inputs/comparison.lock"
 flock -n 8||{ echo 'Model reviewer is already running' >&2;exit 1; }
 sha256sum -c "$MODEL/inputs/source_parameter_checksums.txt" > "$MODEL/inputs/source_verification.log"
 INPUTS=("$MODEL/inputs/manifest.tsv" "$PILOT/comparison/20mm/100GeV/paired_events.tsv" "$BASE/integrated/comparison/20mm/100GeV/events.tsv" "$BASE/100GeV/baseline/reference.tsv")
 while IFS=$'\t' read -r task case_name energy file_index old_raw old_cal truth new_raw new_cal;do
  [[ $task == task ]]&&continue
  [[ -f $MODEL/mc/$case_name/100GeV/${new_raw##*/}.COMPLETE ]]||{ echo "Model output incomplete: $new_raw" >&2;exit 1; }
  INPUTS+=("$old_raw" "$old_cal" "$truth" "$new_raw" "$new_cal")
 done < "$MODEL/inputs/manifest.tsv"
 while IFS=$'\t' read -r task energy file_index old_raw old_cal truth new_raw new_cal;do [[ $task == task || $energy != 100 ]]&&continue;INPUTS+=("$new_raw" "$new_cal");done < "$PILOT/inputs/manifest.tsv"
 INPUTS+=("$PILOT/inputs/reviewer_parameters/"*.root)
 stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$MODEL/inputs/comparison_stats_before.tsv"
 sha256sum -- "${INPUTS[@]}" > "$MODEL/inputs/comparison_sha256_before.txt"
 SNAP=$MODEL/inputs/reviewer_source
 mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan"
 cp -- "$HERE/DigitizationImpactReview.cc" "$HERE/run_digitization_impact_review.sh" "$HERE/ResolutionImpact.cc" "$HERE/CurrentParameters.hh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
 cp -- "$PILOT/inputs/reviewer_source/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
 cp -- "$PILOT/inputs/reviewer_source/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
 sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* > "$MODEL/inputs/reviewer_code_checksums.txt"
 printf 'bash %q %q modelscan\n' "$HERE/run_digitization_impact_review.sh" "$LABEL" > "$MODEL/inputs/review_command.sh"
 mkdir -p -- "$PROJECT/.scratch"
 TASK_TMP=$(mktemp -d "$PROJECT/.scratch/digitization_model_review.XXXXXXXX")
 trap 'rm -rf -- "$TASK_TMP"' EXIT
 read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
 read -r -a ROOT_LIBS <<< "$(root-config --libs)"
 g++ -std=c++20 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/DigitizationImpactReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/DigitizationImpactReview" > "$MODEL/logs/reviewer_build.log" 2>&1
 for case_name in pixels14400 pixels19600 threshold_plus5 threshold_plus10;do
  CASE_OUT=$MODEL/comparison/$case_name
  [[ ! -e $CASE_OUT/COMPARISON_COMPLETE ]]||{ echo "Case comparison exists: $case_name" >&2;exit 1; }
  mkdir -p -- "$CASE_OUT/inputs"
  ln -s -- "$PILOT/inputs/reviewer_parameters" "$CASE_OUT/inputs/reviewer_parameters"
  ln -s -- "$SNAP" "$CASE_OUT/inputs/reviewer_source"
  awk -F '\t' -v name="$case_name" 'BEGIN{OFS="\t";print "task","energy_GeV","file_index","old_decode","old_calib","truth","new_decode","new_calib"} NR>1&&$2==name{print $1,$3,$4,$5,$6,$7,$8,$9}' "$MODEL/inputs/manifest.tsv" > "$CASE_OUT/inputs/manifest.tsv"
  awk -F '\t' -v name="$case_name" 'NR==1||$1==name' "$MODEL/inputs/cases.tsv" > "$CASE_OUT/inputs/model_condition.tsv"
  "$TASK_TMP/DigitizationImpactReview" "$BASE" "$CASE_OUT" 100 > "$MODEL/logs/review_$case_name.log" 2>&1
  touch "$CASE_OUT/COMPARISON_COMPLETE"
 done
 "$TASK_TMP/DigitizationImpactReview" "$BASE" "$MODEL" --model-summary > "$MODEL/logs/review_summary.log" 2>&1
 stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$MODEL/inputs/comparison_stats_after.tsv"
 sha256sum -- "${INPUTS[@]}" > "$MODEL/inputs/comparison_sha256_after.txt"
 cmp -s "$MODEL/inputs/comparison_stats_before.tsv" "$MODEL/inputs/comparison_stats_after.tsv"
 cmp -s "$MODEL/inputs/comparison_sha256_before.txt" "$MODEL/inputs/comparison_sha256_after.txt"
 touch "$MODEL/MODEL_COMPARISON_COMPLETE"
 printf 'Completed bounded100GeV model comparison: %s\n' "$MODEL"
 exit 0
fi
[[ ! -e $OUT/COMPARISON_COMPLETE ]]||{ echo "Completed comparison exists; choose a new label: $OUT" >&2;exit 1; }
MANIFEST=$OUT/inputs/manifest.tsv
[[ -f $MANIFEST ]]||{ echo "Missing frozen digitization manifest: $MANIFEST" >&2;exit 1; }
exec 8>"$OUT/inputs/comparison.lock"
flock -n 8||{ echo 'Another process is writing the paired comparison' >&2;exit 1; }
sha256sum -c "$OUT/inputs/source_checksums.txt" > "$OUT/inputs/frozen_production_source_verification.log"
INPUTS=("$MANIFEST")
while IFS=$'\t' read -r task energy file_index old_raw old_cal truth new_raw new_cal;do
 [[ $task == task ]]&&continue
 for path in "$old_raw" "$old_cal" "$truth" "$new_raw" "$new_cal";do [[ -f $path ]]||{ echo "Missing pilot input: $path" >&2;exit 1; };INPUTS+=("$path");done
 dest=${new_raw%/decode/*}
 [[ -f $dest/${new_raw##*/}.COMPLETE ]]||{ echo "Production pilot generation incomplete: $new_raw" >&2;exit 1; }
done < "$MANIFEST"
for energy in 5 20 100;do INPUTS+=("$BASE/integrated/comparison/20mm/${energy}GeV/events.tsv" "$BASE/${energy}GeV/baseline/reference.tsv");done
PARAM_BASE=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share
PARAM_FILES=("$PARAM_BASE/all_hl_electron2023_hlratio_v3.root" "$PARAM_BASE/all_auto_muon_v4_trackfit.root" "$PARAM_BASE/pedestal2023_SPS.root" "$PARAM_BASE/20230504_1425_Calib_LED_all.root" "$PROJECT/analysis/result/threshold/threshold.root")
INPUTS+=("${PARAM_FILES[@]}")
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_stats_before.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_sha256_before.txt"
SNAP=$OUT/inputs/reviewer_source
mkdir -p -- "$SNAP/investigations/e_ssa_threshold_review_20260930" "$SNAP/analysis/impact_studies" "$SNAP/analysis/cog_channel_scan" "$OUT/inputs/reviewer_parameters"
cp -- "$HERE/DigitizationImpactReview.cc" "$HERE/run_digitization_impact_review.sh" "$HERE/ResolutionImpact.cc" "$HERE/CurrentParameters.hh" "$SNAP/investigations/e_ssa_threshold_review_20260930/"
cp -- "$PROJECT/analysis/impact_studies/StudyCommon.hh" "$SNAP/analysis/impact_studies/"
cp -- "$PROJECT/analysis/cog_channel_scan/SingleEnergyAnalysis.cc" "$SNAP/analysis/cog_channel_scan/"
cp -- "${PARAM_FILES[@]}" "$OUT/inputs/reviewer_parameters/"
sha256sum "$SNAP/investigations/e_ssa_threshold_review_20260930/"* "$SNAP/analysis/impact_studies/"* "$SNAP/analysis/cog_channel_scan/"* "$OUT/inputs/reviewer_parameters/"* > "$OUT/inputs/reviewer_source_parameter_checksums.txt"
printf 'bash %q %q\n' "$HERE/run_digitization_impact_review.sh" "$LABEL" > "$OUT/inputs/comparison_command.sh"
mkdir -p -- "$PROJECT/.scratch"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/digitization_impact_review.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++17 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$SNAP/investigations/e_ssa_threshold_review_20260930/DigitizationImpactReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/DigitizationImpactReview" > "$OUT/comparison_build.log" 2>&1
"$TASK_TMP/DigitizationImpactReview" "$BASE" "$OUT" > "$OUT/comparison_run.log" 2>&1
stat -L --printf='%n\t%s\t%Y\n' -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_stats_after.tsv"
sha256sum -- "${INPUTS[@]}" > "$OUT/inputs/comparison_input_sha256_after.txt"
cmp -s -- "$OUT/inputs/comparison_input_stats_before.tsv" "$OUT/inputs/comparison_input_stats_after.tsv"
cmp -s -- "$OUT/inputs/comparison_input_sha256_before.txt" "$OUT/inputs/comparison_input_sha256_after.txt"
touch "$OUT/COMPARISON_COMPLETE"
printf 'Completed two-file production digitization impact comparison: %s\n' "$OUT"
