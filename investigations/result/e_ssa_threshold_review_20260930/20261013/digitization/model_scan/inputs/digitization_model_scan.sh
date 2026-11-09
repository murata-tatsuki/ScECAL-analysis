#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}
ACTION=${2:-submit}
[[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||{ echo 'Invalid label' >&2;exit 1; }
PILOT=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL/digitization
OUT=$PILOT/model_scan
if [[ $ACTION == submit ]];then
 [[ -f $PILOT/COMPARISON_COMPLETE && ! -e $OUT/inputs/jobs.tsv ]]||{ echo 'Baseline incomplete or model scan already submitted' >&2;exit 1; }
 mkdir -p -- "$OUT/inputs/source/src" "$OUT/inputs/source/include" "$OUT/inputs/binaries" "$OUT/logs"
 sha256sum -c "$PILOT/inputs/source_checksums.txt" > "$OUT/inputs/frozen_base_verification.log"
 cp -- "$PILOT/inputs/digitizer_source/main.cxx" "$PILOT/inputs/digitizer_source/CMakeLists.txt" "$OUT/inputs/source/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.cxx "$OUT/inputs/source/src/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.h "$OUT/inputs/source/include/"
 perl -0pi -e 's/\A/#ifndef STUDY_NPIXEL10\n#define STUDY_NPIXEL10 10000\n#endif\n#ifndef STUDY_THRESHOLD_SHIFT\n#define STUDY_THRESHOLD_SHIFT 0.0\n#endif\n/;s/int nPixel10 = 10000;/int nPixel10 = STUDY_NPIXEL10;/;s/ApplyThreshold\(hADC, ex.thresholdValue\[layerNo\]\[chipNo\]\[chanNo\],/ApplyThreshold(hADC, ex.thresholdValue[layerNo][chipNo][chanNo] + STUDY_THRESHOLD_SHIFT,/' "$OUT/inputs/source/src/MCDigi.cxx"
 grep -q 'nPixel10 = STUDY_NPIXEL10' "$OUT/inputs/source/src/MCDigi.cxx"
 grep -q 'chanNo] + STUDY_THRESHOLD_SHIFT' "$OUT/inputs/source/src/MCDigi.cxx"
 cp -- "$HERE/digitization_model_scan.sh" "$OUT/inputs/"
 printf 'case\teffective_pixels_10um\tthreshold_shift_HG_ADC\n' > "$OUT/inputs/cases.tsv"
 printf 'task\tcase\tenergy_GeV\tfile_index\told_decode\told_calib\ttruth\tnew_decode\tnew_calib\n' > "$OUT/inputs/manifest.tsv"
 read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
 read -r -a ROOT_LIBS <<< "$(root-config --libs)"
 task=0
 for case_name in pixels14400 pixels19600 threshold_plus5 threshold_plus10;do
  pixels=10000;shift_adc=0
  case $case_name in pixels14400)pixels=14400;;pixels19600)pixels=19600;;threshold_plus5)shift_adc=5;;threshold_plus10)shift_adc=10;;esac
  printf '%s\t%s\t%s\n' "$case_name" "$pixels" "$shift_adc" >> "$OUT/inputs/cases.tsv"
  g++ -std=c++20 -O2 "${ROOT_CFLAGS[@]}" -DSTUDY_NPIXEL10="$pixels" -DSTUDY_THRESHOLD_SHIFT="$shift_adc" -I"$OUT/inputs/source/include" "$OUT/inputs/source/main.cxx" "$OUT/inputs/source/src/"{MCDigi,Extract}.cxx "${ROOT_LIBS[@]}" -o "$OUT/inputs/binaries/MCAnalysis_$case_name" > "$OUT/logs/build_$case_name.log" 2>&1
  while IFS=$'\t' read -r index energy file_index old_raw old_cal truth raw cal;do
   [[ $index == task || $energy != 100 ]]&&continue
   name=${truth##*/}
   dest=$OUT/mc/$case_name/100GeV
   mkdir -p -- "$dest/decode" "$dest/calib"
   printf '%s\t%s\t100\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$case_name" "$file_index" "$old_raw" "$old_cal" "$truth" "$dest/decode/$name" "$dest/calib/$name" >> "$OUT/inputs/manifest.tsv"
   task=$((task+1))
  done < "$PILOT/inputs/manifest.tsv"
 done
 [[ $task == 8 ]]||{ echo 'Expected8modeltasks' >&2;exit 1; }
 sha256sum "$OUT/inputs/source/"*.* "$OUT/inputs/source/src/"* "$OUT/inputs/source/include/"* "$OUT/inputs/binaries/"* "$OUT/inputs/digitization_model_scan.sh" "$OUT/inputs/cases.tsv" "$PILOT/inputs/calib" "$PILOT/inputs/reviewer_parameters/"* > "$OUT/inputs/source_parameter_checksums.txt"
 printf 'Diagnostic100GeVmodels:10um effectivepixels120x120/140x140 represent recovery/nonlinearity proxy, notliteralphysicalpixelcount;threshold+5/+10ADCforbothSiPMtypes;CTand(1-p)normalizationretained;calibrationunchanged;RNGSetSeed(0)independentstochasticrealizations;2truthfilesx1000events.\n' > "$OUT/inputs/model_definitions.txt"
 printf 'bash %q %q submit\n' "$HERE/digitization_model_scan.sh" "$LABEL" > "$OUT/inputs/command.sh"
 job=$(sbatch --parsable --partition=all --array=0-7 --mem=4G --time=02:00:00 --cpus-per-task=1 --ntasks=1 --job-name=ct_model_scan --output="$OUT/logs/slurm-%A_%a.out" --error="$OUT/logs/slurm-%A_%a.err" --wrap="bash '$HERE/digitization_model_scan.sh' '$LABEL' worker")
 printf '%s\t%s\n' "$(date -Is)" "$job" > "$OUT/inputs/jobs.tsv"
 printf 'Submitted bounded100GeV diagnostic model array %s\n' "$job"
elif [[ $ACTION == worker ]];then
 index=${SLURM_ARRAY_TASK_ID:?Array task required}
 line=$(sed -n "$((index+2))p" "$OUT/inputs/manifest.tsv")
 IFS=$'\t' read -r task case_name energy file_index old_raw old_cal truth raw cal <<< "$line"
 [[ $task == "$index" && $energy == 100 && ! -e $raw && ! -e $cal ]]||{ echo 'Bad task or existing output' >&2;exit 1; }
 PARAM=$PILOT/inputs/reviewer_parameters
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$old_raw" "$old_cal" > "$OUT/logs/task_${task}_inputs_before.tsv"
 "$OUT/inputs/binaries/MCAnalysis_$case_name" "$truth" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$PARAM/20230504_1425_Calib_LED_all.root" "$PARAM/threshold.root" "$raw" 0
 "$PILOT/inputs/calib" "$raw" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$cal"
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$old_raw" "$old_cal" > "$OUT/logs/task_${task}_inputs_after.tsv"
 cmp -s "$OUT/logs/task_${task}_inputs_before.tsv" "$OUT/logs/task_${task}_inputs_after.tsv"
 touch "$OUT/mc/$case_name/100GeV/${raw##*/}.COMPLETE"
else
 echo 'Usage: digitization_model_scan.sh LABEL submit|worker' >&2;exit 1
fi
