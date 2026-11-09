#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
LABEL=${1:-20261013}; ACTION=${2:-submit}
[[ $LABEL =~ ^[A-Za-z0-9_-]+$ ]]||{ echo 'Invalid label' >&2;exit 1; }
BASE=$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$LABEL
PILOT=$BASE/digitization; OUT=$BASE/joint_response/digitization
if [[ $ACTION == submit ]];then
 [[ -f $PILOT/COMPARISON_COMPLETE && ! -e $OUT/inputs/jobs.tsv ]]||{ echo 'Baseline incomplete or joint study already submitted' >&2;exit 1; }
 mkdir -p -- "$OUT/inputs/source/src" "$OUT/inputs/source/include" "$OUT/inputs/binaries" "$OUT/logs"
 sha256sum -c "$PILOT/inputs/source_checksums.txt" > "$OUT/inputs/frozen_base_verification.log"
 cp -- "$PILOT/inputs/digitizer_source/main.cxx" "$PILOT/inputs/digitizer_source/CMakeLists.txt" "$OUT/inputs/source/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.cxx "$OUT/inputs/source/src/"
 cp -- "$PILOT/inputs/digitizer_source/"{MCDigi,Extract}.h "$OUT/inputs/source/include/"
 perl -0pi -e 's/\A/#ifndef STUDY_NPIXEL10\n#define STUDY_NPIXEL10 10000\n#endif\n#ifndef STUDY_THRESHOLD_SHIFT\n#define STUDY_THRESHOLD_SHIFT 0.0\n#endif\n#ifndef STUDY_YIELD10\n#define STUDY_YIELD10 1.10\n#endif\n#ifndef STUDY_YIELD15\n#define STUDY_YIELD15 1.05\n#endif\n/;s/int nPixel10 = 10000;/int nPixel10 = STUDY_NPIXEL10;/;s/double lightYieldperMIP = \(1-pct_SiPMType\) \*/double lightYieldperMIP = (SiPMType == 0 ? STUDY_YIELD10 : STUDY_YIELD15) * (1-pct_SiPMType) */;s/ApplyThreshold\(hADC, ex.thresholdValue\[layerNo\]\[chipNo\]\[chanNo\],/ApplyThreshold(hADC, ex.thresholdValue[layerNo][chipNo][chanNo] + STUDY_THRESHOLD_SHIFT,/' "$OUT/inputs/source/src/MCDigi.cxx"
 rg -q 'lightYieldperMIP = .*STUDY_YIELD10' "$OUT/inputs/source/src/MCDigi.cxx"
 cp -- "$HERE/joint_digitization_study.sh" "$OUT/inputs/"
 printf 'case\tyield_scale_10um\tyield_scale_15um\teffective_pixels_10um\tthreshold_shift_HG_ADC\tbeam_energies_GeV\n' > "$OUT/inputs/cases.tsv"
 printf 'task\tcase\tenergy_GeV\tfile_index\told_decode\told_calib\ttruth\tnew_decode\tnew_calib\n' > "$OUT/inputs/manifest.tsv"
 read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"; read -r -a ROOT_LIBS <<< "$(root-config --libs)"
 task=0
 for case_name in muon_norm norm_threshold10 norm_pixels19600 norm_pixels19600_threshold10;do
  pixels=10000; shift_adc=0; energies=100
  case $case_name in muon_norm)energies=5,20,100;;norm_threshold10)shift_adc=10;;norm_pixels19600)pixels=19600;;norm_pixels19600_threshold10)pixels=19600;shift_adc=10;;esac
  printf '%s\t1.10\t1.05\t%s\t%s\t%s\n' "$case_name" "$pixels" "$shift_adc" "$energies" >> "$OUT/inputs/cases.tsv"
  g++ -std=c++20 -O2 "${ROOT_CFLAGS[@]}" -DSTUDY_NPIXEL10="$pixels" -DSTUDY_THRESHOLD_SHIFT="$shift_adc" -I"$OUT/inputs/source/include" "$OUT/inputs/source/main.cxx" "$OUT/inputs/source/src/"{MCDigi,Extract}.cxx "${ROOT_LIBS[@]}" -o "$OUT/inputs/binaries/MCAnalysis_$case_name" > "$OUT/logs/build_$case_name.log" 2>&1
  while IFS=$'\t' read -r index energy file_index old_raw old_cal truth raw cal;do
   [[ $index == task ]]&&continue
   [[ $case_name != muon_norm && $energy != 100 ]]&&continue
   name=${truth##*/}; dest=$OUT/mc/$case_name/${energy}GeV
   mkdir -p -- "$dest/decode" "$dest/calib"
   printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$case_name" "$energy" "$file_index" "$old_raw" "$old_cal" "$truth" "$dest/decode/$name" "$dest/calib/$name" >> "$OUT/inputs/manifest.tsv"
   task=$((task+1))
  done < "$PILOT/inputs/manifest.tsv"
 done
 [[ $task == 12 ]]||{ echo 'Expected 12 bounded joint tasks' >&2;exit 1; }
 sha256sum "$OUT/inputs/source/"*.* "$OUT/inputs/source/src/"* "$OUT/inputs/source/include/"* "$OUT/inputs/binaries/"* "$OUT/inputs/joint_digitization_study.sh" "$OUT/inputs/cases.tsv" "$PILOT/inputs/calib" "$PILOT/inputs/reviewer_parameters/"* > "$OUT/inputs/source_parameter_checksums.txt"
 printf 'Frozen CT-normalized digitizer,20C; primary photon yield multiplied by 1.10/1.05 for10/15um. Multipliers informed by completed muon MPV ratios,not exact validated muon closure. Shared calibration MIP,gain,CT,noise unchanged. Effective19600 pixels is recovery/nonlinearity proxy;threshold+10ADC is diagnostic. Independent ROOT SetSeed(0) draws;same two truthfilesx1000events;norm5/20/100GeV,others100GeVonly.\n' > "$OUT/inputs/model_definitions.txt"
 printf 'bash %q %q submit\n' "$HERE/joint_digitization_study.sh" "$LABEL" > "$OUT/inputs/command.sh"
 job=$(sbatch --parsable --partition=all --array=0-11 --mem=4G --time=02:00:00 --cpus-per-task=1 --ntasks=1 --job-name=joint_response --output="$OUT/logs/slurm-%A_%a.out" --error="$OUT/logs/slurm-%A_%a.err" --wrap="bash '$HERE/joint_digitization_study.sh' '$LABEL' worker")
 printf '%s\t%s\n' "$(date -Is)" "$job" > "$OUT/inputs/jobs.tsv"
 printf 'Submitted joint ADC/hit diagnostic array %s\n' "$job"
elif [[ $ACTION == worker ]];then
 index=${SLURM_ARRAY_TASK_ID:?Array task required}; line=$(sed -n "$((index+2))p" "$OUT/inputs/manifest.tsv")
 IFS=$'\t' read -r task case_name energy file_index old_raw old_cal truth raw cal <<< "$line"
 [[ $task == "$index" && ! -e $raw && ! -e $cal ]]||{ echo 'Bad task or existing output' >&2;exit 1; }
 PARAM=$PILOT/inputs/reviewer_parameters
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$old_raw" "$old_cal" > "$OUT/logs/task_${task}_inputs_before.tsv"
 "$OUT/inputs/binaries/MCAnalysis_$case_name" "$truth" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$PARAM/20230504_1425_Calib_LED_all.root" "$PARAM/threshold.root" "$raw" 0
 "$PILOT/inputs/calib" "$raw" "$PARAM/pedestal2023_SPS.root" "$PARAM/all_hl_electron2023_hlratio_v3.root" "$PARAM/all_auto_muon_v4_trackfit.root" "$cal"
 stat -L --printf='%n\t%s\t%Y\n' -- "$truth" "$old_raw" "$old_cal" > "$OUT/logs/task_${task}_inputs_after.tsv"
 cmp -s "$OUT/logs/task_${task}_inputs_before.tsv" "$OUT/logs/task_${task}_inputs_after.tsv"
 touch "$OUT/mc/$case_name/${energy}GeV/${raw##*/}.COMPLETE"
else echo 'Usage: joint_digitization_study.sh LABEL submit|worker' >&2;exit 1;fi
