#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
label=${1:-20261013}
action=${2:-submit}
[[ $label =~ ^[a-zA-Z0-9_-]+$ ]] || { echo 'Invalid label' >&2; exit 1; }
out="$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$label/digitization"
mc=/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new
calib=/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO
params="$mc/Analysis_edit/share"
thre="$PROJECT/analysis/result/threshold/threshold.root"
if [[ $action == submit ]];then
 [[ ! -e $out/inputs/jobs.tsv ]] || { echo "Study already submitted: $out" >&2; exit 1; }
 mkdir -p "$out/inputs/digitizer_source" "$out/inputs/calibration_source"
 cp -- "$mc/Analysis_edit/build/MCAnalysis" "$out/inputs/MCAnalysis"
 cp -- "$calib/bin/calib" "$out/inputs/calib"
 cp -- "$mc/Analysis_edit/"{main.cxx,CMakeLists.txt} "$out/inputs/digitizer_source/"
 cp -- "$mc/Analysis_edit/src/"{MCDigi,Extract}.cxx "$out/inputs/digitizer_source/"
 cp -- "$mc/Analysis_edit/include/"{MCDigi,Extract}.h "$out/inputs/digitizer_source/"
 cp -- "$calib/"calib.cxx "$calib/src/"{Calibration,EBUdecode}.cxx "$out/inputs/calibration_source/"
 cp -- "$HERE/ct_digitization_study.sh" "$out/inputs/"
 sha256sum "$out/inputs/"{MCAnalysis,calib} "$out/inputs/digitizer_source/"* "$out/inputs/calibration_source/"* "$params/"{pedestal2023_SPS,all_hl_electron2023_hlratio_v3,all_auto_muon_v4_trackfit,20230504_1425_Calib_LED_all}.root "$thre" > "$out/inputs/source_checksums.txt"
 printf 'task\tenergy_GeV\tfile_index\told_decode\told_calib\ttruth\tnew_decode\tnew_calib\n' > "$out/inputs/manifest.tsv"
 task=0
 for energy in 5 20 100;do
  for file_index in 0 5;do
   line=$(sed -n "$((file_index+1))p" "$PROJECT/investigations/result/e_ssa_threshold_review_20260930/$label/${energy}GeV/baseline/inputs/mc.tsv")
   IFS=$'\t' read -r old_raw old_cal <<< "$line"
   name=${old_raw##*/}
   truth="$mc/Result_MC/generate/e-/e_ssa/${energy}GeV/$name"
   [[ -f $truth && -f $old_raw && -f $old_cal ]] || { echo "Missing study input: $truth" >&2; exit 1; }
   dest="$out/mc/ct_input/${energy}GeV"
   mkdir -p "$dest/"{decode,calib}
   printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$energy" "$file_index" "$old_raw" "$old_cal" "$truth" "$dest/decode/$name" "$dest/calib/$name" >> "$out/inputs/manifest.tsv"
   stat -L --printf='%n\t%s\t%Y\n' "$old_raw" "$old_cal" "$truth" >> "$out/inputs/input_stats_before.tsv"
   task=$((task+1))
  done
 done
 printf 'ROOT digitizer calls SetSeed(0); generation is stochastic. Frozen binaries/sources and generated ROOTs are retained for exact analysis reproduction. Calibration MIP is unchanged. Two existing truth files (indices 0 and 5) per energy; 5/20/100 GeV only.\n' > "$out/inputs/generation_conditions.txt"
 job=$(sbatch --parsable --partition=all --array=0-5 --mem=4G --time=02:00:00 --cpus-per-task=1 --ntasks=1 --job-name=ct_digitization --output="$out/slurm-%A_%a.out" --error="$out/slurm-%A_%a.err" --wrap="bash '$HERE/ct_digitization_study.sh' '$label' worker")
 printf '%s\t%s\n' "$(date -Is)" "$job" > "$out/inputs/jobs.tsv"
 printf 'Submitted CT-normalized diagnostic digitization array %s\nOutput: %s\n' "$job" "$out"
elif [[ $action == worker ]];then
 index=${SLURM_ARRAY_TASK_ID:?Array task ID required}
 line=$(sed -n "$((index+2))p" "$out/inputs/manifest.tsv")
 IFS=$'\t' read -r task energy file_index old_raw old_cal truth raw cal <<< "$line"
 [[ $task == "$index" && -n $cal ]] || { echo 'Bad task row' >&2; exit 1; }
 [[ ! -e $raw && ! -e $cal ]] || { echo 'Diagnostic output already exists' >&2; exit 1; }
 dest=${raw%/decode/*}
 exec > >(tee "$dest/${raw##*/}.run.log") 2>&1
 date -Is
 stat -L --printf='%n\t%s\t%Y\n' "$truth" "$old_raw" "$old_cal" > "$dest/${raw##*/}.input_before.tsv"
 "$out/inputs/MCAnalysis" "$truth" "$params/pedestal2023_SPS.root" "$params/all_hl_electron2023_hlratio_v3.root" "$params/all_auto_muon_v4_trackfit.root" "$params/20230504_1425_Calib_LED_all.root" "$thre" "$raw" 0
 "$out/inputs/calib" "$raw" "$params/pedestal2023_SPS.root" "$params/all_hl_electron2023_hlratio_v3.root" "$params/all_auto_muon_v4_trackfit.root" "$cal"
 stat -L --printf='%n\t%s\t%Y\n' "$truth" "$old_raw" "$old_cal" > "$dest/${raw##*/}.input_after.tsv"
 cmp "$dest/${raw##*/}.input_before.tsv" "$dest/${raw##*/}.input_after.tsv"
 touch "$dest/${raw##*/}.COMPLETE"
 date -Is
else
 echo 'Usage: ct_digitization_study.sh LABEL submit|worker' >&2;exit 1
fi
