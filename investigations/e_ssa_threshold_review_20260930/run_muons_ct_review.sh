#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
PROJECT=$(cd -- "$HERE/../.." && pwd -P)
OUT=${1:-"$PROJECT/investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct"}
[[ $# -le 1 ]]||{ echo 'Usage: bash run_muons_ct_review.sh [OUTPUT_DIR]' >&2;exit 1; }
OUT=$(realpath -m -- "$OUT")
mkdir -p -- "$OUT"
for tag in threshold threshold_MPVwoCT;do
 if [[ $tag == threshold ]];then label=old_threshold;else label=new_threshold_MPVwoCT;fi
 collection=$OUT/$label
 if [[ -f $collection/COMPLETE && -f $collection/flat.root ]];then
  echo "Use completed collection: $collection"
  [[ -f $collection/inputs/input_stats_before.tsv && -f $collection/inputs/input_stats_after.tsv ]]
  cmp -s -- "$collection/inputs/input_stats_before.tsv" "$collection/inputs/input_stats_after.tsv"
  while IFS=$'\t' read -r input_path input_size input_mtime;do
   [[ $(stat -L --printf='%s' -- "$input_path") == "$input_size" && $(stat -L --printf='%Y' -- "$input_path") == "$input_mtime" ]]||{ echo "Input changed since collection: $input_path" >&2;exit 1; }
  done < "$collection/inputs/input_stats_after.tsv"
 else
  bash "$PROJECT/analysis/channel_response/muon_mip_response/run.sh" --runs 44 --max-events 20000 --mc-tag "$tag" --collect-only --label "$label" --result-dir "$OUT"
 fi
done
mkdir -p -- "$PROJECT/.scratch"
TASK_TMP=$(mktemp -d "$PROJECT/.scratch/muons_ct_review.XXXXXXXX")
trap 'rm -rf -- "$TASK_TMP"' EXIT
read -r -a ROOT_CFLAGS <<< "$(root-config --cflags)"
read -r -a ROOT_LIBS <<< "$(root-config --libs)"
g++ -std=c++17 -O2 -Wall -Wextra "${ROOT_CFLAGS[@]}" "$HERE/MuonsCTReview.cc" "${ROOT_LIBS[@]}" -o "$TASK_TMP/MuonsCTReview" > "$OUT/build.log" 2>&1
sha256sum "$HERE/MuonsCTReview.cc" "$HERE/run_muons_ct_review.sh" "$PROJECT/analysis/channel_response/muon_mip_response/CurrentCalibration.hh" "$PROJECT/analysis/channel_response/muon_mip_response/MuonMIPCommon.hh" > "$OUT/reviewer_code_checksums.txt"
printf 'bash %q %q\n' "$HERE/run_muons_ct_review.sh" "$OUT" > "$OUT/command.sh"
"$TASK_TMP/MuonsCTReview" "$OUT/old_threshold" "$OUT/new_threshold_MPVwoCT" "$OUT" > "$OUT/review.log" 2>&1
printf 'Completed focused muon comparison: %s\n' "$OUT"
