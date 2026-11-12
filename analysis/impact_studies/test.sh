#!/bin/bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
test_dir=$(mktemp -d /tmp/impact-studies-test.XXXXXX)
trap 'rm -rf -- "$test_dir"' EXIT
"$HERE/TestHLIntercept" fixtures "$test_dir/input"
args=(--energy 1 --data-decode "$test_dir/input/raw.root" --data-calib "$test_dir/input/cal.root" --mc-decode "$test_dir/input/raw.root" --mc-calib "$test_dir/input/cal.root" --hl "$test_dir/input/hl.root" --mip "$test_dir/input/mip.root" --pedestal "$test_dir/input/ped.root" --result-dir "$test_dir/result")
bash "$HERE/run.sh" "${args[@]}" --dry-run
[[ ! -e $test_dir/result ]]
bash "$HERE/run.sh" "${args[@]}"
"$HERE/TestHLIntercept" check "$test_dir/result/hl_intercept/threshold/default/1GeV"
"$HERE/TestHLIntercept" binning "$test_dir"
study=$test_dir/result/hl_intercept/threshold/default/1GeV
cp "$study/figures/event_energy.png" "$test_dir/original.png"
cp "$study/summary.tsv" "$test_dir/original.tsv"
bash "$HERE/replot.sh" --energy 1 --result-dir "$test_dir/result"
cmp "$study/figures/event_energy.png" "$test_dir/original.png"
cmp "$study/summary.tsv" "$test_dir/original.tsv"
# Legacy outputs have no data-event tree/cache; exercise calibrated-data fallback.
"$HERE/TestHLIntercept" legacy "$study"
rm "$study/event_histograms.root"
bash "$HERE/replot.sh" --energy all --result-dir "$test_dir/result"
cmp "$study/figures/event_energy.png" "$test_dir/original.png"
cmp "$study/summary.tsv" "$test_dir/original.tsv"
# Further redraws can use the new data-event cache without reopening Calib_Hit.
mv "$test_dir/input/cal.root" "$test_dir/input/cal.hidden"
bash "$HERE/replot.sh" --energy 1 --result-dir "$test_dir/result"
mv "$test_dir/input/cal.hidden" "$test_dir/input/cal.root"
cmp "$study/figures/event_energy.png" "$test_dir/original.png"
bash "$HERE/summarize.sh" --result-dir "$test_dir/result"
if bash "$HERE/run.sh" "${args[@]}" > "$test_dir/repeat.log" 2>&1;then echo 'Overwrite was accepted' >&2;exit 1;fi
"$HERE/TestHLIntercept" bad "$test_dir/input"
if bash "$HERE/run.sh" "${args[@]}" --label bad > "$test_dir/bad.log" 2>&1;then echo 'Incorrect baseline was accepted' >&2;exit 1;fi
rg -q 'Baseline calibration mismatch' "$test_dir/bad.log"
[[ ! -e $test_dir/result/hl_intercept/threshold/bad/1GeV/study.root ]]
# Exercise the actual submission wrapper using a fake sbatch; no cluster jobs.
mkdir -p "$test_dir/bin"
cat > "$test_dir/bin/sbatch" <<'EOF'
#!/bin/bash
printf '%s\n' "$*" >> "$IMPACT_TEST_CALLS"
printf '12345\n'
EOF
chmod +x "$test_dir/bin/sbatch"
export IMPACT_TEST_CALLS=$test_dir/submissions.txt
PATH="$test_dir/bin:$PATH" bash "$HERE/submit.sh" "${args[@]}" --label submitted --log-dir "$test_dir/jobs" --no-plots
[[ $(wc -l < "$IMPACT_TEST_CALLS") -eq 1 ]]
rg -q -- '--energy 1' "$IMPACT_TEST_CALLS"
rg -q -- '--label submitted' "$IMPACT_TEST_CALLS"
echo 'All impact study tests passed'
