#!/bin/bash
# Shell orchestration tests only; no ROOT processing or real cluster submission.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
test_dir=$(mktemp -d /tmp/impact-all-test.XXXXXX)
trap 'rm -rf -- "$test_dir"' EXIT
mkdir -p "$test_dir/driver" "$test_dir/bin"
cp "$HERE/energies.sh" "$HERE/run_all.sh" "$HERE/submit.sh" "$test_dir/driver/"
export IMPACT_TEST_LOG=$test_dir/calls
cat > "$test_dir/driver/run.sh" <<'EOF'
#!/bin/bash
set -euo pipefail
dry=0 energy=''
while (($#));do case $1 in --dry-run) dry=1;shift;; --energy) energy=$2;shift 2;; *) shift;; esac;done
echo "run $energy $dry" >> "$IMPACT_TEST_LOG"
if [[ ${IMPACT_TEST_FAIL:-} == "$energy" && $dry == 0 ]];then exit 1;fi
EOF
cat > "$test_dir/driver/summarize.sh" <<'EOF'
#!/bin/bash
echo "summary $*" >> "$IMPACT_TEST_LOG"
EOF
cat > "$test_dir/bin/sbatch" <<'EOF'
#!/bin/bash
echo "sbatch $*" >> "$IMPACT_TEST_LOG"
count=$(rg -c '^sbatch ' "$IMPACT_TEST_LOG")
echo "$((1000+count));test_cluster"
EOF
cat > "$test_dir/bin/make" <<'EOF'
#!/bin/bash
exit 0
EOF
chmod +x "$test_dir/bin/"*
export PATH="$test_dir/bin:$PATH"
source "$HERE/energies.sh"
bash "$test_dir/driver/submit.sh" --energy all --summarize --label test --log-dir "$test_dir/jobs" --dry-run > "$test_dir/dry.log"
[[ $(rg -c '^run ' "$IMPACT_TEST_LOG") -eq 19 ]]
! rg -q '^sbatch ' "$IMPACT_TEST_LOG"
: > "$IMPACT_TEST_LOG"
bash "$test_dir/driver/submit.sh" --energy all --summarize --label test --log-dir "$test_dir/jobs" > "$test_dir/batch.log"
[[ $(rg -c '^sbatch ' "$IMPACT_TEST_LOG") -eq 20 ]]
dependency=$(seq -s : 1001 1019)
tail -n 1 "$IMPACT_TEST_LOG" | rg -q -- "--dependency=afterok:$dependency"
tail -n 1 "$IMPACT_TEST_LOG" | rg -q -- "--energies $all_energies"
tail -n 1 "$IMPACT_TEST_LOG" | rg -q -- '--label test'
: > "$IMPACT_TEST_LOG"
bash "$test_dir/driver/run_all.sh" --label test > "$test_dir/direct.log"
[[ $(rg -c '^run .* 1$' "$IMPACT_TEST_LOG") -eq 19 ]]
[[ $(rg -c '^run .* 0$' "$IMPACT_TEST_LOG") -eq 19 ]]
tail -n 1 "$IMPACT_TEST_LOG" | rg -q "^summary .*--energies $all_energies"
: > "$IMPACT_TEST_LOG"
if IMPACT_TEST_FAIL=5 bash "$test_dir/driver/run_all.sh" > "$test_dir/fail.log";then echo 'Failed analysis did not stop the scan' >&2;exit 1;fi
! rg -q '^summary ' "$IMPACT_TEST_LOG"
! rg -q '^run 10 0$' "$IMPACT_TEST_LOG"
# The real summary driver must reject an incomplete exact list before invoking C++.
cp "$HERE/summarize.sh" "$test_dir/driver/summarize.sh"
cat > "$test_dir/driver/SummarizeImpact" <<'EOF'
#!/bin/bash
echo "cpp_summary $*" >> "$IMPACT_TEST_LOG"
EOF
chmod +x "$test_dir/driver/SummarizeImpact"
base=$test_dir/results/hl_intercept/threshold/test
mkdir -p "$base/1GeV" "$base/0.5GeV"
touch "$base/1GeV/study.root" "$base/0.5GeV/study.root"
: > "$IMPACT_TEST_LOG"
if bash "$test_dir/driver/summarize.sh" --result-dir "$test_dir/results" --label test --energies 1,2 > "$test_dir/incomplete.log" 2>&1;then exit 1;fi
[[ ! -s $IMPACT_TEST_LOG ]]
bash "$test_dir/driver/summarize.sh" --result-dir "$test_dir/results" --label test --energies 1
rg -q '/1GeV/study.root' "$IMPACT_TEST_LOG"
! rg -q '/0.5GeV/study.root' "$IMPACT_TEST_LOG"
echo 'All-energy shell tests passed'
