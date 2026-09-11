#!/usr/bin/env bash
# Windows diagnose for the "file 15 then file 22" hard abort.
#
# Two stages:
#   1. Minimal testthat repro (control: file 22 alone; repro: files 15,22 same
#      process) -- confirms the bug is still present in this build.
#   2. Parametrized matrix (.github/scripts/repro-15-22-matrix.R) that varies
#      team size, GC, libomp re-warm, and hold_memory to isolate the trigger.
#      Each variant runs in its OWN Rscript process so a crash in one does not
#      poison the others.
#
# Exit code is the matrix baseline variant's exit code (so the workflow fails
# when the bug is still present), but every variant's result is printed in a
# summary table regardless.
set -uo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"

echo "===== sessionInfo ====="
Rscript -e 'print(R.version.string); print(utils::packageVersion("adlaplace")); print(utils::packageVersion("RCppAD")); print(sessionInfo())' || true

# Alphabetical indices: 15=test-cppad-teardown.R, 22=test-fun-obj-fdfh.R
echo "::group::control: file 22 alone"
echo "+ Rscript run-adlaplace-testthat-suite.R --only=22"
set +e
Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R" --only=22
ctrl_ec=$?
set -e
echo "exit_code=${ctrl_ec}"
echo "::endgroup::"

echo "::group::repro: files 15 then 22 (same process)"
echo "+ Rscript run-adlaplace-testthat-suite.R --only=15,22"
set +e
Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R" --only=15,22
repro_ec=$?
set -e
echo "exit_code=${repro_ec}"
echo "::endgroup::"

echo "===== testthat diagnose: control=${ctrl_ec} repro_15_then_22=${repro_ec} ====="

# ---- matrix ---------------------------------------------------------------
# Each line: NAME|ENV0=VAL ENV1=VAL ...|ADLAPLACE_REPRO_*
# The script runs Rscript repro-15-22-matrix.R with the given env overrides.
# exit_code 0 = ok, 1 = bad values (caught), 127 = hard abort (the bug).

repro="${root}/.github/scripts/repro-15-22-matrix.R"

run_variant() {
  name="$1"; shift
  env_overrides="$1"; shift
  echo "::group::matrix: ${name}"
  echo "+ ${env_overrides} Rscript repro-15-22-matrix.R"
  set +e
  env ${env_overrides} Rscript "${repro}" 2>&1
  ec=$?
  set -e
  echo "exit_code=${ec}"
  echo "::endgroup::"
  printf '%s|%s|%d\n' "${name}" "${env_overrides}" "${ec}" >> "${root}/.matrix-results.txt"
}

# Fresh results file (relative to repo root so the workflow can upload it).
: > "${root}/.matrix-results.txt"

# BASELINE: A=4 threads/20 iters, B=2 threads. Expect 127 (the bug).
run_variant "BASELINE_A4_B2" \
  "ADLAPLACE_REPRO_PHASE_A_THREADS=4 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=none"

# SKIP_A: control -- phase B only. Expect 0.
run_variant "SKIP_A_B2" \
  "ADLAPLACE_REPRO_SKIP_A=1 ADLAPLACE_REPRO_PHASE_B_THREADS=2"

# SAME_SIZE: A=2 (not 4), B=2. Does the team-size CHANGE (4->2) matter?
run_variant "SAME_SIZE_A2_B2" \
  "ADLAPLACE_REPRO_PHASE_A_THREADS=2 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=none"

# INTERLUDE_GC: A=4 then gc() then B=2. Does finalizing phase A's ad_pack fix it?
run_variant "INTERLUDE_GC_A4_B2" \
  "ADLAPLACE_REPRO_PHASE_A_THREADS=4 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=gc"

# INTERLUDE_WARM: A=4 then warm_openmp_runtime() then B=2. libomp TLS re-init?
run_variant "INTERLUDE_WARM_A4_B2" \
  "ADLAPLACE_REPRO_PHASE_A_THREADS=4 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=warm"

# INTERLUDE_GCWARM: A=4 then gc()+warm then B=2. Both mitigations.
run_variant "INTERLUDE_GCWARM_A4_B2" \
  "ADLAPLACE_REPRO_PHASE_A_THREADS=4 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=gcwarm"

# HOLD0: A=4 then B=2 with ADLAPLACE_HOLD_MEMORY=0. Is hold_memory the trigger?
# (env override placed before the ADLAPLACE_REPRO_* vars so the C++ side sees it.)
run_variant "HOLD0_A4_B2" \
  "ADLAPLACE_HOLD_MEMORY=0 ADLAPLACE_REPRO_PHASE_A_THREADS=4 ADLAPLACE_REPRO_PHASE_A_ITERS=20 ADLAPLACE_REPRO_PHASE_B_THREADS=2 ADLAPLACE_REPRO_INTERLUDE=none"

echo "===== matrix summary ====="
echo "variant|env|exit_code"
cat "${root}/.matrix-results.txt"
echo "===== matrix done ====="

# Exit on the baseline variant's code so the workflow fails while the bug
# is still present, but the summary above shows every variant's result.
baseline_ec=$(grep '^BASELINE_A4_B2|' "${root}/.matrix-results.txt" | cut -d'|' -f3)
baseline_ec=${baseline_ec:-0}
exit "${baseline_ec}"
