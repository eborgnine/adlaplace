#!/usr/bin/env bash
# Windows diagnose for R-CMD-check test failures after the team-size latch.
#
# Main CI now gets past the old hard abort (15 then 22) but fails expectations
# in test-log-lik-deriv-parallel.R and test-reorder-shards.R that assert
#   owners == (0:(n-1)) %% 2
# after earlier tests have raised the process high-water mark to 4.
#
# Stages:
#   1. Regression: file 22 alone + files 15 then 22 (abort should stay gone).
#   2. testthat isolation of the new failures (files 28 / 42 alone vs after 15).
#   3. Minimal owner-modulo matrix (warmup 4 then probe 2), with/without clamp.
#
# Alphabetical indices:
#   15 = test-cppad-teardown.R   (4-thread warmup in suite)
#   22 = test-fun-obj-fdfh.R
#   28 = test-log-lik-deriv-parallel.R
#   42 = test-reorder-shards.R
#
# Exit code is non-zero if the latch-poisoned owner expectations still fail
# (AFTER4_* or PROBE_AFTER4), so the workflow stays red while the issue is open.
# The summary table always prints every variant.
set -uo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"

suite="${root}/.github/scripts/run-adlaplace-testthat-suite.R"
owner_repro="${root}/.github/scripts/repro-owner-modulo.R"
results="${root}/.matrix-results.txt"
: > "${results}"

echo "===== sessionInfo ====="
Rscript -e 'print(R.version.string); print(utils::packageVersion("adlaplace")); print(utils::packageVersion("RCppAD")); print(sessionInfo())' || true

VARIANT_TIMEOUT="${VARIANT_TIMEOUT:-180}"

run_with_timeout() {
  if command -v timeout >/dev/null 2>&1; then
    timeout -k 10 "${VARIANT_TIMEOUT}" "$@"
  else
    "$@"
  fi
}

record() {
  name="$1"
  env_label="$2"
  ec="$3"
  printf '%s|%s|%d\n' "${name}" "${env_label}" "${ec}" >> "${results}"
}

run_suite() {
  name="$1"
  only="$2"
  env_overrides="${3:-}"
  echo "::group::suite: ${name}"
  echo "+ ${env_overrides} Rscript run-adlaplace-testthat-suite.R --only=${only} (timeout ${VARIANT_TIMEOUT}s)"
  set +e
  if [ -n "${env_overrides}" ]; then
    run_with_timeout env ${env_overrides} Rscript "${suite}" --only="${only}" 2>&1
  else
    run_with_timeout Rscript "${suite}" --only="${only}" 2>&1
  fi
  ec=$?
  set -e
  echo "exit_code=${ec}"
  echo "::endgroup::"
  record "${name}" "only=${only} ${env_overrides}" "${ec}"
}

run_owner() {
  name="$1"
  env_overrides="$2"
  echo "::group::owner-modulo: ${name}"
  echo "+ ${env_overrides} Rscript repro-owner-modulo.R (timeout ${VARIANT_TIMEOUT}s)"
  set +e
  run_with_timeout env ${env_overrides} Rscript "${owner_repro}" 2>&1
  ec=$?
  set -e
  echo "exit_code=${ec}"
  echo "::endgroup::"
  record "${name}" "${env_overrides}" "${ec}"
}

# ---- 1. Old abort regression ---------------------------------------------
run_suite "CONTROL_22" "22"
run_suite "REPRO_15_THEN_22" "15,22"

# ---- 2. New expectation failures (testthat) ------------------------------
# Alone: should pass (no prior high-water mark in a fresh process).
run_suite "ALONE_28_LOG_LIK" "28"
run_suite "ALONE_42_REORDER" "42"

# After file 15 (4-thread work): with latch on, expect owner %% 2 failures.
run_suite "AFTER4_15_THEN_28" "15,28"
run_suite "AFTER4_15_THEN_42" "15,42"

# Same sequence with clamp disabled: owners should match %% 2, OR hard-abort.
run_suite "AFTER4_15_THEN_28_NOCLAMP" "15,28" "ADLAPLACE_CLAMP_TEAM_THREADS=0"
run_suite "AFTER4_15_THEN_42_NOCLAMP" "15,42" "ADLAPLACE_CLAMP_TEAM_THREADS=0"

# ---- 3. Minimal owner-modulo matrix --------------------------------------
run_owner "PROBE_ALONE_2" \
  "ADLAPLACE_REPRO_SKIP_WARMUP=1 ADLAPLACE_REPRO_PROBE_THREADS=2"

run_owner "PROBE_AFTER4" \
  "ADLAPLACE_REPRO_WARMUP_THREADS=4 ADLAPLACE_REPRO_PROBE_THREADS=2"

run_owner "PROBE_AFTER4_NOCLAMP" \
  "ADLAPLACE_CLAMP_TEAM_THREADS=0 ADLAPLACE_REPRO_WARMUP_THREADS=4 ADLAPLACE_REPRO_PROBE_THREADS=2"

run_owner "PROBE_SAME_SIZE_4" \
  "ADLAPLACE_REPRO_WARMUP_THREADS=4 ADLAPLACE_REPRO_PROBE_THREADS=4"

echo "===== diagnose summary ====="
echo "variant|env|exit_code"
cat "${results}"
echo "===== interpret ====="
echo "0 = ok; 1 = testthat/owner mismatch; 127 = hard abort"
echo "Expect ALONE_* and PROBE_ALONE_2 / PROBE_SAME_SIZE_4 = 0."
echo "If AFTER4_* / PROBE_AFTER4 = 1 with clamp on, latch is breaking %%2 owner tests."
echo "If AFTER4_*_NOCLAMP / PROBE_AFTER4_NOCLAMP = 127, old decrease-abort remains."
echo "===== diagnose done ====="

# Fail the job while the latch still breaks owner expectations after a
# larger parallel team (the open main-CI failure mode).
after28=$(grep '^AFTER4_15_THEN_28|' "${results}" | cut -d'|' -f3)
after42=$(grep '^AFTER4_15_THEN_42|' "${results}" | cut -d'|' -f3)
probe=$(grep '^PROBE_AFTER4|' "${results}" | cut -d'|' -f3)
after28=${after28:-0}
after42=${after42:-0}
probe=${probe:-0}

if [ "${after28}" != "0" ] || [ "${after42}" != "0" ] || [ "${probe}" != "0" ]; then
  exit 1
fi
exit 0
