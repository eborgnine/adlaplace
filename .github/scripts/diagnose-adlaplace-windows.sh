#!/usr/bin/env bash
# Minimal Windows repro: file 15 (cppad-teardown) then 22 (fun-obj-fdfh)
# in one process. File 22 alone is the control (should pass).
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

echo "===== diagnose finished: control=${ctrl_ec} repro_15_then_22=${repro_ec} ====="
exit "${repro_ec}"
