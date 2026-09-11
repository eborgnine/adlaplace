#!/usr/bin/env bash
# Extra Windows diagnostics for the early adlaplace abort.
# Continues across steps so one hard abort does not hide later signals.
set -uo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"

run() {
  local title=$1
  shift
  echo "::group::${title}"
  echo "+ $*"
  "$@"
  local ec=$?
  echo "exit_code=${ec}"
  echo "::endgroup::"
  return 0
}

echo "===== sessionInfo (installed adlaplace / RCppAD) ====="
Rscript -e 'print(R.version.string); print(utils::packageVersion("adlaplace")); print(utils::packageVersion("RCppAD")); print(sessionInfo())' || true

run "repro ad_pack_ptr breadcrumbs" \
  Rscript "${root}/.github/scripts/repro-ad-pack-ptr.R"

# Early files that produce the five gamSim prints before the R-universe abort.
for f in \
  adlaplace/tests/testthat/test-collect-terms-constructors.R \
  adlaplace/tests/testthat/test-format-parameters.R \
  adlaplace/tests/testthat/test-fun-obj-fdfh.R
do
  run "test_file $(basename "$f")" \
    Rscript "${root}/.github/scripts/run-adlaplace-test-files.R" "$f"
done

echo "===== diagnose-adlaplace-windows.sh finished ====="
