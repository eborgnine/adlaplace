#!/usr/bin/env bash
# Windows-only diagnostics for adlaplace test_check hard aborts.
# 1) Quick isolated repro (known-good path smoke)
# 2) Same-process full suite with Location reporter (order/state crash)
# 3) On abort, binary-search prefixes in fresh processes
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
  return "$ec"
}

echo "===== sessionInfo ====="
Rscript -e 'print(R.version.string); print(utils::packageVersion("adlaplace")); print(utils::packageVersion("RCppAD")); print(sessionInfo())' || true

run "repro ad_pack_ptr breadcrumbs" \
  Rscript "${root}/.github/scripts/repro-ad-pack-ptr.R" || true

echo "::group::same-process full testthat suite"
echo "+ Rscript run-adlaplace-testthat-suite.R"
set +e
Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R"
suite_ec=$?
set -e
echo "exit_code=${suite_ec}"
echo "::endgroup::"

if [ "${suite_ec}" -eq 0 ]; then
  echo "===== diagnose finished: full suite OK ====="
  exit 0
fi

# Soft failures and hard aborts both non-zero. Bisect only helps if R actually
# ran tests (look for n_files in the group above; "No test files" is fixed).
echo "===== full suite failed/aborted (exit=${suite_ec}); binary-search prefixes ====="
n_files=$(find adlaplace/tests/testthat -name 'test-*.R' | wc -l | tr -d ' ')
lo=1
hi=$n_files
# Find smallest prefix 1:k that fails (last END FILE before abort ≈ culprit zone).
while [ "$lo" -lt "$hi" ]; do
  mid=$(( (lo + hi) / 2 ))
  echo "::group::prefix 1:${mid}/${n_files}"
  set +e
  Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R" --from=1 --to="$mid"
  ec=$?
  set -e
  echo "exit_code=${ec}"
  echo "::endgroup::"
  if [ "$ec" -eq 0 ]; then
    lo=$(( mid + 1 ))
  else
    hi=$mid
  fi
done

echo "===== first failing prefix ends at file index ${lo} (of ${n_files}) ====="
echo "::group::re-run failing prefix 1:${lo}"
set +e
Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R" --from=1 --to="$lo"
ec=$?
set -e
echo "exit_code=${ec}"
echo "::endgroup::"

if [ "$lo" -gt 1 ]; then
  echo "::group::re-run single file at index ${lo}"
  set +e
  Rscript "${root}/.github/scripts/run-adlaplace-testthat-suite.R" --from="$lo" --to="$lo"
  ec=$?
  set -e
  echo "exit_code=${ec}"
  echo "::endgroup::"
fi

exit 1
