#!/usr/bin/env bash
# Windows vignette timing: knit adlaplace vignettes at 1 vs 4 OpenMP threads.
#
# Requires a local OpenMP build of adlaplace (has_openmp() == TRUE). Sets
# NOT_CRAN=true for full vignette mode and ADLAPLACE_VIGNETTE_NUM_THREADS for
# multi-thread ad_pack / fit calls (see vignettes/vignette_full.R).
#
# Writes .vignette-timing.txt and exits non-zero if a knit fails or OpenMP is off.
set -uo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root"

vignette_dir="${root}/adlaplace/vignettes"
results="${root}/.vignette-timing.txt"
out_dir="${root}/.vignette-out"
: > "${results}"
mkdir -p "${out_dir}"

THREADS_A="${ADLAPLACE_TIMING_THREADS_A:-1}"
THREADS_B="${ADLAPLACE_TIMING_THREADS_B:-4}"

echo "===== sessionInfo ====="
Rscript -e 'print(R.version.string); print(utils::packageVersion("adlaplace")); print(utils::packageVersion("RCppAD")); cat("has_openmp=", adlaplace::has_openmp(), "\n", sep=""); print(sessionInfo())' || true

echo "===== OpenMP gate ====="
set +e
Rscript -e 'if (!isTRUE(adlaplace::has_openmp())) quit(save="no", status=2L)'
omp_ec=$?
set -e
if [ "${omp_ec}" -ne 0 ]; then
  echo "ERROR: adlaplace was not built with OpenMP; cannot compare thread timings."
  exit 1
fi

echo "===== knit vignettes (${THREADS_A} then ${THREADS_B} threads) ====="
set +e
NOT_CRAN=true \
ADLAPLACE_TIMING_THREADS_A="${THREADS_A}" \
ADLAPLACE_TIMING_THREADS_B="${THREADS_B}" \
ADLAPLACE_VIGNETTE_DIR="${vignette_dir}" \
ADLAPLACE_VIGNETTE_OUT="${out_dir}" \
ADLAPLACE_VIGNETTE_RESULTS="${results}" \
Rscript "${root}/.github/scripts/time-adlaplace-vignettes.R"
ec=$?
set -e

echo "===== vignette timing summary ====="
if [ -f "${results}" ]; then
  echo "vignette|threads|elapsed_sec|exit_code"
  cat "${results}"
fi

exit "${ec}"
