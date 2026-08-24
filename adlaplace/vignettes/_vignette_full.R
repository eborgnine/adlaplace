# Shared vignette mode switch (sourced from setup chunks; not a vignette itself).
# Full examples: CI (NOT_CRAN=true) and interactive knits.
# Abbreviated: R CMD check on CRAN (_R_CHECK_PACKAGE_NAME_ set, NOT_CRAN unset).
vignette_full <- identical(Sys.getenv("NOT_CRAN"), "true") ||
  !nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))

if (!vignette_full) {
  message(
    "Abbreviated vignette for R CMD check; full HTML at ",
    "https://eborgnine.github.io/adlaplace/"
  )
}
