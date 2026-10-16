# Shared vignette mode switch (sourced from setup chunks; not a vignette itself).
# Name must not start with "_" — R CMD build excludes such files from the tarball.
# Full examples: CI (NOT_CRAN=true) and interactive knits.
# Abbreviated: R CMD check on CRAN (_R_CHECK_PACKAGE_NAME_ set, NOT_CRAN unset).
vignette_full <- identical(Sys.getenv("NOT_CRAN"), "true") ||
  !nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))

# Parallel team size for multi-thread ad_pack / fit calls in vignettes.
# Default: 8 for full knits, 2 for abbreviated R CMD check.
# Override with ADLAPLACE_VIGNETTE_NUM_THREADS (used by Windows timing CI).
vignette_num_threads <- {
  raw <- Sys.getenv("ADLAPLACE_VIGNETTE_NUM_THREADS", unset = "")
  if (!nzchar(raw)) {
    if (vignette_full) 8L else 2L
  } else {
    nt <- suppressWarnings(as.integer(raw))
    if (is.na(nt) || nt < 1L) {
      stop(
        "ADLAPLACE_VIGNETTE_NUM_THREADS must be a positive integer, got: ",
        raw,
        call. = FALSE
      )
    }
    nt
  }
}

if (!vignette_full) {
  message(
    "Abbreviated vignette for R CMD check; full HTML at ",
    "https://eborgnine.github.io/adlaplaceFem/"
  )
}
