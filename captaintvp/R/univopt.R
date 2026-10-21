#' Estimate a trend noise-variance ratio and an AR perturbation
#'
#' @description
#' `univopt()` performs CAPTAIN's sequential spectral decomposition. It smooths
#' an initial trend at `nvr0`, fits the requested AR lags to the observed
#' trend residuals, and rescales the trend NVR by the ratio of the initial
#' concentrated disturbance variance to the AR innovation variance. This keeps
#' the AR model from absorbing the trend.
#'
#' @param y Numeric series. Missing values are retained for the trend smoother
#'   and omitted from the AR regression.
#' @param p Lags in the perturbation AR. `0` or `NULL` fits no AR component.
#'   A vector such as `c(1:3, 12, 13)` is a subset polynomial; excluded lags
#'   are fixed at zero.
#' @param TVP `0` for a random walk and `1` (default) for an integrated random
#'   walk.
#' @param nvr0 Initial trend noise-variance ratio. The default is `1`.
#' @param tar Lags of an AR polynomial for an IAR or DIAR trend. Only `0` is
#'   implemented.
#' @param out `1` (default) prints the coefficient table. `0` is silent.
#' @param Int Reserved variance-intervention indicator. Only `0` is implemented.
#'
#' @return A list in CAPTAIN output order: `nvr`, `ARp`, `ARt`, `ARpse`, and
#'   `ARtse`. `ARp` excludes the leading unity coefficient.
#' @export
#'
#' @examples
#' y <- as.numeric(AirPassengers)
#' est <- univopt(y, 1:2, TVP = 1, nvr0 = 1e-4, out = 0)
#' est$nvr
univopt <- function(y, p = 0, TVP = 1, nvr0 = 1, tar = 0, out = 1, Int = 0) {
  reject_tar(tar)
  parsed <- as_captain_series(y)
  reject_intervention(Int, length(parsed$values))
  nvr0 <- validate_nvr(nvr0)
  lags <- validate_lags(p)
  TVP <- validate_tvp(TVP)
  if (!isTRUE(out == 0) && !isTRUE(out == 1) && !isTRUE(out == 0L) && !isTRUE(out == 1L)) {
    stop("out must be 0 or 1", call. = FALSE)
  }

  initial <- smooth_univ(parsed$values, ARp = NULL, TVP = TVP, nvr = nvr0, sm = 1)
  pert <- parsed$values - as.numeric(initial$trend)
  pert[!is.finite(parsed$values)] <- NA_real_
  ar <- ar_subset(pert, lags)
  nvr <- if (!length(lags)) {
    nvr0
  } else {
    nvr0 * initial$sigma2 / ar$sigma2
  }

  result <- list(
    nvr = nvr,
    ARp = ar$ARp,
    ARt = numeric(),
    ARpse = ar$ARpse,
    ARtse = numeric(),
    sigma2 = initial$sigma2,
    sigma2_ar = ar$sigma2
  )
  if (out == 1 || out == 1L) {
    print_univopt(result, lags, TVP)
  }
  result
}

#' @noRd
print_univopt <- function(result, lags, TVP) {
  cat("ESTIMATION OF TREND+AR MODEL\n")
  cat("AR Model for Perturbations:\n")
  cat("==========================\n")
  if (!length(lags)) {
    cat("none\n")
  } else {
    cat(sprintf("%4s %10s %10s %10s\n", "AR", "ARp", "S.E.", "T"))
    for (lag in lags) {
      estimate <- result$ARp[lag]
      se <- result$ARpse[lag]
      cat(sprintf(
        "%4d %10.4f %10.4f %10.4f\n",
        lag, estimate, se, estimate / se
      ))
    }
  }
  cat(sprintf("Final trend NVR estimate: %.4e\n", result$nvr))
  cat(sprintf("Integration order of trend: %d\n", TVP + 1L))
}
