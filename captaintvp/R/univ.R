#' Trend plus autoregressive unobserved-components model
#'
#' @description
#' `univ()` estimates a random-walk or integrated-random-walk trend and,
#' optionally, a stochastic AR perturbation. Missing values are interpolated
#' from both directions. Missing values beyond either end of the observed
#' samples are backcasts or forecasts. Standard errors come from the Kalman
#' filter when `sm = 0` and from the fixed-interval smoother when `sm = 1`.
#'
#' The AR polynomial omits its leading 1. A coefficient of `-0.9` is a first-order
#' pole at `0.9`, because the model is `z_t + a1 z_{t-1} + ... = e_t`.
#' `TVP = 0` is a random walk and `TVP = 1` is an integrated random walk
#' (integration order 2).
#'
#' @param y Numeric series. `NA` marks samples to interpolate or forecast.
#' @param ARp AR polynomial for the perturbation, without the leading 1.
#'   `NULL` or `numeric()` estimates a trend plus white noise.
#' @param TVP `0` for RW and `1` (default) for IRW.
#' @param nvr Noise-variance ratio of the trend. `NULL` omits the trend and
#'   estimates the AR component alone.
#' @param ARt Reserved AR polynomial for an IAR or DIAR trend. Only the empty
#'   polynomial is implemented.
#' @param Int Reserved variance-intervention indicator. Only `0` is implemented.
#' @param sm `1` (default) smooths; `0` returns filtered estimates.
#'
#' @return A list with CAPTAIN's output names `fit`, `fitse`, `trend`,
#'   `trendse`, `comp`, and `y0`, followed by the concentrated scale `sigma2`,
#'   smoothed or filtered `state`, and `state_se`. `y0` is `y` with missing
#'   samples replaced by `fit`.
#' @export
#'
#' @examples
#' y <- fcast(as.numeric(AirPassengers), c(132, 144))
#' trend_only <- univ(y, numeric(), 1, 1e-4)
univ <- function(y, ARp = NULL, TVP = 1, nvr = NULL, ARt = NULL, Int = 0, sm = 1) {
  reject_trend_ar(ARt)
  parsed <- as_captain_series(y)
  reject_intervention(Int, length(parsed$values))
  smooth_univ(parsed$values, ARp = ARp, TVP = TVP, nvr = nvr, sm = sm)[c(
    "fit", "fitse", "trend", "trendse", "comp", "y0", "sigma2", "state", "state_se"
  )]
}
