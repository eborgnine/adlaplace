#' Prepare a series for forecasting, backcasting, or interpolation
#'
#' @description
#' `fcast()` replaces selected samples of a series with `NA`, matching the
#' CAPTAIN `fcast` convention. Modelling functions treat those missing samples
#' as interpolation, forecast, or backcast targets.
#'
#' @param y Numeric time series.
#' @param n Mask specification, defaulting to a 10-step forecast. A scalar
#'   `p` or row `c(0, p)` appends `p` missing values. A row `c(p, 0)` prepends
#'   `p` missing values. A row `c(p, q)` with `p > 0` and `q > 0` replaces
#'   samples `p` through `q` of the original series. A two-column matrix
#'   combines any of these operations. Sample numbers are 1-based and the
#'   endpoints are included.
#'
#' @return `y` with the requested samples replaced or extended by `NA`. A `ts`
#'   input keeps its frequency, and its start is moved backward when values
#'   are prepended.
#' @export
#'
#' @examples
#' y <- 1:30
#' x <- fcast(y, rbind(c(14, 19), c(22, 22), c(0, 10), c(5, 0)))
#' which(is.na(x))
fcast <- function(y, n = 10) {
  parsed <- as_captain_series(y)
  spec <- fcast_spec(n)
  values <- parsed$values
  original_n <- length(values)

  masked <- values
  extend_to <- original_n
  for (i in seq_len(nrow(spec$interpolate))) {
    p <- spec$interpolate[i, 1L]
    q <- spec$interpolate[i, 2L]
    if (q > extend_to) {
      masked <- c(masked, rep(NA_real_, q - extend_to))
      extend_to <- q
    }
    masked[p:q] <- NA_real_
  }

  if (spec$backcast > 0L) {
    masked <- c(rep(NA_real_, spec$backcast), masked)
  }
  if (spec$forecast > 0L) {
    masked <- c(masked, rep(NA_real_, spec$forecast))
  }
  restore_series(masked, parsed$tsp, -spec$backcast)
}

#' Decode the CAPTAIN `fcast` interval argument.
#' @noRd
fcast_spec <- function(n) {
  if (is.null(n) || length(n) == 0L) {
    stop("n must specify a forecast, backcast, or interpolation interval", call. = FALSE)
  }
  if (is.null(dim(n)) && length(n) == 1L) {
    n <- matrix(c(0, n), nrow = 1L)
  } else if (is.null(dim(n)) && length(n) == 2L) {
    n <- matrix(n, nrow = 1L)
  } else {
    n <- as.matrix(n)
  }
  if (ncol(n) != 2L) {
    stop("n must be a scalar, a length-2 vector, or a two-column matrix", call. = FALSE)
  }
  if (any(!is.finite(n)) || any(abs(n - round(n)) > 1e-8)) {
    stop("fcast sample numbers must be finite integers", call. = FALSE)
  }
  n <- matrix(as.integer(round(n)), ncol = 2L)
  if (any(n < 0L)) {
    stop("fcast sample numbers must be non-negative", call. = FALSE)
  }

  forecast <- 0L
  backcast <- 0L
  interpolate <- matrix(integer(), ncol = 2L)
  for (i in seq_len(nrow(n))) {
    p <- n[i, 1L]
    q <- n[i, 2L]
    if (p == 0L && q == 0L) {
      next
    } else if (p == 0L) {
      forecast <- forecast + q
    } else if (q == 0L) {
      backcast <- backcast + p
    } else if (p > q) {
      stop("interpolation interval p:q requires p <= q", call. = FALSE)
    } else {
      interpolate <- rbind(interpolate, c(p, q))
    }
  }
  list(forecast = forecast, backcast = backcast, interpolate = interpolate)
}

#' Split a vector or `ts` into values and an optional `tsp` attribute.
#' @noRd
as_captain_series <- function(y) {
  tsp_attr <- NULL
  if (stats::is.ts(y)) {
    tsp_attr <- stats::tsp(y)
    values <- as.numeric(y)
  } else if (is.matrix(y) || is.data.frame(y)) {
    if (ncol(y) != 1L) {
      stop("y must be a numeric vector or a single-column series", call. = FALSE)
    }
    values <- as.numeric(y[, 1L])
  } else if (is.numeric(y) && is.null(dim(y))) {
    values <- as.numeric(y)
  } else {
    stop("y must be a numeric vector or a univariate time series", call. = FALSE)
  }
  if (!length(values)) {
    stop("y is empty", call. = FALSE)
  }
  if (any(is.infinite(values))) {
    stop("y contains infinite values; use NA for missing samples", call. = FALSE)
  }
  list(values = values, tsp = tsp_attr)
}

#' Restore a series, optionally shifting its time origin by `start_shift` samples.
#' @noRd
restore_series <- function(values, tsp_attr, start_shift = 0) {
  values <- as.numeric(values)
  if (is.null(tsp_attr)) {
    return(values)
  }
  stats::ts(
    values,
    start = tsp_attr[1L] + start_shift / tsp_attr[3L],
    frequency = tsp_attr[3L]
  )
}
