#' Akaike information criterion for autoregressions
#'
#' @description
#' `aic()` selects an autoregressive order by minimizing
#' `n log(RSS / n) + 2 p` over conditional least-squares fits. The returned
#' polynomial includes the leading unity term, so `length(aic(y)) - 1` is the
#' selected order.
#'
#' @param y Numeric series. Missing values are omitted from every candidate
#'   fit. Models in one search share the sample defined by the largest order.
#' @param mod Maximum order, searched from order 1, or a length-2 vector
#'   `c(minimum, maximum)`. The default maximum is `min(floor(N / 2), 32)`.
#' @param out If `1`, print the criterion table and plot it. If `0`, the
#'   default, return silently.
#'
#' @return The selected polynomial `c(1, a1, ..., ap)`. The residual series is
#'   stored as attribute `eef`, and the full criterion table as attribute `aic`.
#' @export
#'
#' @examples
#' set.seed(1)
#' y <- arima.sim(n = 200, list(ar = c(0.6, -0.2)), sd = 1)
#' poly <- aic(y, 6)
#' length(poly) - 1
aic <- function(y, mod = NULL, out = 0) {
  if (!isTRUE(out == 0) && !isTRUE(out == 1) && !isTRUE(out == 0L) && !isTRUE(out == 1L)) {
    stop("out must be 0 or 1", call. = FALSE)
  }
  parsed <- as_captain_series(y)
  values <- parsed$values
  range <- aic_range(values, mod)
  table <- ar_aic_table(values, range[1L], range[2L])
  best <- which.min(table$aic)
  order <- table$order[best]
  coefs <- table$coefficients[[best]]
  arpoly <- c(1, coefs)
  residuals <- rep(NA_real_, length(values))
  residuals[table$index] <- table$residuals[[best]]

  if (out == 1 || out == 1L) {
    print_aic(table, order)
  }

  attr(arpoly, "eef") <- restore_series(residuals, parsed$tsp)
  attr(arpoly, "aic") <- data.frame(
    order = table$order,
    aic = table$aic,
    sigma2 = table$sigma2
  )
  arpoly
}

#' @noRd
aic_range <- function(y, mod) {
  n <- length(y)
  if (is.null(mod) || length(mod) == 0L) {
    c(1L, min(floor(n / 2), 32L))
  } else if (length(mod) == 1L) {
    c(1L, as.integer(mod))
  } else if (length(mod) == 2L) {
    as.integer(mod)
  } else {
    stop("mod must be a maximum order or c(minimum, maximum)", call. = FALSE)
  }
}

#' Conditional least-squares AIC on a sample common to every order.
#' @noRd
ar_aic_table <- function(y, min_order, max_order) {
  if (any(!is.finite(c(min_order, max_order))) || min_order < 1L || max_order < min_order) {
    stop("AR orders must satisfy 1 <= minimum <= maximum", call. = FALSE)
  }
  index <- complete_ar_index(y, seq_len(max_order))
  if (length(index) <= max_order) {
    stop("not enough complete observations for the requested AR orders", call. = FALSE)
  }
  response <- y[index]
  orders <- seq.int(min_order, max_order)
  coefficients <- vector("list", length(orders))
  residuals <- vector("list", length(orders))
  aic <- numeric(length(orders))
  sigma2 <- numeric(length(orders))
  n <- length(index)
  for (i in seq_along(orders)) {
    lags <- seq_len(orders[i])
    X <- ar_design(y, index, lags)
    fit <- stats::lm.fit(X, response)
    resid <- response - as.numeric(X %*% fit$coefficients)
    rss <- sum(resid^2)
    sigma2[i] <- rss / n
    aic[i] <- n * log(sigma2[i]) + 2 * orders[i]
    coefficients[[i]] <- -as.numeric(fit$coefficients)
    residuals[[i]] <- resid
  }
  list(
    order = orders,
    aic = aic,
    sigma2 = sigma2,
    coefficients = coefficients,
    residuals = residuals,
    index = index
  )
}

#' Least squares for a subset AR polynomial.
#'
#' The polynomial is returned without its leading 1 and is expanded to length
#' `max(lags)`, with structural zeros at excluded lags.
#' @noRd
ar_subset <- function(y, lags) {
  lags <- validate_lags(lags)
  if (!length(lags)) {
    return(empty_ar_fit())
  }
  index <- complete_ar_index(y, lags)
  if (length(index) <= length(lags)) {
    stop("not enough complete observations to estimate the AR polynomial", call. = FALSE)
  }
  response <- y[index]
  X <- ar_design(y, index, lags)
  fit <- stats::lm.fit(X, response)
  phi <- as.numeric(fit$coefficients)
  if (any(!is.finite(phi))) {
    stop("the AR least-squares fit is rank deficient", call. = FALSE)
  }
  resid <- response - as.numeric(X %*% phi)
  n <- length(response)
  df <- n - length(lags)
  sigma2 <- sum(resid^2) / n
  sigma2_unbiased <- sum(resid^2) / df
  p <- ncol(X)
  R <- fit$qr$qr[seq_len(p), seq_len(p), drop = FALSE]
  Rinv <- backsolve(R, diag(p))
  xtx_inv <- tcrossprod(Rinv)
  pivot <- fit$qr$pivot
  xtx_inv <- xtx_inv[order(pivot), order(pivot), drop = FALSE]
  se_phi <- sqrt(diag(xtx_inv) * sigma2_unbiased)

  order <- max(lags)
  ARp <- numeric(order)
  ARpse <- rep(NA_real_, order)
  ARp[lags] <- -phi
  ARpse[lags] <- se_phi
  list(
    ARp = ARp,
    ARpse = ARpse,
    lags = lags,
    sigma2 = sigma2,
    n = n,
    second_moment = mean(response^2),
    residual_moment = mean(resid^2)
  )
}

#' @noRd
empty_ar_fit <- function() {
  list(
    ARp = numeric(),
    ARpse = numeric(),
    lags = integer(),
    sigma2 = NA_real_,
    n = 0L,
    second_moment = NA_real_,
    residual_moment = NA_real_
  )
}

#' @noRd
validate_lags <- function(p) {
  if (is.null(p) || length(p) == 0L) {
    return(integer())
  }
  p <- as.numeric(p)
  if (length(p) == 1L && p == 0) {
    return(integer())
  }
  if (any(!is.finite(p)) || any(p < 1) || any(abs(p - round(p)) > 1e-8)) {
    stop("p must be 0 or positive integer AR lags", call. = FALSE)
  }
  sort(unique(as.integer(round(p))))
}

#' Rows whose response and requested lags are all observed.
#' @noRd
complete_ar_index <- function(y, lags) {
  n <- length(y)
  start <- max(lags) + 1L
  if (start > n) {
    return(integer())
  }
  index <- seq.int(start, n)
  ok <- is.finite(y[index])
  for (lag in lags) {
    ok <- ok & is.finite(y[index - lag])
  }
  index[ok]
}

#' @noRd
ar_design <- function(y, index, lags) {
  X <- vapply(lags, function(lag) y[index - lag], numeric(length(index)))
  if (length(lags) == 1L) {
    X <- matrix(X, ncol = 1L)
  }
  X
}

#' @noRd
print_aic <- function(table, selected) {
  cat("AKAIKE INFORMATION CRITERIUM\n")
  cat(sprintf("%8s %14s\n", "Order", "AIC"))
  for (i in seq_along(table$order)) {
    cat(sprintf("%8d %14.4f\n", table$order[i], table$aic[i]))
  }
  cat(sprintf("Selected order: %d\n", selected))
  graphics::plot(
    table$order,
    table$aic,
    type = "b",
    xlab = "AR order",
    ylab = "AIC",
    main = "AIC for autoregressive models"
  )
  graphics::abline(v = selected, lty = 2)
}
