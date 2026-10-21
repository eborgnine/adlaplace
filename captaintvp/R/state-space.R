#' State-space system for a CAPTAIN RW/IRW trend plus AR perturbation.
#'
#' Variances are normalized so that the concentrated observational or AR
#' innovation variance equals 1. Trend disturbance variances equal `nvr`.
#' The initial covariance follows CAPTAIN's default `1e5` on nonstationary
#' trend states. Stationary AR states use their stationary covariance.
#'
#' @param ARp AR polynomial excluding the leading 1, or `NULL`.
#' @param TVP `0` for a random walk and `1` for an integrated random walk.
#' @param nvr Trend noise-variance ratio, or `NULL` for an AR-only model.
#' @param P0 Initial variance of each nonstationary trend state.
#' @return A list of KFAS system matrices and state positions.
#' @noRd
trend_ar_system <- function(ARp, TVP, nvr, P0 = CAPTAIN_P0, diffuse = FALSE) {
  TVP <- validate_tvp(TVP)
  has_trend <- !is.null(nvr)
  if (has_trend) {
    nvr <- validate_nvr(nvr)
  }
  ARp <- validate_arp(ARp)
  has_ar <- length(ARp) > 0L
  if (!has_trend && !has_ar) {
    stop("specify a trend nvr, an AR polynomial, or both", call. = FALSE)
  }
  if (!is.finite(P0) || P0 <= 0) {
    stop("P0 must be a positive initial variance", call. = FALSE)
  }

  blocks <- list()
  if (has_trend) {
    blocks$trend <- trend_block(TVP, nvr, P0, diffuse = diffuse)
  }
  if (has_ar) {
    blocks$ar <- ar_block(ARp, P0)
  }

  m <- sum(vapply(blocks, function(block) nrow(block$T), integer(1L)))
  k <- sum(vapply(blocks, function(block) ncol(block$R), integer(1L)))
  Tmat <- matrix(0, m, m)
  Z <- matrix(0, 1L, m)
  R <- matrix(0, m, k)
  Q <- matrix(0, k, k)
  P1 <- matrix(0, m, m)
  P1inf <- matrix(0, m, m)
  names <- character(m)
  cursor <- 0L
  noise <- 0L
  positions <- list()
  for (nm in names(blocks)) {
    block <- blocks[[nm]]
    idx <- cursor + seq_len(nrow(block$T))
    qidx <- noise + seq_len(ncol(block$R))
    Tmat[idx, idx] <- block$T
    Z[1L, idx] <- block$Z
    R[idx, qidx] <- block$R
    Q[qidx, qidx] <- block$Q
    P1[idx, idx] <- block$P1
    P1inf[idx, idx] <- block$P1inf
    names[idx] <- block$names
    positions[[nm]] <- idx
    cursor <- cursor + length(idx)
    noise <- noise + length(qidx)
  }

  list(
    Z = Z,
    T = Tmat,
    R = R,
    Q = Q,
    P1 = P1,
    P1inf = P1inf,
    a1 = matrix(0, m, 1L),
    H = if (has_ar) 0 else 1,
    names = names,
    trend = positions$trend,
    ar = positions$ar,
    nvr = if (has_trend) nvr else NA_real_,
    ARp = ARp
  )
}

#' @noRd
trend_block <- function(TVP, nvr, P0, diffuse = FALSE) {
  n_state <- if (TVP == 0L) 1L else 2L
  P1 <- if (diffuse) diag(0, n_state) else diag(P0, n_state)
  P1inf <- if (diffuse) diag(1, n_state) else diag(0, n_state)
  if (TVP == 0L) {
    list(
      T = matrix(1),
      Z = matrix(1),
      R = matrix(1),
      Q = matrix(nvr),
      P1 = P1,
      P1inf = P1inf,
      names = "level"
    )
  } else {
    list(
      T = rbind(c(1, 1), c(0, 1)),
      Z = matrix(c(1, 0), nrow = 1L),
      R = matrix(c(0, 1), ncol = 1L),
      Q = matrix(nvr),
      P1 = P1,
      P1inf = P1inf,
      names = c("level", "slope")
    )
  }
}

#' Companion form of `z_t + a1 z_{t-1} + ... + ap z_{t-p} = e_t`.
#' @noRd
ar_block <- function(ARp, P0) {
  p <- length(ARp)
  Tmat <- matrix(0, p, p)
  Tmat[1L, ] <- -ARp
  if (p > 1L) {
    Tmat[2:p, 1:(p - 1L)] <- diag(p - 1L)
  }
  R <- matrix(0, p, 1L)
  R[1L, 1L] <- 1
  Q <- matrix(1)
  process_var <- tcrossprod(R %*% Q, R) # R Q R'
  poles <- abs(eigen(Tmat, symmetric = FALSE, only.values = TRUE)$values)
  P1 <- if (max(poles) < 1 - 1e-8) {
    lyapunov(Tmat, process_var)
  } else {
    diag(P0, p)
  }
  list(
    T = Tmat,
    Z = matrix(c(1, rep(0, p - 1L)), nrow = 1L),
    R = R,
    Q = Q,
    P1 = P1,
    P1inf = diag(0, p),
    names = paste0("ar", seq_len(p))
  )
}

#' Solve `P = T P T' + Q` for a stable transition.
#' @noRd
lyapunov <- function(Tmat, Q) {
  m <- nrow(Tmat)
  lhs <- diag(m * m) - kronecker(Tmat, Tmat)
  vec_p <- solve(lhs, as.numeric(Q))
  P <- matrix(vec_p, m, m)
  (P + t(P)) / 2
}

#' Filter or smooth a trend plus AR model and rescale variances.
#' @noRd
smooth_univ <- function(y, ARp, TVP, nvr, sm = 1, P0 = CAPTAIN_P0,
                         diffuse = FALSE) {
  parsed <- as_captain_series(y)
  values <- parsed$values
  system <- trend_ar_system(ARp, TVP, nvr, P0 = P0, diffuse = diffuse)
  model <- SSModel(
    values ~ SSMcustom(
      Z = system$Z,
      T = system$T,
      R = system$R,
      Q = system$Q,
      a1 = system$a1,
      P1 = system$P1,
      P1inf = system$P1inf,
      state_names = system$names
    ) - 1,
    H = system$H
  )
  sm <- validate_sm(sm)
  kfs <- KFS(
    model,
    simplify = FALSE,
    smoothing = if (sm == 1L) c("state", "mean") else "none"
  )
  alpha <- if (sm == 1L) kfs$alphahat else kfs$att
  V <- if (sm == 1L) kfs$V else kfs$Ptt
  alpha <- as.matrix(alpha)
  if (nrow(alpha) != length(values) && ncol(alpha) == length(values)) {
    alpha <- t(alpha)
  }
  V <- normalize_covariance(V, ncol(alpha), nrow(alpha))
  sigma2 <- concentrated_scale(kfs, is.finite(values))

  level <- if (is.null(system$trend)) {
    rep(0, nrow(alpha))
  } else {
    alpha[, system$trend[1L]]
  }
  comp <- if (is.null(system$ar)) {
    rep(0, nrow(alpha))
  } else {
    alpha[, system$ar[1L]]
  }
  fit <- as.numeric(alpha %*% t(system$Z))
  fit_var <- signal_variance(V, system$Z)
  trend_var <- if (is.null(system$trend)) {
    rep(0, nrow(alpha))
  } else {
    V[system$trend[1L], system$trend[1L], ]
  }

  fitse <- sqrt(pmax(fit_var, 0) * sigma2)
  trendse <- sqrt(pmax(as.numeric(trend_var), 0) * sigma2)
  state_se <- sqrt(pmax(apply(V, 3L, diag), 0) * sigma2)
  if (is.null(dim(state_se))) {
    state_se <- matrix(state_se, nrow = 1L)
  }
  state_se <- t(state_se)
  colnames(alpha) <- system$names
  colnames(state_se) <- system$names

  y0 <- values
  missing <- !is.finite(values)
  y0[missing] <- fit[missing]

  list(
    fit = restore_series(fit, parsed$tsp),
    fitse = restore_series(fitse, parsed$tsp),
    trend = restore_series(level, parsed$tsp),
    trendse = restore_series(trendse, parsed$tsp),
    comp = restore_series(comp, parsed$tsp),
    y0 = restore_series(y0, parsed$tsp),
    sigma2 = sigma2,
    state = alpha,
    state_se = state_se,
    system = system
  )
}

#' @noRd
signal_variance <- function(V, Z) {
  vapply(seq_len(dim(V)[3L]), function(t) {
    as.numeric(Z %*% V[, , t] %*% t(Z))
  }, numeric(1L))
}

#' @noRd
normalize_covariance <- function(V, m, n) {
  if (length(dim(V)) == 3L) {
    return(V)
  }
  if (m == 1L) {
    return(array(as.numeric(V), dim = c(1L, 1L, n)))
  }
  stop("unexpected state-covariance dimensions from KFAS", call. = FALSE)
}

#' Concentrated scale `mean(v^2 / F)` for normalized system matrices.
#' @noRd
concentrated_scale <- function(kfs, observed) {
  innov <- as.numeric(kfs$v)
  variance <- as.numeric(kfs$F)
  ok <- observed & is.finite(innov) & is.finite(variance) & variance > 0
  if (!any(ok)) {
    stop("the Kalman filter produced no finite innovations", call. = FALSE)
  }
  sum(innov[ok]^2 / variance[ok]) / sum(ok)
}

#' @noRd
validate_tvp <- function(TVP) {
  if (length(TVP) != 1L || is.na(TVP) || !(as.numeric(TVP) %in% c(0, 1))) {
    stop("TVP must be 0 (random walk) or 1 (integrated random walk)", call. = FALSE)
  }
  as.integer(TVP)
}

#' @noRd
validate_nvr <- function(nvr) {
  if (length(nvr) != 1L || !is.finite(nvr) || nvr < 0) {
    stop("nvr must be a single non-negative noise-variance ratio", call. = FALSE)
  }
  as.numeric(nvr)
}

#' @noRd
validate_sm <- function(sm) {
  if (length(sm) != 1L || is.na(sm) || !(as.numeric(sm) %in% c(0, 1))) {
    stop("sm must be 1 (smooth) or 0 (filter only)", call. = FALSE)
  }
  as.integer(sm)
}

#' @noRd
validate_arp <- function(ARp) {
  if (is.null(ARp) || length(ARp) == 0L) {
    return(numeric())
  }
  ARp <- as.numeric(ARp)
  if (any(!is.finite(ARp))) {
    stop("ARp must contain finite AR polynomial coefficients", call. = FALSE)
  }
  ARp
}

#' @noRd
reject_trend_ar <- function(ARt) {
  if (!(is.null(ARt) || length(ARt) == 0L)) {
    stop(
      "IAR and DIAR trend polynomials (ARt) are not implemented",
      call. = FALSE
    )
  }
}

#' @noRd
reject_tar <- function(tar) {
  if (is.null(tar) || length(tar) == 0L) {
    return(invisible())
  }
  if (length(tar) == 1L && tar == 0) {
    return(invisible())
  }
  stop("IAR and DIAR trend lags (tar) are not implemented", call. = FALSE)
}

#' @noRd
reject_intervention <- function(Int, n) {
  if (is.null(Int) || length(Int) == 0L) {
    return(invisible())
  }
  Int <- as.numeric(Int)
  if (length(Int) == 1L && Int == 0) {
    return(invisible())
  }
  if (length(Int) == n && all(Int == 0 | is.na(Int))) {
    return(invisible())
  }
  if (any(Int > 0, na.rm = TRUE)) {
    stop("variance intervention (Int) is not implemented", call. = FALSE)
  }
  stop("Int must be 0 or a vector of zeros with the same length as y", call. = FALSE)
}
