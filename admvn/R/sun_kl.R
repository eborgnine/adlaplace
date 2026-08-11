#' Merge user options over defaults (named lists)
#' @keywords internal
.merge_opts <- function(defaults, opts = list()) {
  if (length(opts)) {
    utils::modifyList(defaults, as.list(opts))
  } else {
    defaults
  }
}

#' Merge user \code{optim_opts} over L-BFGS-B defaults for SUN KL fits
#' @keywords internal
.merge_optim_opts <- function(optim_opts = list(), maxit = 200, verbose = FALSE,
                              parscale = NULL) {
  out <- list(
    maxit = maxit,
    trace = if (isTRUE(verbose)) 1 else 0,
    REPORT = 10,
    factr = 1e3,
    pgtol = 1e-8
  )
  if (!is.null(parscale)) {
    out$parscale <- parscale
  }
  .merge_opts(out, optim_opts)
}

.evaluate_sun_slice <- function(par, quad,
                                family = c("22", "33", "33hs", "32hs", "43hs", "44"),
                                n_points = NULL,
                                n_shifts = NULL,
                                n_threads = 1L) {
  family <- match.arg(family)
  par <- as.numeric(par)
  w <- as.numeric(quad$w)
  if (!length(w) || any(!is.finite(w)) || sum(w) <= 0) {
    return(c(kl = Inf, neg_Elq = Inf))
  }
  w <- w / sum(w)
  x <- as.matrix(quad$x)

  if (family == "22") {
    if (length(par) != 10L) stop("par must have length 10 for SUN(2,2)")
    if (is.null(n_points)) n_points <- 64L
    if (is.null(n_shifts)) n_shifts <- 2L
    tape <- tryCatch(
      dsun22_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  } else if (family == "33") {
    if (length(par) != 21L) stop("par must have length 21 for SUN(3,3)")
    if (is.null(n_points)) n_points <- 64L
    if (is.null(n_shifts)) n_shifts <- 2L
    tape <- tryCatch(
      dsun_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  } else if (family == "33hs") {
    if (length(par) != 21L) {
      stop("par must have length 21 for SUN(3,3) hyperspherical")
    }
    if (is.null(n_points)) n_points <- 64L
    if (is.null(n_shifts)) n_shifts <- 2L
    tape <- tryCatch(
      dsun_hs_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  } else if (family == "32hs") {
    if (length(par) != 16L) {
      stop("par must have length 16 for SUN(3,2) hyperspherical")
    }
    if (is.null(n_points)) n_points <- 64L
    if (is.null(n_shifts)) n_shifts <- 2L
    tape <- tryCatch(
      dsun32_hs_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  } else if (family == "43hs") {
    if (length(par) != 29L) {
      stop("par must have length 29 for SUN(4,3) hyperspherical")
    }
    if (is.null(n_points)) n_points <- 64L
    if (is.null(n_shifts)) n_shifts <- 2L
    tape <- tryCatch(
      dsun43_hs_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  } else {
    if (length(par) != 36L) stop("par must have length 36 for SUN(4,4)")
    if (is.null(n_points)) n_points <- 256L
    if (is.null(n_shifts)) n_shifts <- 4L
    tape <- tryCatch(
      dsun44_fun(
        x = x,
        par_seed = par,
        weights = w,
        n_points = as.integer(n_points),
        n_shifts = as.integer(n_shifts),
        n_threads = as.integer(n_threads)
      ),
      error = function(e) NULL
    )
  }
  if (is.null(tape)) {
    return(c(kl = Inf, neg_Elq = Inf))
  }
  ## Weighted sum ∑ w log q via admvn (never sn::dsun)
  v <- tryCatch(
    tape$eval(par, log = TRUE, deriv = 0L)$value,
    error = function(e) NA_real_
  )
  if (!is.finite(v)) {
    return(c(kl = Inf, neg_Elq = Inf))
  }
  neg_Elq <- -as.numeric(v)
  if (is.null(quad$logp)) {
    return(c(kl = NA_real_, neg_Elq = neg_Elq))
  }
  c(
    kl = sum(w * as.numeric(quad$logp)) - as.numeric(v),
    neg_Elq = neg_Elq
  )
}

#' Slice KL / cross-entropy of a SUN(2,2) vs a weighted target
#'
#' Uses [dsun22_fun()] (not \pkg{sn}) so the score matches the AD optim objective.
#'
#' @param par Length-10 SUN(2,2) parameter vector.
#' @param quad List with \code{x}, \code{w}, and optionally \code{logp}
#'   (target log-density at each row of \code{x}).
#' @param n_points,n_shifts Orthant QMC settings (defaults \code{64}, \code{2}).
#' @param n_threads OpenMP threads for the density tape.
#' @return Named vector \code{c(kl, neg_Elq)}. \code{kl} is \code{NA} if
#'   \code{logp} is missing.
#' @export
evaluate_sun22_slice <- function(par, quad,
                                 n_points = 64L,
                                 n_shifts = 2L,
                                 n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "22",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

#' Slice KL / cross-entropy of a SUN(3,3) vs a weighted target
#'
#' Uses [dsun_fun()] (not \pkg{sn}) so the score matches the AD optim objective.
#'
#' @param par Length-21 SUN(3,3) parameter vector.
#' @param quad List with \code{x}, \code{w}, and optionally \code{logp}
#'   (target log-density at each row of \code{x}).
#' @param n_points,n_shifts Orthant QMC settings (defaults \code{64}, \code{2}).
#' @param n_threads OpenMP threads for the density tape.
#' @return Named vector \code{c(kl, neg_Elq)}. \code{kl} is \code{NA} if
#'   \code{logp} is missing.
#' @export
evaluate_sun33_slice <- function(par, quad,
                                 n_points = 64L,
                                 n_shifts = 2L,
                                 n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "33",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

#' Slice KL / cross-entropy of a SUN(4,4) vs a weighted target
#'
#' Uses [dsun44_fun()] (not \pkg{sn}).
#'
#' @inheritParams evaluate_sun33_slice
#' @param par Length-36 SUN(4,4) parameter vector.
#' @param n_points,n_shifts Defaults \code{256}, \code{4}.
#' @export
evaluate_sun44_slice <- function(par, quad,
                                 n_points = 256L,
                                 n_shifts = 4L,
                                 n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "44",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

.fit_sun_quad <- function(quad,
                          start,
                          lower,
                          upper,
                          obj_scale,
                          mc.cores,
                          verbose,
                          mvquad_opts,
                          optim_opts,
                          family = c("22", "33", "33hs", "32hs", "43hs", "44")) {
  family <- match.arg(family)
  start <- as.numeric(start)
  npar <- length(start)
  if (family == "22") {
    if (npar != 10L) stop("start must have length 10 for SUN(2,2)")
    make_dp <- make_sun22_params
    mq_default <- list(n_points = 64L, n_shifts = 2L)
    dsun_fun_create <- dsun22_fun
  } else if (family == "33") {
    if (npar != 21L) stop("start must have length 21 for SUN(3,3)")
    make_dp <- make_sun_params
    mq_default <- list(n_points = 64L, n_shifts = 2L)
    dsun_fun_create <- dsun_fun
  } else if (family == "33hs") {
    if (npar != 21L) {
      stop("start must have length 21 for SUN(3,3) hyperspherical")
    }
    make_dp <- make_sun_hs_params
    mq_default <- list(n_points = 64L, n_shifts = 2L)
    dsun_fun_create <- dsun_hs_fun
  } else if (family == "32hs") {
    if (npar != 16L) {
      stop("start must have length 16 for SUN(3,2) hyperspherical")
    }
    make_dp <- make_sun32_hs_params
    mq_default <- list(n_points = 64L, n_shifts = 2L)
    dsun_fun_create <- dsun32_hs_fun
  } else if (family == "43hs") {
    if (npar != 29L) {
      stop("start must have length 29 for SUN(4,3) hyperspherical")
    }
    make_dp <- make_sun43_hs_params
    mq_default <- list(n_points = 64L, n_shifts = 2L)
    dsun_fun_create <- dsun43_hs_fun
  } else {
    if (npar != 36L) stop("start must have length 36 for SUN(4,4)")
    make_dp <- make_sun44_params
    mq_default <- list(n_points = 256L, n_shifts = 4L)
    dsun_fun_create <- dsun44_fun
  }

  if (is.null(lower) || is.null(upper)) {
    b <- switch(family,
      "22" = sun22_bounds(),
      "33" = sun33_bounds(),
      "33hs" = sun33_hs_bounds(),
      "32hs" = sun32_hs_bounds(),
      "43hs" = sun43_hs_bounds(),
      sun44_bounds()
    )
    if (is.null(lower)) lower <- b$lower
    if (is.null(upper)) upper <- b$upper
  }

  mq <- .merge_opts(mq_default, mvquad_opts)
  mc.cores <- max(1L, as.integer(mc.cores))
  n_workers <- mc.cores
  gr_last <- NULL

  ## Normalize weights so optim objective matches evaluate_sun*_slice / KL.
  w <- as.numeric(quad$w)
  if (!length(w) || any(!is.finite(w)) || sum(w) <= 0) {
    stop("quad$w must be finite and sum to a positive value")
  }
  w <- w / sum(w)

  tape <- dsun_fun_create(
    x = as.matrix(quad$x),
    par_seed = start,
    weights = w,
    n_points = as.integer(mq$n_points),
    n_shifts = as.integer(mq$n_shifts),
    n_threads = mc.cores
  )
  ## Score via the same admvn tape used for optim (never sn::dsun).
  evaluate_tape <- function(par) {
    v <- tryCatch(
      tape$eval(as.numeric(par), log = TRUE, deriv = 0L)$value,
      error = function(e) NA_real_
    )
    if (!is.finite(v)) {
      return(c(kl = Inf, neg_Elq = Inf))
    }
    neg_Elq <- -as.numeric(v)
    if (is.null(quad$logp)) {
      return(c(kl = NA_real_, neg_Elq = neg_Elq))
    }
    c(
      kl = sum(w * as.numeric(quad$logp)) - as.numeric(v),
      neg_Elq = neg_Elq
    )
  }
  value_start <- evaluate_tape(start)

  ofs <- tape$optim_fns(n_threads = 0L)
  fn <- function(par) {
    v <- ofs$fn(par)
    if (!is.finite(v)) return(1e20)
    obj_scale * (-v)
  }
  gr <- function(par) {
    g <- obj_scale * (-ofs$gr(par))
    gr_last <<- g
    g
  }
  opt_control <- .merge_optim_opts(
    optim_opts, verbose = verbose,
    parscale = pmax(abs(start), 0.1)
  )
  fit <- stats::optim(
    par = start,
    fn = fn,
    gr = gr,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = opt_control
  )

  gradient <- tryCatch(gr(fit$par), error = function(e) {
    if (!is.null(gr_last)) gr_last else rep(NA_real_, npar)
  })
  if (!is.null(names(start))) {
    names(gradient) <- names(start)
  }

  value_fit <- evaluate_tape(fit$par)

  list(
    fit = fit,
    par = fit$par,
    dp = make_dp(fit$par),
    start = start,
    gradient = gradient,
    quad = quad,
    value = unname(value_fit[["neg_Elq"]]),
    value_start = unname(value_start[["neg_Elq"]]),
    kl = unname(value_fit[["kl"]]),
    kl_start = unname(value_start[["kl"]]),
    obj_scale = obj_scale,
    mc.cores = n_workers,
    engine = "admvn"
  )
}

#' Fit SUN(2,2) by minimising KL on a weighted quadrature
#'
#' Minimises the weighted cross-entropy \eqn{-\sum_i w_i \log q(x_i)} with
#' analytic gradients from [dsun22_fun()].
#'
#' @param quad List with \code{x} (\eqn{n\times 2}), \code{w}, and optionally
#'   \code{logp} (for returned KL).
#' @param start Length-10 start (see [make_sun22_start_from_normal()]).
#' @param lower,upper Box constraints (defaults from [sun22_bounds()]).
#' @param obj_scale Scale applied to the optim objective (default \code{1e4}).
#' @param mc.cores OpenMP threads for the density tape.
#' @param verbose If \code{TRUE}, set \code{optim} \code{trace = 1}.
#' @param mvquad_opts Named list; used keys \code{n_points}, \code{n_shifts}
#'   for orthant QMC (defaults \code{64}, \code{2}).
#' @param optim_opts Named list merged into \code{stats::optim} L-BFGS-B
#'   control.
#'
#' @return A list with \code{fit}, \code{par}, \code{dp}, \code{start},
#'   \code{gradient} (objective-scale at \code{par}), KL / cross-entropy
#'   summaries, and \code{engine}.
#' @export
fit_sun22_quad <- function(quad,
                           start,
                           lower = NULL,
                           upper = NULL,
                           obj_scale = 1e4,
                           mc.cores = 1L,
                           verbose = FALSE,
                           mvquad_opts = list(),
                           optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "22"
  )
}

#' Fit SUN(3,3) by minimising KL on a weighted quadrature
#'
#' Minimises the weighted cross-entropy \eqn{-\sum_i w_i \log q(x_i)} with
#' analytic gradients from [dsun_fun()].
#'
#' @param quad List with \code{x} (\eqn{n\times 3}), \code{w}, and optionally
#'   \code{logp} (for returned KL).
#' @param start Length-21 start (see [make_sun33_start_from_normal()]).
#' @param lower,upper Box constraints (defaults from [sun33_bounds()]).
#' @param obj_scale Scale applied to the optim objective (default \code{1e4}).
#' @param mc.cores OpenMP threads for the density tape.
#' @param verbose If \code{TRUE}, set \code{optim} \code{trace = 1}.
#' @param mvquad_opts Named list; used keys \code{n_points}, \code{n_shifts}
#'   for orthant QMC (defaults \code{64}, \code{2}).
#' @param optim_opts Named list merged into \code{stats::optim} L-BFGS-B
#'   control.
#'
#' @return A list with \code{fit}, \code{par}, \code{dp}, \code{start},
#'   \code{gradient} (objective-scale at \code{par}), KL / cross-entropy
#'   summaries, and \code{engine}.
#' @export
fit_sun33_quad <- function(quad,
                           start,
                           lower = NULL,
                           upper = NULL,
                           obj_scale = 1e4,
                           mc.cores = 1L,
                           verbose = FALSE,
                           mvquad_opts = list(),
                           optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "33"
  )
}

#' Slice KL / cross-entropy of a SUN(3,3) hyperspherical vs a weighted target
#'
#' Like [evaluate_sun33_slice()], but uses [dsun_hs_fun()] and
#' [make_sun_hs_params()].
#'
#' @inheritParams evaluate_sun33_slice
#' @param par Length-21 hyperspherical parameter vector
#'   (see [make_sun_hs_params()]).
#' @export
evaluate_sun33_hs_slice <- function(par, quad,
                                    n_points = 64L,
                                    n_shifts = 2L,
                                    n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "33hs",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

#' Slice KL / cross-entropy of a SUN(3,2) hyperspherical fit
#' @inheritParams evaluate_sun33_hs_slice
#' @export
evaluate_sun32_hs_slice <- function(par, quad,
                                    n_points = 64L,
                                    n_shifts = 2L,
                                    n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "32hs",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

#' Slice KL / cross-entropy of a SUN(4,3) hyperspherical fit
#' @inheritParams evaluate_sun33_hs_slice
#' @export
evaluate_sun43_hs_slice <- function(par, quad,
                                    n_points = 64L,
                                    n_shifts = 2L,
                                    n_threads = 1L) {
  .evaluate_sun_slice(
    par, quad,
    family = "43hs",
    n_points = n_points,
    n_shifts = n_shifts,
    n_threads = n_threads
  )
}

#' Fit SUN(3,3) hyperspherical by minimising KL on a weighted quadrature
#'
#' Like [fit_sun33_quad()], but uses [dsun_hs_fun()] and defaults from
#' [sun33_hs_bounds()] / [make_sun33_hs_start_from_normal()].
#'
#' @inheritParams fit_sun33_quad
#' @param start Length-21 hyperspherical start
#'   (see [make_sun33_hs_start_from_normal()]).
#' @param lower,upper Box constraints (defaults from [sun33_hs_bounds()]).
#' @export
fit_sun33_hs_quad <- function(quad,
                              start,
                              lower = NULL,
                              upper = NULL,
                              obj_scale = 1e4,
                              mc.cores = 1L,
                              verbose = FALSE,
                              mvquad_opts = list(),
                              optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "33hs"
  )
}

#' Fit SUN(3,2) hyperspherical by minimising KL on a weighted quadrature
#'
#' @inheritParams fit_sun33_hs_quad
#' @param start Length-16 hyperspherical start
#'   (see [make_sun32_hs_start_from_normal()]).
#' @param lower,upper Box constraints (defaults from [sun32_hs_bounds()]).
#' @export
fit_sun32_hs_quad <- function(quad,
                              start,
                              lower = NULL,
                              upper = NULL,
                              obj_scale = 1e4,
                              mc.cores = 1L,
                              verbose = FALSE,
                              mvquad_opts = list(),
                              optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "32hs"
  )
}

#' Fit SUN(4,3) hyperspherical by minimising KL on a weighted quadrature
#'
#' @inheritParams fit_sun33_hs_quad
#' @param start Length-29 hyperspherical start
#'   (see [make_sun43_hs_start_from_normal()]).
#' @param lower,upper Box constraints (defaults from [sun43_hs_bounds()]).
#' @export
fit_sun43_hs_quad <- function(quad,
                              start,
                              lower = NULL,
                              upper = NULL,
                              obj_scale = 1e4,
                              mc.cores = 1L,
                              verbose = FALSE,
                              mvquad_opts = list(),
                              optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "43hs"
  )
}

#' Fit SUN(4,4) by minimising KL on a weighted quadrature
#'
#' Like [fit_sun33_quad()] for SUN(4,4) via [dsun44_fun()]. Default QMC
#' budget is \code{n_points = 256}, \code{n_shifts = 4}.
#'
#' @inheritParams fit_sun33_quad
#' @param start Length-36 start (see [make_sun44_start_from_normal()]).
#' @export
fit_sun44_quad <- function(quad,
                           start,
                           lower = NULL,
                           upper = NULL,
                           obj_scale = 1e4,
                           mc.cores = 1L,
                           verbose = FALSE,
                           mvquad_opts = list(),
                           optim_opts = list()) {
  .fit_sun_quad(
    quad = quad,
    start = start,
    lower = lower,
    upper = upper,
    obj_scale = obj_scale,
    mc.cores = mc.cores,
    verbose = verbose,
    mvquad_opts = mvquad_opts,
    optim_opts = optim_opts,
    family = "44"
  )
}

#' Print SUN KL fit diagnostics (including gradient)
#'
#' @param fit Result from [fit_sun22_quad()] / [fit_sun33_quad()] /
#'   [fit_sun33_hs_quad()] / [fit_sun44_quad()]
#'   (or a toenail wrapper that preserves the same fields).
#' @param label Optional header string.
#' @return \code{fit}, invisibly.
#' @export
print_sun_kl_fit <- function(fit, label = NULL) {
  if (!is.null(label) && nzchar(label)) {
    cat(label, "\n", sep = "")
  }
  opt <- fit$fit
  conv <- if (!is.null(opt$convergence)) opt$convergence else NA_integer_
  counts <- if (!is.null(opt$counts)) {
    paste(names(opt$counts), as.integer(opt$counts), sep = "=", collapse = ", ")
  } else {
    NA_character_
  }
  g <- as.numeric(fit$gradient)
  g_inf <- if (length(g) && all(is.finite(g))) max(abs(g)) else NA_real_
  cat(
    sprintf(
      "  engine=%s  convergence=%s  counts={%s}\n",
      if (!is.null(fit$engine)) fit$engine else "?",
      as.character(conv),
      counts
    )
  )
  if (!is.null(opt$message) && nzchar(opt$message)) {
    cat("  message:", opt$message, "\n")
  }
  cat(sprintf(
    "  neg_Elq: start=%.6g  fit=%.6g\n",
    as.numeric(fit$value_start), as.numeric(fit$value)
  ))
  if (!is.null(fit$kl_start) && is.finite(fit$kl_start)) {
    cat(sprintf(
      "  kl:      start=%.6g  fit=%.6g\n",
      as.numeric(fit$kl_start), as.numeric(fit$kl)
    ))
  }
  cat(sprintf("  ||grad||_inf=%.6g  (objective scale", g_inf))
  if (!is.null(fit$obj_scale)) {
    cat(sprintf(", obj_scale=%g", as.numeric(fit$obj_scale)))
  }
  cat(")\n")
  if (length(g)) {
    cat("  gradient:\n")
    print(g)
  }
  invisible(fit)
}
