#' Near-monotone mGP FEM helpers
#'
#' Ported from Ziang Zhang's BayesGP near-monotone FEM
#' (`Compute_Prec`, `Compute_Prec_rev`, `Compute_Design`) with Gauss-Legendre
#' quadrature in place of a left-point Riemann sum, and
#' [splines::splineDesign()] in place of **fda**.
#'
#' @name mgp-fem
#' @keywords internal
#' @noRd
NULL

#' Default original-scale shift so that x + c > 0 on the domain.
#' @noRd
mgp_shift <- function(c, domain) {
  domain <- as.numeric(domain)
  if (!length(domain) || any(!is.finite(domain))) {
    stop("mgp domain must be finite numeric", call. = FALSE)
  }
  if (!is.null(c)) {
    if (length(c) != 1L || !is.finite(c)) {
      stop("c must be a single finite value", call. = FALSE)
    }
    return(as.numeric(c))
  }
  lower <- min(domain)
  if (lower > 0) {
    return(0)
  }
  scale <- max(1, diff(range(domain)), max(abs(domain)))
  -lower + 1e-6 * scale
}

#' Resolve reference location on a domain.
#' @noRd
mgp_reference <- function(x, ref_value) {
  domain <- range(as.numeric(x), na.rm = TRUE)
  if (is.character(ref_value)) {
    if (length(ref_value) != 1L) {
      stop("ref_value must be a scalar", call. = FALSE)
    }
    key <- tolower(ref_value)
    if (identical(key, "left")) {
      return(domain[[1L]])
    }
    if (identical(key, "middle")) {
      return(mean(domain))
    }
    if (identical(key, "right")) {
      return(domain[[2L]])
    }
    stop("ref_value must be numeric or one of 'left', 'middle', 'right'", call. = FALSE)
  }
  if (!is.numeric(ref_value) || length(ref_value) != 1L || !is.finite(ref_value)) {
    stop("ref_value must be numeric or one of 'left', 'middle', 'right'", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  if (ref_value < domain[[1L]] - tol || ref_value > domain[[2L]] + tol) {
    stop("ref_value must lie inside the mgp region", call. = FALSE)
  }
  as.numeric(ref_value)
}

#' Require x + c > 0 over the region (original scale).
#' @noRd
mgp_check_domain <- function(region, c) {
  region <- range(as.numeric(region), na.rm = TRUE)
  if (any(region + c <= 0)) {
    stop(
      "mgp requires x + c > 0 over the region; got min(x) + c = ",
      region[[1L]] + c,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Normalized boundary (Box-Cox) transform of centred coordinate u.
#'
#' With \eqn{\lambda = (a-1)/a} and internal shift `cs`,
#' \eqn{m(u) = ((u+cs)^\lambda - cs^\lambda)/(\lambda cs^{\lambda-1})},
#' or \eqn{cs\log((u+cs)/cs)} when \eqn{\lambda = 0}. For `a = Inf`,
#' returns `u`.
#' @noRd
mgp_boxcox <- function(u, a, cs) {
  u <- as.numeric(u)
  if (isTRUE(all.equal(a, Inf))) {
    return(u)
  }
  if (!is.finite(a) || a == 0) {
    stop("a must be finite and nonzero, or Inf", call. = FALSE)
  }
  if (!is.finite(cs) || cs <= 0) {
    stop("internal shift cs must be positive", call. = FALSE)
  }
  arg <- u + cs
  if (any(arg <= 0)) {
    stop("all evaluation points must satisfy u + cs > 0", call. = FALSE)
  }
  lambda <- (a - 1) / a
  if (isTRUE(all.equal(lambda, 0))) {
    return(cs * log(arg / cs))
  }
  ((arg^lambda) - (cs^lambda)) / (lambda * cs^(lambda - 1))
}

#' Open cubic knot vector on \code{[0, region_max]} with fda-style k bases, drop 1:2.
#'
#' fda `create.bspline.basis(nbasis = k, norder = 4)` uses `k - 2` breakpoints
#' (including endpoints), i.e. `k - 4` interior knots. Dropping the first two
#' bases enforces `f(0) = f'(0) = 0`.
#' @noRd
mgp_side_basis <- function(region_max, k) {
  k <- as.integer(k)
  if (!is.finite(k) || k < 4L) {
    stop("k must be an integer >= 4", call. = FALSE)
  }
  region_max <- as.numeric(region_max)[1L]
  if (!is.finite(region_max) || region_max <= 0) {
    stop("side region length must be positive", call. = FALSE)
  }
  n_interior <- k - 4L
  knots <- open_knot_vector(c(0, region_max), degree = 3L, n_interior = n_interior)
  list(knots = knots, degree = 3L, drop = c(1L, 2L), n_kept = k - 2L)
}

#' Alpha coefficient on a side: positive uses -1/(a(u+cs)), reverse 1/(a(cs-v)).
#' @noRd
mgp_alpha <- function(u, a, cs, side = c("positive", "reverse")) {
  side <- match.arg(side)
  if (isTRUE(all.equal(a, Inf))) {
    return(rep(0, length(u)))
  }
  if (identical(side, "positive")) {
    return(-1 / (a * (u + cs)))
  }
  1 / (a * (cs - u))
}

#' FEM precision G - M + C on one side via Gauss-Legendre quadrature.
#' @noRd
mgp_side_precision <- function(knots, a, cs, side = c("positive", "reverse"),
                               n_quad = 10L, drop = c(1L, 2L)) {
  side <- match.arg(side)
  degree <- 3L
  n_quad <- as.integer(n_quad)
  if (!is.finite(n_quad) || n_quad < 2L) {
    stop("n_quad must be an integer >= 2", call. = FALSE)
  }
  gl <- gauss_legendre(n_quad)
  spans <- knot_spans(knots)
  nb <- n_basis_knots(knots, degree)
  G <- Matrix::Matrix(0, nb, nb, sparse = TRUE)
  M <- Matrix::Matrix(0, nb, nb, sparse = TRUE)
  C <- Matrix::Matrix(0, nb, nb, sparse = TRUE)

  for (s in seq_along(spans$left)) {
    ab <- gl_on_interval(spans$left[s], spans$right[s], gl)
    B1 <- bspline_eval(knots, ab$x, degree, 1L)
    B2 <- bspline_eval(knots, ab$x, degree, 2L)
    alpha <- mgp_alpha(ab$x, a = a, cs = cs, side = side)
    w <- ab$w
    # G_ij = ∫ (α B_i') (α B_j') dx
    B1a <- Matrix::Diagonal(x = alpha) %*% B1
    G <- G + Matrix::crossprod(Matrix::Diagonal(x = w) %*% B1a, B1a)
    # C_ij = ∫ B_i'' B_j'' dx
    C <- C + Matrix::crossprod(Matrix::Diagonal(x = w) %*% B2, B2)
    # M = ss(∫ (α B_i') B_j'' dx)
    cross <- Matrix::crossprod(Matrix::Diagonal(x = w) %*% B1a, B2)
    M <- M + cross + Matrix::t(cross)
  }

  Q <- G - M + C
  if (length(drop)) {
    keep <- setdiff(seq_len(nb), as.integer(drop))
    Q <- Q[keep, keep, drop = FALSE]
  }
  as_dgc_matrix(Matrix::drop0(Matrix::forceSymmetric(Q)))
}

#' Design matrix on one side (centred coordinates), dropping boundary columns.
#' @noRd
mgp_side_design <- function(u, knots, drop = c(1L, 2L)) {
  B <- bspline_eval(knots, u, degree = 3L, derivs = 0L)
  B <- methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix")
  if (length(drop)) {
    keep <- setdiff(seq_len(ncol(B)), as.integer(drop))
    B <- B[, keep, drop = FALSE]
  }
  B
}

#' Assemble one- or two-sided mGP FEM (knots, Q, log_det, rank).
#' @noRd
mgp_fem <- function(domain, x0, a, cs, k, n_quad = 10L) {
  domain <- range(as.numeric(domain), na.rm = TRUE)
  x0 <- as.numeric(x0)[1L]
  pos_len <- max(0, domain[[2L]] - x0)
  neg_len <- max(0, x0 - domain[[1L]])
  two_sided <- neg_len > 0 && pos_len > 0

  blocks_Q <- list()
  side_info <- list()

  if (pos_len > 0) {
    basis <- mgp_side_basis(pos_len, k)
    Q_pos <- mgp_side_precision(
      basis$knots, a = a, cs = cs, side = "positive",
      n_quad = n_quad, drop = basis$drop
    )
    blocks_Q[[length(blocks_Q) + 1L]] <- Q_pos
    side_info$positive <- basis
  }
  if (neg_len > 0) {
    basis <- mgp_side_basis(neg_len, k)
    Q_neg <- mgp_side_precision(
      basis$knots, a = a, cs = cs, side = "reverse",
      n_quad = n_quad, drop = basis$drop
    )
    blocks_Q[[length(blocks_Q) + 1L]] <- Q_neg
    side_info$negative <- basis
  }

  if (!length(blocks_Q)) {
    stop("mgp region must have positive length on at least one side of x0", call. = FALSE)
  }
  Q <- if (length(blocks_Q) == 1L) {
    blocks_Q[[1L]]
  } else {
    as_dgc_matrix(Matrix::bdiag(blocks_Q))
  }
  Q <- Matrix::forceSymmetric(Q, uplo = "U")
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  if (!is.finite(log_det)) {
    stop("mgp precision Q is not positive definite (log_det not finite)", call. = FALSE)
  }
  list(
    Q = Q,
    log_det = log_det,
    rank = ncol(Q),
    sides = side_info,
    pos_len = pos_len,
    neg_len = neg_len,
    two_sided = two_sided,
    x0 = x0,
    a = a,
    cs = cs,
    k = as.integer(k),
    n_quad = as.integer(n_quad)
  )
}

#' Design at original-scale x using a cached mgp_fem() result.
#' @noRd
mgp_fem_design <- function(x, fem) {
  x <- as.numeric(x)
  u <- x - fem$x0
  n <- length(u)
  blocks <- list()

  fill_side <- function(u_side, basis) {
    B <- matrix(0, n, basis$n_kept)
    active <- which(u_side > 0)
    if (length(active)) {
      B[active, ] <- as.matrix(
        mgp_side_design(u_side[active], basis$knots, drop = basis$drop)
      )
    }
    as_dgc_matrix(Matrix::Matrix(B, sparse = TRUE))
  }

  if (!is.null(fem$sides$positive)) {
    blocks[[length(blocks) + 1L]] <- fill_side(u, fem$sides$positive)
  }
  if (!is.null(fem$sides$negative)) {
    blocks[[length(blocks) + 1L]] <- fill_side(-u, fem$sides$negative)
  }

  if (length(blocks) == 1L) {
    return(blocks[[1L]])
  }
  as_dgc_matrix(cbind(blocks[[1L]], blocks[[2L]]))
}

#' Unit predictive SD via Green's-function quadrature (centred coordinates).
#'
#' Matches BayesGP `sqrt(SS_cov(u + h, u, c = cs, a = a))` for supported `a`.
#' @noRd
mgp_psd_unit <- function(u, h, a, cs) {
  u <- as.numeric(u)
  h <- as.numeric(h)
  n <- max(length(u), length(h))
  u <- rep_len(u, n)
  h <- rep_len(h, n)
  vapply(seq_len(n), function(i) {
    ui <- u[[i]]
    hi <- h[[i]]
    if (!is.finite(ui) || !is.finite(hi)) {
      stop("u and h must be finite", call. = FALSE)
    }
    if (hi == 0) {
      return(0)
    }
    if (hi < 0) {
      # SS_cov assumes s > t; swap endpoints.
      return(mgp_psd_unit(ui + hi, -hi, a = a, cs = cs))
    }
    lo <- ui
    hi_end <- ui + hi
    if (lo + cs <= 0 || hi_end + cs <= 0) {
      stop("mgp PSD requires u + cs > 0 and u + h + cs > 0", call. = FALSE)
    }
    var_inc <- mgp_green_var(lo, hi_end, a = a, cs = cs)
    if (var_inc < 0 && var_inc > -1e-12) {
      var_inc <- 0
    }
    if (var_inc < 0) {
      stop("mgp PSD integral returned a negative variance", call. = FALSE)
    }
    sqrt(var_inc)
  }, numeric(1))
}

#' ∫_lo^hi ((g(hi) - g(τ)) / w(τ))^2 dτ
#' @noRd
mgp_green_var <- function(lo, hi, a, cs) {
  if (isTRUE(all.equal(a, Inf))) {
    h <- hi - lo
    return((h^3) / 3)
  }
  w <- function(t) (t + cs)^(-1 / a)
  # g(hi) - g(τ) = ∫_τ^hi w
  g_diff <- function(tau) {
    if (isTRUE(all.equal(a, 1))) {
      log((hi + cs) / (tau + cs))
    } else {
      # ∫ (t+cs)^{-1/a} dt = a/(a-1) (t+cs)^{(a-1)/a}
      coef <- a / (a - 1)
      expn <- (a - 1) / a
      coef * ((hi + cs)^expn - (tau + cs)^expn)
    }
  }
  integrand <- function(tau) {
    ((g_diff(tau) / w(tau))^2)
  }
  out <- stats::integrate(integrand, lower = lo, upper = hi, rel.tol = 1e-8)
  as.numeric(out$value)
}
