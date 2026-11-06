#' Seasonal B-spline FEM for the seasonal Gaussian process
#'
#' Ported from Ziang Zhang's BayesGP `Compute_Q_sB` / `Compute_B_sB`
#' (https://github.com/Bayes-GP/BayesGP, master). The precision is the
#' least-squares FEM matrix \eqn{Q_{ij} = \int (L\phi_i)(L\phi_j)\,du} for
#' \eqn{L = d^2/du^2 + \omega^2} on the seasonal B-spline basis
#' \eqn{(b_j\cos(\omega u), b_j\sin(\omega u), b_j)}, which is the expanded
#' form \eqn{\omega^4 G + C + \omega^2(M+M^\top)}.
#' Gauss-Legendre quadrature replaces the Riemann sum, and
#' [splines::splineDesign()] replaces **fda**.
#'
#' @name sgp-fem
#' @keywords internal
#' @noRd
NULL

#' Open cubic basis on the shifted interval \code{[0, L]}, drop the first two.
#'
#' `knots` are original-scale breakpoints, including the domain ends.
#' Dropping bases 1 and 2 enforces \eqn{f(0) = f'(0) = 0}.
#' fda `nbasis = k` corresponds to `k - 2` breakpoints; see [sgp_axis_from_k()].
#' @noRd
sgp_basis <- function(knots) {
  axis <- sort(unique(as.numeric(knots)))
  if (length(axis) < 2L || any(!is.finite(axis))) {
    stop("knots must have at least two finite positions", call. = FALSE)
  }
  x0 <- axis[[1L]]
  local <- axis - x0
  tol <- sqrt(.Machine$double.eps)
  local[abs(local) < tol] <- 0
  if (local[[length(local)]] <= tol) {
    stop("sgp knot range must have positive length", call. = FALSE)
  }
  open <- axis_to_open_knots(local, degree = 3L)
  nb <- n_basis_knots(open, degree = 3L)
  if (nb < 3L) {
    stop(
      "sgp needs enough knots for at least one basis after dropping the boundary",
      call. = FALSE
    )
  }
  list(
    knots = open,
    degree = 3L,
    drop = c(1L, 2L),
    n_kept = nb - 2L,
    axis = axis,
    x0 = x0
  )
}

#' Equally spaced breakpoints for fda `nbasis = k`, `norder = 4`.
#'
#' That basis uses `k - 2` breaks, including the ends of `region`.
#' @noRd
sgp_axis_from_k <- function(region, k) {
  k <- as.integer(k)[1L]
  if (!is.finite(k) || k < 4L) {
    stop("k must be an integer >= 4", call. = FALSE)
  }
  region <- as.numeric(region)
  if (length(region) != 2L || any(!is.finite(region)) || region[[2L]] <= region[[1L]]) {
    stop("region must be a finite length-2 interval with max > min", call. = FALSE)
  }
  seq(region[[1L]], region[[2L]], length.out = k - 2L)
}

#' One-harmonic precision on a cached [sgp_basis()] result.
#'
#' Column order matches BayesGP: cos block, sin block, plain block.
#' @noRd
sgp_harmonic_precision <- function(basis, omega, n_quad = 16L) {
  omega <- as.numeric(omega)[1L]
  if (!is.finite(omega) || omega <= 0) {
    stop("frequency omega must be positive", call. = FALSE)
  }
  n_quad <- as.integer(n_quad)[1L]
  if (!is.finite(n_quad) || n_quad < 2L) {
    stop("n_quad must be an integer >= 2", call. = FALSE)
  }
  degree <- basis$degree
  gl <- gauss_legendre(n_quad)
  spans <- knot_spans(basis$knots)
  nb <- n_basis_knots(basis$knots, degree)
  keep <- setdiff(seq_len(nb), basis$drop)
  n <- 3L * length(keep)
  Q <- Matrix::Matrix(0, n, n, sparse = TRUE)

  for (s in seq_along(spans$left)) {
    ab <- gl_on_interval(spans$left[s], spans$right[s], gl)
    B0 <- bspline_eval(basis$knots, ab$x, degree, 0L)[, keep, drop = FALSE]
    B1 <- bspline_eval(basis$knots, ab$x, degree, 1L)[, keep, drop = FALSE]
    B2 <- bspline_eval(basis$knots, ab$x, degree, 2L)[, keep, drop = FALSE]
    cu <- cos(omega * ab$x)
    su <- sin(omega * ab$x)
    # L(b cos) = b'' cos - 2 ω b' sin
    # L(b sin) = b'' sin + 2 ω b' cos
    # L(b)     = b'' + ω^2 b
    Lc <- Matrix::Diagonal(x = cu) %*% B2 -
      (2 * omega) * (Matrix::Diagonal(x = su) %*% B1)
    Ls <- Matrix::Diagonal(x = su) %*% B2 +
      (2 * omega) * (Matrix::Diagonal(x = cu) %*% B1)
    Lp <- B2 + (omega^2) * B0
    Lphi <- cbind(Lc, Ls, Lp)
    Q <- Q + Matrix::crossprod(Matrix::Diagonal(x = ab$w) %*% Lphi, Lphi)
  }
  as_dgc_matrix(Matrix::drop0(Matrix::forceSymmetric(Q)))
}

#' Block-diagonal seasonal FEM (one block per harmonic) with log-determinant.
#' @noRd
sgp_fem <- function(knots, a, m = 1L, n_quad = 16L) {
  basis <- sgp_basis(knots)
  a <- as.numeric(a)[1L]
  m <- as.integer(m)[1L]
  if (!is.finite(a) || a <= 0) {
    stop("a must be positive", call. = FALSE)
  }
  if (!is.finite(m) || m < 1L) {
    stop("m must be an integer >= 1", call. = FALSE)
  }
  n_quad <- as.integer(n_quad)[1L]
  blocks <- lapply(seq_len(m), function(i) {
    sgp_harmonic_precision(basis, omega = i * a, n_quad = n_quad)
  })
  Q <- if (length(blocks) == 1L) {
    blocks[[1L]]
  } else {
    as_dgc_matrix(Matrix::bdiag(blocks))
  }
  Q <- Matrix::forceSymmetric(Q, uplo = "U")
  # B-splines that nearly reproduce sin and cos are a numerical kernel of
  # L = d^2/du^2 + ω^2. Those eigenvalues sit about 1e-12 below the rest, so
  # the Laplace Hessian loses positive definiteness and conjugate gradients
  # stall. A diagonal shift of 1e-4 of the largest diagonal entry is far
  # under the 1e-3 relative tolerance against BayesGP's Q. The shifted
  # directions have almost no design weight, so the seasonal curve is unchanged.
  d <- Matrix::diag(Q)
  scale <- max(d)
  if (is.finite(scale) && scale > 0) {
    Matrix::diag(Q) <- d + 1e-4 * scale
    Q <- Matrix::forceSymmetric(Q, uplo = "U")
  }
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  if (!is.finite(log_det)) {
    stop(
      "sgp precision Q is not positive definite (log_det not finite)",
      call. = FALSE
    )
  }
  c(
    basis,
    list(
      Q = Q,
      log_det = log_det,
      rank = ncol(Q),
      a = a,
      m = m,
      n_quad = n_quad
    )
  )
}

#' Seasonal B-spline design at original-scale `x`.
#'
#' Columns are harmonic-major, and within a harmonic: cos, sin, plain.
#' @noRd
sgp_fem_design <- function(x, fem) {
  x <- as.numeric(x)
  u <- x - fem$x0
  B <- bspline_eval(fem$knots, u, fem$degree, 0L)
  B <- methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix")
  keep <- setdiff(seq_len(ncol(B)), fem$drop)
  B <- B[, keep, drop = FALSE]
  blocks <- vector("list", fem$m)
  for (i in seq_len(fem$m)) {
    omega <- i * fem$a
    blocks[[i]] <- cbind(
      Matrix::Diagonal(x = cos(omega * u)) %*% B,
      Matrix::Diagonal(x = sin(omega * u)) %*% B,
      B
    )
  }
  out <- if (fem$m == 1L) blocks[[1L]] else do.call(cbind, blocks)
  as_dgc_matrix(out)
}

#' Squared unit-SD factor \eqn{(1/\omega^2)(h/2 - \sin(2\omega h)/(4\omega))}.
#' @noRd
sgp_step_var <- function(h, omega) {
  h <- abs(as.numeric(h))
  inside <- h / 2 - sin(2 * omega * h) / (4 * omega)
  pmax(inside, 0) / (omega^2)
}

#' Unit predictive SD of m independent harmonics with a shared SD of 1.
#'
#' For `m = 1` this is BayesGP `compute_d_step_sgpsd`. For `m > 1` the
#' harmonics are independent, so the factors combine in variance.
#' @noRd
sgp_psd_unit <- function(h, a, m = 1L) {
  h <- as.numeric(h)
  a <- as.numeric(a)[1L]
  m <- as.integer(m)[1L]
  acc <- rep(0, length(h))
  for (i in seq_len(m)) {
    acc <- acc + sgp_step_var(h, i * a)
  }
  sqrt(pmax(acc, 0))
}
