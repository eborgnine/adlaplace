#' Gauss-Legendre nodes and weights on the reference interval
#' @keywords internal
gauss_legendre <- function(n) {
  n <- as.integer(n)
  if (n < 1L) stop("n must be positive")
  if (n == 1L) {
    return(list(nodes = 0, weights = 2))
  }
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  eg <- eigen(J, symmetric = TRUE)
  nodes <- eg$values
  weights <- 2 * eg$vectors[1, ]^2
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = as.numeric(weights[ord]))
}

#' Map GL nodes from the reference interval to an arbitrary interval
#' @keywords internal
gl_on_interval <- function(a, b, gl) {
  mid <- 0.5 * (a + b)
  half <- 0.5 * (b - a)
  list(x = mid + half * gl$nodes, w = half * gl$weights)
}

#' Unique increasing knot breakpoints for an open knot vector
#' @keywords internal
knot_spans <- function(knots) {
  uk <- sort(unique(as.numeric(knots)))
  if (length(uk) < 2L) {
    stop("need at least two distinct knots")
  }
  list(left = uk[-length(uk)], right = uk[-1L])
}

#' Evaluate B-spline basis or derivatives via splineDesign (sparse)
#' @keywords internal
bspline_eval <- function(knots, x, degree, derivs = 0L) {
  ord <- as.integer(degree) + 1L
  splines::splineDesign(
    knots = knots,
    x = x,
    ord = ord,
    derivs = as.integer(derivs),
    outer.ok = TRUE,
    sparse = TRUE
  )
}

#' 1D Gram matrix of B-spline derivatives
#'
#' Assembles \eqn{M^{pq}_{ik} = \int B_i^{(p)} B_k^{(q)}\,dx} by Gauss-Legendre
#' quadrature on each knot span.
#' @keywords internal
gram_1d <- function(knots, degree, deriv_a = 0L, deriv_b = 0L, n_quad = NULL) {
  degree <- as.integer(degree)
  deriv_a <- as.integer(deriv_a)
  deriv_b <- as.integer(deriv_b)
  if (deriv_a > degree || deriv_b > degree) {
    nb <- n_basis_knots(knots, degree)
    return(Matrix::Matrix(0, nb, nb, sparse = TRUE))
  }
  # Products of piecewise polynomials of degree (degree - deriv) need
  # ceil((2*degree - deriv_a - deriv_b)/2) + 1 nodes; use a safe default.
  if (is.null(n_quad)) {
    n_quad <- as.integer(degree - min(deriv_a, deriv_b) + 2L)
    n_quad <- max(n_quad, 2L)
  }
  gl <- gauss_legendre(n_quad)
  spans <- knot_spans(knots)
  nb <- n_basis_knots(knots, degree)
  acc <- Matrix::Matrix(0, nb, nb, sparse = TRUE)

  for (s in seq_along(spans$left)) {
    ab <- gl_on_interval(spans$left[s], spans$right[s], gl)
    Ba <- bspline_eval(knots, ab$x, degree, deriv_a)
    Bb <- if (deriv_b == deriv_a) Ba else bspline_eval(knots, ab$x, degree, deriv_b)
    Ba_w <- Matrix::Diagonal(x = ab$w) %*% Ba
    acc <- acc + Matrix::crossprod(Ba_w, Bb)
  }
  methods::as(Matrix::drop0(acc), "dgCMatrix")
}
