#' @noRd
log_mvn_density <- function(x, center, chol_cov) {
  delta <- x - center
  quad <- sum((backsolve(chol_cov, delta, transpose = TRUE))^2)
  log_det <- 2 * sum(log(diag(chol_cov)))
  ngamma <- length(center)
  -0.5 * (ngamma * log(2 * pi) + log_det + quad)
}

#' @noRd
full_parameters_from_gamma <- function(gamma_row, parameters, n_beta, n_theta) {
  c(
    parameters[seq_len(n_beta)],
    gamma_row,
    parameters[n_beta + seq_len(n_theta)]
  )
}

#' Posterior mean of random effects via multivariate quadrature
#'
#' Approximates \eqn{E[\gamma \mid y]} by importance-weighted quadrature over
#' \eqn{\gamma}, using a Gaussian proposal \eqn{N(\hat\gamma, H^{-1})} (Laplace
#' covariance) and correcting with \code{\link{joint_log_dens}}.
#'
#' @param parameters Outer \eqn{(\beta, \theta)} at the Laplace fit.
#' @param mode Inner mode \eqn{\hat\gamma} (e.g. \code{laplace$opt$solution}).
#' @param cov \eqn{H^{-1}} for the inner block (e.g. \code{laplace$extra$hessian$H_inv}).
#' @param ad_fun \code{ad_fun} object for \code{joint_log_dens}.
#' @param n Quadrature level passed to \pkg{mvQuad} \code{createNIGrid}.
#' @param type mvQuad rule (default \code{"nHN"}).
#' @param ndConstruction mvQuad grid construction (default \code{"sparse"}).
#' @param dec.type Cholesky decomposition flag for \code{mvQuad::rescale}.
#'
#' @return Numeric vector of length \code{n_gamma}.
#'
#' @details
#' Builds a sparse \code{nHN} grid with \pkg{mvQuad}, rescales to the Laplace
#' proposal, then forms
#' \eqn{\sum_i w_i \gamma_i \exp(\ell(\gamma_i) - \log q(\gamma_i)) /
#' \sum_i w_i \exp(\ell(\gamma_i) - \log q(\gamma_i))}
#' with log-sum-exp stabilization.
#'
#' @seealso \code{\link{log_lik_laplace}}, \code{\link{inner_opt}},
#'   \code{\link{joint_log_dens}}, \code{\link{rmvnldl}}
#' @export
mean_mvquad <- function(parameters,
                        mode,
                        cov,
                        ad_fun,
                        n = 3L,
                        type = "nHN",
                        ndConstruction = "sparse",
                        dec.type = 2L) {
  if (!requireNamespace("mvQuad", quietly = TRUE)) {
    stop("Package 'mvQuad' is required for mean_mvquad()", call. = FALSE)
  }
  if (!is(ad_fun, "ad_fun")) {
    stop("ad_fun must be an ad_fun object", call. = FALSE)
  }
  if (missing(parameters) || is.null(parameters)) {
    stop("parameters is required (outer beta and theta)", call. = FALSE)
  }
  if (missing(mode) || is.null(mode)) {
    stop("mode is required (inner gamma at Laplace mode)", call. = FALSE)
  }
  if (missing(cov) || is.null(cov)) {
    stop("cov is required (inner Laplace covariance H_inv)", call. = FALSE)
  }

  sz <- ad_fun@sizes
  n_beta <- as.integer(sz["beta"])
  n_gamma <- as.integer(sz["gamma"])
  n_theta <- as.integer(sz["theta"])

  mode <- as.numeric(mode)
  parameters <- as.numeric(parameters)

  if (length(mode) != n_gamma) {
    stop(
      "length(mode) (", length(mode), ") must equal ad_fun@sizes['gamma'] (",
      n_gamma, ")",
      call. = FALSE
    )
  }
  if (length(parameters) != n_beta + n_theta) {
    stop(
      "length(parameters) (", length(parameters),
      ") must equal n_beta + n_theta (", n_beta + n_theta, ")",
      call. = FALSE
    )
  }

  # H_inv is a symmetric CSC matrix but often stores only the lower triangle;
  # as.matrix() would zero the upper part and fail mvQuad::rescale's symmetry check.
  cov_mat <- if (inherits(cov, "Matrix")) {
    as.matrix(Matrix::forceSymmetric(cov, uplo = "L"))
  } else {
    cm <- as.matrix(cov)
    if (!isSymmetric(cm, check.attributes = FALSE)) {
      cm <- (cm + t(cm)) / 2
    }
    cm
  }
  if (!is.matrix(cov_mat) || nrow(cov_mat) != n_gamma || ncol(cov_mat) != n_gamma) {
    stop(
      "cov must be an n_gamma x n_gamma matrix (n_gamma = ", n_gamma, ")",
      call. = FALSE
    )
  }

  n <- as.integer(n)
  if (n < 1L) {
    stop("n (quadrature level) must be >= 1", call. = FALSE)
  }

  # mvQuad uses data.table::rbindlist with mismatched column names for sparse
  # Smolyak grids; harmless messages only (see data.table v1.12.2 news item 5).
  grid <- suppressMessages({
    g <- mvQuad::createNIGrid(
      dim = n_gamma,
      type = type,
      level = n,
      ndConstruction = ndConstruction
    )
    mvQuad::rescale(g, m = mode, C = cov_mat, dec.type = dec.type)
    g
  })

  nodes <- mvQuad::getNodes(grid)
  weights <- as.numeric(mvQuad::getWeights(grid))
  if (is.vector(nodes)) {
    nodes <- matrix(nodes, nrow = 1L)
  }
  n_pts <- nrow(nodes)
  if (length(weights) != n_pts) {
    stop("mvQuad grid: length(weights) must equal nrow(nodes)", call. = FALSE)
  }

  chol_cov <- chol(cov_mat)
  log_ratio <- vapply(
    seq_len(n_pts),
    function(i) {
      gamma_row <- nodes[i, , drop = TRUE]
      full <- full_parameters_from_gamma(
        gamma_row,
        parameters,
        n_beta,
        n_theta
      )
      log_p <- joint_log_dens(ad_fun, full, negative = FALSE)
      log_q <- log_mvn_density(gamma_row, mode, chol_cov)
      log_p - log_q
    },
    numeric(1)
  )

  log_ratio <- log_ratio - max(log_ratio)
  rw <- weights * exp(log_ratio)
  denom <- sum(rw)
  if (!is.finite(denom) || denom <= 0) {
    stop("mean_mvquad: non-positive or non-finite importance weights", call. = FALSE)
  }

  colSums(nodes * rw) / denom
}
