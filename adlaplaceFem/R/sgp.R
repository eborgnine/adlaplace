#' Seasonal Gaussian process FEM term
#'
#' Cubic seasonal B-spline FEM for the seasonal Gaussian process
#' \eqn{L g = g'' + \alpha^2 g = \sigma \xi} on \eqn{u = x - x_0}, with
#' \eqn{g(x_0) = g'(x_0) = 0}. The null space
#' \eqn{\{\cos(i \alpha u), \sin(i \alpha u)\}_{i=1}^{m}} is returned separately
#' as an [`sgp_boundary-class`] term. Harmonics are independent and share one
#' SD. \eqn{x_0} is the left endpoint of the knots.
#'
#' Pass original-scale knot positions (including the domain ends), or an
#' integer `k` and `region`. `k` matches BayesGP's equally spaced cubic basis
#' (`nbasis = k`, `norder = 4`), which uses `k - 2` breakpoints.
#'
#' The prior density uses adlaplace's `random_mult` kernel with a fixed
#' precision `Q` (one SD theta). Use [adlaplace::psd()] /
#' [adlaplace::psd_prior()] for a PC prior on the predictive SD. The predictive
#' SD does not depend on location. For `m > 1` it is the SD of the sum of
#' independent harmonics.
#'
#' FEM construction follows Ziang Zhang's BayesGP seasonal B-splines
#' (arXiv:2305.09914), with Gauss-Legendre quadrature and
#' [splines::splineDesign()] (no **fda** dependency at runtime).
#'
#' @slot a Frequency \eqn{\alpha = 2\pi/\mathrm{period}}.
#' @slot m Number of harmonics. Frequencies are \eqn{\alpha, 2\alpha, \ldots, m\alpha}.
#' @slot n_quad Gauss-Legendre nodes per knot span.
#' @slot fem Cached list from the internal FEM assembler.
#'
#' @seealso [adlaplace::psd()], [adlaplace::psd_prior()], [sgp_boundary-class]
#' @exportClass sgp
#' @rdname sgp-class
setClass(
  "sgp",
  slots = list(
    a = "numeric",
    m = "integer",
    n_quad = "integer",
    fem = "list"
  ),
  contains = "model_term",
  prototype = prototype(
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_mult",
    ad_kind = "random",
    package = "adlaplace",
    p.order = 2L,
    a = NA_real_,
    m = 1L,
    n_quad = 16L,
    fem = list()
  )
)

#' Sin/cos boundary companion of an [`sgp`] term
#'
#' One cosine and one sine column per harmonic, in harmonic order:
#' \eqn{\cos(i \alpha u), \sin(i \alpha u)}. Default `sd = Inf` is a flat prior
#' (no precision shard). A finite `boundary_sd` uses `random_diagonal` with
#' that known SD on every column.
#'
#' @slot sd Prior SD for the boundary coefficients (`Inf` = flat).
#' @slot a Frequency (same as the parent [`sgp-class`]).
#' @slot m Number of harmonics (same as the parent).
#' @exportClass sgp_boundary
#' @rdname sgp_boundary-class
setClass(
  "sgp_boundary",
  slots = list(
    sd = "numeric",
    a = "numeric",
    m = "integer"
  ),
  contains = "model_term",
  prototype = prototype(
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random",
    package = "adlaplace",
    p.order = 1L,
    sd = Inf,
    a = NA_real_,
    m = 1L
  )
)

#' @noRd
.sgp_frequency <- function(period, a) {
  has_period <- !is.null(period)
  has_a <- !is.null(a)
  if (has_period == has_a) {
    stop("supply exactly one of period or a", call. = FALSE)
  }
  if (has_period) {
    period <- as.numeric(period)
    if (length(period) != 1L || !is.finite(period) || period <= 0) {
      stop("period must be a single positive value", call. = FALSE)
    }
    return(2 * pi / period)
  }
  a <- as.numeric(a)
  if (length(a) != 1L || !is.finite(a) || a <= 0) {
    stop("a must be a single positive value", call. = FALSE)
  }
  a
}

#' @noRd
.sgp_resolve_knots <- function(knots, k, region) {
  has_knots <- !is.null(knots)
  has_k <- !is.null(k)
  if (has_knots && has_k) {
    stop("supply knots or k, not both", call. = FALSE)
  }
  if (!has_knots && !has_k) {
    stop("supply knots, or both k and region", call. = FALSE)
  }
  if (has_k) {
    if (is.null(region)) {
      stop("region is required when k is supplied", call. = FALSE)
    }
    return(sgp_axis_from_k(region, k))
  }
  knots
}

#' Seasonal GP constructor
#'
#' @param x Variable name.
#' @param period Period of the first harmonic. Supply this or `a`, not both.
#'   The frequency is \eqn{\alpha = 2\pi/\mathrm{period}}.
#' @param a Frequency of the first harmonic. Supply this or `period`, not both.
#' @param knots Breakpoints on the original scale, including the domain ends.
#'   Supply `knots`, or `k` and `region`.
#' @param k Number of cubic B-spline bases before the boundary drop, in the
#'   BayesGP sense (`nbasis = k`). Requires `region`. Must be at least 4.
#' @param region Length-2 interval used with `k`.
#' @param m Number of harmonics (default 1).
#' @param boundary_sd Prior SD for the sin/cos columns (`Inf` = flat).
#'   BayesGP's precision 0.001 is `boundary_sd = sqrt(1000)`.
#' @param n_quad Gauss-Legendre nodes per knot span.
#' @param init,lower,upper,parscale,log SD hyperparameter controls.
#' @return Named list with the [`sgp-class`] term and its
#'   [`sgp_boundary-class`] companion.
#' @examples
#' sgp("year", period = 1, k = 10, region = c(0, 10))
#' @export
#' @rdname sgp-class
sgp <- function(
    x,
    period = NULL,
    a = NULL,
    knots = NULL,
    k = NULL,
    region = NULL,
    m = 1L,
    boundary_sd = Inf,
    n_quad = 16L,
    init = 0.1,
    lower = 1e-9,
    upper = Inf,
    parscale = 1,
    log = TRUE) {
  if (length(boundary_sd) != 1L ||
    !(is.finite(boundary_sd) || is.infinite(boundary_sd)) ||
    boundary_sd <= 0) {
    stop("boundary_sd must be a single positive value (or Inf)", call. = FALSE)
  }
  m <- as.integer(m)[1L]
  if (!is.finite(m) || m < 1L) {
    stop("m must be an integer >= 1", call. = FALSE)
  }
  n_quad <- as.integer(n_quad)[1L]
  if (!is.finite(n_quad) || n_quad < 2L) {
    stop("n_quad must be an integer >= 2", call. = FALSE)
  }
  alpha <- .sgp_frequency(period, a)
  axis <- .sgp_resolve_knots(knots, k, region)
  x <- adlaplace::strip_term_name(as.character(x))
  if (length(x) != 1L || !nzchar(x)) {
    stop("x must be a single variable name", call. = FALSE)
  }
  fem <- sgp_fem(knots = axis, a = alpha, m = m, n_quad = n_quad)

  the_f <- stats::as.formula(paste0("~ 0 + ", x), env = new.env())
  sgp_name <- paste(c(x, "sgp"), collapse = "_")
  boundary_name <- paste(c(x, "sgp_boundary"), collapse = "_")

  sgp_term <- methods::new(
    "sgp",
    name = x,
    label = sgp_name,
    formula = the_f,
    p.order = 2L,
    ref_value = fem$x0,
    knots = fem$axis,
    init = as.numeric(init)[1L],
    lower = as.numeric(lower)[1L],
    upper = as.numeric(upper)[1L],
    parscale = as.numeric(parscale)[1L],
    log = as.logical(log)[1L],
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_mult",
    ad_kind = "random",
    package = "adlaplace",
    a = alpha,
    m = m,
    n_quad = n_quad,
    fem = fem
  )

  boundary_term <- methods::new(
    "sgp_boundary",
    name = x,
    label = boundary_name,
    formula = the_f,
    p.order = 1L,
    ref_value = fem$x0,
    knots = fem$axis,
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random",
    package = "adlaplace",
    sd = as.numeric(boundary_sd)[1L],
    a = alpha,
    m = m
  )

  result <- list()
  result[[sgp_name]] <- sgp_term
  result[[boundary_name]] <- boundary_term
  result
}

.sgp_basis_names <- function(term, n) {
  pad <- formatC(
    seq_len(n),
    width = max(1L, ceiling(log10(n + 1e-6))),
    flag = "0"
  )
  paste0(term@name, "_sgp_k", pad)
}

.sgp_boundary_labels <- function(term) {
  m <- as.integer(term@m)
  labs <- character(2L * m)
  for (i in seq_len(m)) {
    labs[[2L * i - 1L]] <- paste0(term@label, "_c", i)
    labs[[2L * i]] <- paste0(term@label, "_s", i)
  }
  labs
}

.sgp_check_x <- function(term, data, what = "sgp") {
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  x <- as.numeric(data[[term@name]])
  if (any(!is.finite(x))) {
    stop(what, " covariate '", term@name, "' has non-finite values", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  kn <- range(term@knots)
  if (any(x < kn[[1L]] - tol | x > kn[[2L]] + tol)) {
    stop(
      what, " covariate '", term@name, "' has values outside knots [",
      kn[[1L]], ", ", kn[[2L]], "]",
      call. = FALSE
    )
  }
  x
}

#' @describeIn sgp-class Sparse seasonal B-spline design at `data[[x]]`.
#' @param term An `sgp` term.
#' @param data Model data frame.
#' @export
setMethod("design", "sgp", function(term, data) {
  .sgp_check_x(term, data, what = "sgp")
  B <- sgp_fem_design(data[[term@name]], term@fem)
  colnames(B) <- .sgp_basis_names(term, ncol(B))
  B
})

#' @describeIn sgp-class Fixed `random_mult` precision payload.
#' @export
setMethod("precision", "sgp", function(term, data) {
  list(
    Q = term@fem$Q,
    log_det = term@fem$log_det,
    rank = as.integer(term@fem$rank)
  )
})

#' @describeIn sgp-class One SD theta.
#' @export
setMethod("theta_info", "sgp", function(term) {
  data.frame(
    term = term@name,
    model = "sgp",
    label = term@label,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    model_role = term@model_role,
    log = term@log,
    stringsAsFactors = FALSE
  )
})

#' @describeIn sgp-class No fixed effects.
#' @export
setMethod("beta_info", "sgp", function(term, data) {
  NULL
})

#' @describeIn sgp-class One row per FEM weight.
#' @export
setMethod("random_info", "sgp", function(term, data) {
  n <- as.integer(term@fem$rank)
  result <- data.frame(
    term = term@name,
    model = "sgp",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = seq_len(n),
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result$gamma_label <- .sgp_basis_names(term, n)
  result
})

#' @describeIn sgp-class No companion log-det density (`random_mult` includes it).
#' @export
setMethod("extra_density", "sgp", function(term) {
  NULL
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "sgp", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  x <- methods::callNextMethod(x, map)
  x@label <- paste(c(x@name, "sgp"), collapse = "_")
  x
})

#' @describeIn sgp-class Predictive SD of the seasonal GP (location-free).
#' @param h Step size(s). The absolute value is used.
#' @param x Ignored. The sGP predictive SD does not depend on location.
#' @param sd SD hyperparameter scale (default 1).
#' @param ... Unused; for S4 generic compatibility.
#' @export
setMethod("psd", "sgp", function(term, h, x = NULL, sd = 1, ...) {
  h <- as.numeric(h)
  if (!length(h) || any(!is.finite(h))) {
    stop("h must be finite", call. = FALSE)
  }
  sd <- as.numeric(sd)[1L]
  if (!is.finite(sd) || sd < 0) {
    stop("sd must be a non-negative finite value", call. = FALSE)
  }
  abs(sd) * sgp_psd_unit(h, a = term@a, m = term@m)
})

# --- sgp_boundary methods ----------------------------------------------------

#' @describeIn sgp_boundary-class Cosine and sine columns, one pair per harmonic.
#' @param term An `sgp_boundary` term.
#' @param data Model data frame.
#' @export
setMethod("design", "sgp_boundary", function(term, data) {
  x <- .sgp_check_x(term, data, what = "sgp boundary")
  u <- x - term@ref_value
  m <- as.integer(term@m)
  D <- matrix(0, nrow = length(u), ncol = 2L * m)
  for (i in seq_len(m)) {
    D[, 2L * i - 1L] <- cos(i * term@a * u)
    D[, 2L * i] <- sin(i * term@a * u)
  }
  colnames(D) <- .sgp_boundary_labels(term)
  D
})

#' @describeIn sgp_boundary-class Diagonal precision, or `NULL` when `sd = Inf`.
#' @export
setMethod("precision", "sgp_boundary", function(term, data) {
  if (is.infinite(term@sd)) {
    return(NULL)
  }
  labs <- .sgp_boundary_labels(term)
  precision_mat <- Matrix::Diagonal(n = length(labs), x = 1 / (term@sd^2))
  dimnames(precision_mat) <- list(labs)[c(1L, 1L)]
  precision_mat
})

#' @describeIn sgp_boundary-class One random coefficient per sin/cos column.
#' @export
setMethod("random_info", "sgp_boundary", function(term, data) {
  labs <- .sgp_boundary_labels(term)
  data.frame(
    term = term@name,
    model = "sgp_boundary",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = seq_along(labs),
    order = 1L,
    gamma_label = labs,
    stringsAsFactors = FALSE
  )
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "sgp_boundary", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  x <- methods::callNextMethod(x, map)
  x@label <- paste(c(x@name, "sgp_boundary"), collapse = "_")
  x
})
