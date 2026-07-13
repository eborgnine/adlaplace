#' @include 000.R
#' @include crds.R
#' @include grf_bspline.R
#' @include fem_precision.R
NULL

#' Matérn FEM spatial random field term
#'
#' Formula term for a 2D Matérn GRF via tensor-product B-spline FEM, wired to
#' adlaplace kernels `random_fem_2` (α = 2) or `random_fem_3` (α = 3).
#' Observation locations are taken from `x` via [crds()]; pass points with a
#' named argument (`matern(x = sites, ...)`) so formula parsing does not
#' stringify the object.
#'
#' @name matern-class
#' @aliases matern
#' @exportClass matern
#' @export
setClass(
  "matern",
  slots = c(
    points = "ANY",
    knots_list = "list",
    degree = "integer",
    alpha = "integer",
    fem = "ANY"
  ),
  contains = "model",
  prototype = prototype(
    ref_value = numeric(0),
    p.order = as.integer(0),
    knots = numeric(0),
    type = factor("random", levels = adlaplace:::.type_factor_levels),
    ad_fun = "random_fem_2",
    ad_kind = "random",
    package = "adlaplace",
    points = NULL,
    knots_list = list(),
    degree = 2L,
    alpha = 2L,
    fem = NULL
  )
)

#' @param x Observation points (matrix/data.frame or terra SpatVector points).
#' @param knots Knot-line list `list(x = ..., y = ...)` for [grf_bspline()].
#' @param degree B-spline degree (`>= 2`); use `3` for α = 3.
#' @param alpha SPDE order; default equals `degree` when `degree` is 2 or 3.
#' @param init Initial values for `(kappa, tau)` on the natural scale.
#' @param lower,upper Bounds for `(kappa, tau)`.
#' @param parscale Optimization scales for `(kappa, tau)`.
#' @return A `matern` term object.
#' @rdname matern-class
#' @export
matern <- function(x,
                   knots,
                   degree = 2L,
                   alpha = NULL,
                   init = c(1, 1),
                   lower = c(1e-9, 1e-9),
                   upper = c(Inf, Inf),
                   parscale = c(1, 1)) {
  degree <- as.integer(degree)[1L]
  if (is.na(degree) || degree < 2L) {
    stop("degree must be an integer >= 2", call. = FALSE)
  }
  if (is.null(alpha)) {
    alpha <- degree
  }
  alpha <- as.integer(alpha)[1L]
  if (!alpha %in% c(2L, 3L)) {
    stop("alpha must be 2 or 3", call. = FALSE)
  }
  if (alpha == 3L && degree < 3L) {
    stop("alpha = 3 requires degree >= 3", call. = FALSE)
  }
  if (!is.list(knots) || is.null(knots$x) || is.null(knots$y)) {
    stop("knots must be list(x = ..., y = ...)", call. = FALSE)
  }
  # Validate points early (also rejects polygons) and cache FEM ingredients
  xy <- crds(x)
  sites <- data.frame(x = xy[, 1L], y = xy[, 2L])
  fem <- grf_bspline(sites, knots, degree = degree)

  ad_fun <- if (identical(alpha, 3L)) "random_fem_3" else "random_fem_2"
  methods::new(
    "matern",
    term = "matern",
    label = "matern",
    formula = ~0,
    init = as.numeric(init),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    parscale = as.numeric(parscale),
    points = x,
    knots_list = knots,
    degree = degree,
    alpha = alpha,
    ad_fun = ad_fun,
    ad_kind = "random",
    package = "adlaplace",
    fem = fem
  )
}

#' @keywords internal
ensure_matern_fem <- function(term) {
  if (!is.null(term@fem)) {
    return(term@fem)
  }
  xy <- crds(term@points)
  sites <- data.frame(x = xy[, 1L], y = xy[, 2L])
  grf_bspline(sites, term@knots_list, degree = term@degree)
}

#' @describeIn matern-class Design matrix of B-spline basis evaluations at points.
#' @export
setMethod("design", "matern", function(term, data) {
  fem <- ensure_matern_fem(term)
  if (nrow(fem$A) != nrow(data)) {
    stop(
      "matern design has ", nrow(fem$A), " rows but data has ", nrow(data),
      call. = FALSE
    )
  }
  A <- fem$A
  colnames(A) <- paste0(term@label, "_b", seq_len(ncol(A)))
  A
})

#' @describeIn matern-class FEM precision payload for `random_fem_*`.
#' @export
setMethod("precision", "matern", function(term, data) {
  fem <- ensure_matern_fem(term)
  fem_precision_payload(fem, alpha = term@alpha)
})

#' @describeIn matern-class Random-effect metadata (one row per basis weight).
#' @export
setMethod("random_info", "matern", function(term, data) {
  fem <- ensure_matern_fem(term)
  n_basis <- nrow(fem$C)
  basis <- as.character(seq_len(n_basis))
  result <- data.frame(
    term = term@term,
    model = "matern",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = basis,
    order = NA,
    stringsAsFactors = FALSE
  )
  result$gamma_label <- paste0(term@label, "_b", basis)
  result
})

#' @describeIn matern-class Theta info for `kappa` and `tau` (log scale).
#' @export
setMethod("theta_info", "matern", function(term) {
  data.frame(
    term = term@term,
    model = "matern",
    label = paste0(term@label, c("_log_kappa", "_log_tau")),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    transform = TRUE,
    stringsAsFactors = FALSE
  )
})

#' @describeIn matern-class Beta info (none).
#' @export
setMethod("beta_info", "matern", function(term, data) {
  NULL
})
