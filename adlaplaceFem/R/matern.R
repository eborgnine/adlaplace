#' @include 000.R
#' @include matern-geometry.R
#' @include fem_bspline.R
#' @include fem_precision.R
NULL

#' Map Matern shape nu to B-spline degree / SPDE alpha (2D)
#'
#' @param shape Matern smoothness nu (`1` or `2`).
#' @return Integer degree (= alpha = nu + 1 in 2D).
#' @keywords internal
#' @noRd
matern_shape_degree <- function(shape) {
  shape <- as.integer(shape)[1L]
  if (is.na(shape) || !shape %in% c(1L, 2L)) {
    stop("shape must be 1 or 2 (Matern nu)", call. = FALSE)
  }
  shape + 1L
}

#' FEM Gram ingredients without an observation design matrix
#' @keywords internal
#' @noRd
matern_fem_grams <- function(knots, degree) {
  kn <- resolve_knots_list(knots, degree = degree)
  # Grams depend only on knots; use a dummy interior point for assembly
  mid_x <- mean(range(kn$x))
  mid_y <- mean(range(kn$y))
  fem <- fem_bspline_xy(mid_x, mid_y, knots = kn, degree = degree)
  fem$A <- NULL
  fem
}

#' Matern FEM spatial random field term
#'
#' Formula term for a 2D Matern GRF via tensor-product B-spline FEM, wired to
#' adlaplaceFem kernels `random_fem_ssq_2` / `random_fem_det_2` (shape nu = 1)
#' or `random_fem_ssq_3` / `random_fem_det_3` (nu = 2). The quadratic form is a
#' random shard; the log-determinant is a parameters companion via
#' [adlaplace::extra_ad_fun()].
#' The first argument is a **column name** in `data` (like [adlaplace::iid()]):
#' that column holds observation locations as a 2-column matrix, WKT text, or
#' HEX WKB points. Knot lines are supplied at construction; the observation
#' design matrix is built later in [design()] from `data`.
#'
#' In 2D, B-spline degree and SPDE order alpha equal `shape + 1`.
#'
#' @name matern-class
#' @aliases matern
#' @exportClass matern
#' @export
setClass(
  "matern",
  slots = c(
    fem = "ANY"
  ),
  contains = "model",
  prototype = prototype(
    ref_value = numeric(0),
    p.order = 2L,
    knots = numeric(0),
    type = factor("random", levels = adlaplace::.type_factor_levels),
    ad_fun = "random_fem_ssq_2",
    ad_kind = "random",
    package = "adlaplaceFem",
    fem = NULL
  )
)

#' @param x Name of the geometry / coordinate column in `data`.
#' @param knots Knot lines as `list(x = ..., y = ...)` or a terra
#'   `SpatRaster` (extent endpoints plus interior cell centers; see
#'   [knots_from_spatraster()]).
#' @param shape Matern smoothness nu (`1` or `2`). Default `1L` (SPDE alpha = 2,
#'   quadratic B-splines). Use `2` for nu = 2 / cubic B-splines.
#' @param init Initial values for `(range, sd)` on the natural scale, where
#'   `range = sqrt(8*nu)/kappa` is the practical range (same units as
#'   coordinates; matches geostatsp/INLA) and `sd` is the marginal SD of the
#'   Matern field. Internally converted to SPDE `(kappa, tau)`.
#' @param lower,upper Bounds for `(range, sd)`.
#' @param parscale Optimization scales for `(range, sd)`.
#' @param term A [`matern-class`] term.
#' @param data Model data frame (passed by adlaplace generics).
#' @return A `matern` term object.
#' @rdname matern-class
#' @export
matern <- function(x,
                   knots,
                   shape = 1L,
                   init = c(1, 1),
                   lower = c(1e-9, 1e-9),
                   upper = c(Inf, Inf),
                   parscale = c(1, 1)) {
  x <- adlaplace::strip_term_name(as.character(x))
  degree <- matern_shape_degree(shape)
  knots <- matern_knots(knots)
  fem <- matern_fem_grams(knots, degree = degree)

  ad_fun <- if (identical(degree, 3L)) {
    "random_fem_ssq_3"
  } else {
    "random_fem_ssq_2"
  }
  methods::new(
    "matern",
    term = x,
    label = paste(x, "matern", sep = "_"),
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = as.numeric(init),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    parscale = as.numeric(parscale),
    p.order = degree,
    ad_fun = ad_fun,
    ad_kind = "random",
    package = "adlaplaceFem",
    fem = fem
  )
}

#' @describeIn matern-class Companion parameters density for the FEM
#'   log-determinant shard.
#' @export
setMethod("extra_ad_fun", "matern", function(term) {
  sub("^random_fem_ssq_", "random_fem_det_", term@ad_fun, perl = TRUE)
})

#' @keywords internal
#' @noRd
ensure_matern_fem <- function(term) {
  if (is.null(term@fem)) {
    stop("matern term is missing cached FEM ingredients (@fem)", call. = FALSE)
  }
  term@fem
}

#' @describeIn matern-class Design matrix of B-spline basis evaluations at
#'   coordinates in `data[[term@term]]`.
#' @export
setMethod("design", "matern", function(term, data) {
  if (!term@term %in% names(data)) {
    stop(
      "matern geometry column '", term@term, "' not found in data",
      call. = FALSE
    )
  }
  fem <- ensure_matern_fem(term)
  xy <- matern_column_xy(data[[term@term]])
  if (nrow(xy) != nrow(data)) {
    stop(
      "matern geometry has ", nrow(xy), " rows but data has ", nrow(data),
      call. = FALSE
    )
  }
  A <- tensor_design(
    xy[, 1L],
    xy[, 2L],
    fem$knots$x,
    fem$knots$y,
    fem$degree
  )
  colnames(A) <- paste0(term@label, "_b", seq_len(ncol(A)))
  A
})

#' @describeIn matern-class FEM precision payload for `random_fem_*`.
#' @export
setMethod("precision", "matern", function(term, data) {
  fem <- ensure_matern_fem(term)
  fem_precision_payload(fem, alpha = term@p.order)
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

#' @describeIn matern-class Theta info for practical `range = sqrt(8*nu)/kappa`
#'   and field `sd` (log scale when `transform_theta` is TRUE).
#' @export
setMethod("theta_info", "matern", function(term) {
  data.frame(
    term = term@term,
    model = "matern",
    label = paste0(term@label, c("_log_range", "_log_sd")),
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
