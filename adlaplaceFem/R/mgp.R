#' Near-monotone mGP FEM term
#'
#' Cubic B-spline FEM for the near-monotone Gaussian process
#' \eqn{L_a f = f'' - \alpha f' = \sigma\xi} with
#' \eqn{\alpha(x) = -1/(a(x+c))}. The null space is spanned by
#' \eqn{\{1, m(x)\}} where \eqn{m} is the normalized Box-Cox transform of
#' \eqn{x+c}; the FEM weights drop the first two bases so
#' \eqn{f(x_0) = f'(x_0) = 0}, and the slope column is returned separately as
#' an [`mgp_boundary-class`] term (like [`adlaplace::iwp()`] returning `rpoly`).
#'
#' The prior density uses adlaplace's `random_mult` kernel with a fixed
#' precision `Q` (one SD theta). Use [adlaplace::psd()] /
#' [adlaplace::psd_prior()] for a PC prior on the predictive SD.
#'
#' FEM construction follows Ziang Zhang's BayesGP near-monotone FEM, with
#' Gauss-Legendre quadrature and [splines::splineDesign()] (no **fda**
#' dependency at runtime).
#'
#' @slot a Curvature parameter (`Inf` gives IWP2: \eqn{\alpha = 0}, \eqn{m(u)=u}).
#' @slot c Original-scale shift ensuring \eqn{x + c > 0}.
#' @slot shift Internal shift \eqn{c_s = x_0 + c}.
#' @slot region Length-2 numeric domain.
#' @slot k Number of cubic B-spline bases per side before dropping two.
#' @slot n_quad Gauss-Legendre nodes per knot span.
#' @slot fem Cached list from the internal FEM assembler.
#' @slot two_sided Logical; both sides of `ref_value` are active.
#'
#' @seealso [adlaplace::psd()], [adlaplace::psd_prior()], [mgp_boundary-class]
#' @exportClass mgp
#' @rdname mgp-class
setClass(
  "mgp",
  slots = list(
    a = "numeric",
    c = "numeric",
    shift = "numeric",
    region = "numeric",
    k = "integer",
    n_quad = "integer",
    fem = "list",
    two_sided = "logical"
  ),
  contains = "model_term",
  prototype = prototype(
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_mult",
    ad_kind = "random",
    package = "adlaplace",
    p.order = 2L,
    a = 2,
    c = 0,
    shift = 0,
    region = c(0, 1),
    k = 30L,
    n_quad = 10L,
    fem = list(),
    two_sided = FALSE
  )
)

#' Boundary (null-space slope) companion of an [`mgp`] term
#'
#' One-column design given by the normalized Box-Cox transform of
#' \eqn{x - x_0}. Default `sd = Inf` is a flat prior (no precision shard);
#' a finite `boundary_sd` uses `random_diagonal` with known SD.
#'
#' @slot sd Prior SD for the boundary coefficient (`Inf` = flat).
#' @slot a Curvature (same as the parent [`mgp-class`]).
#' @slot c Original-scale shift.
#' @slot shift Internal shift \eqn{c_s = x_0 + c}.
#' @exportClass mgp_boundary
#' @rdname mgp_boundary-class
setClass(
  "mgp_boundary",
  slots = list(
    sd = "numeric",
    a = "numeric",
    c = "numeric",
    shift = "numeric"
  ),
  contains = "model_term",
  prototype = prototype(
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random",
    package = "adlaplace",
    p.order = 1L,
    sd = Inf,
    a = 2,
    c = 0,
    shift = 0
  )
)

#' Near-monotone mGP constructor
#'
#' @param x Variable name.
#' @param a Curvature; `2` gives \eqn{\sqrt{x+c}}, `1` gives `\log`,
#'   `Inf` gives a linear base model (IWP2).
#' @param c Original-scale shift. Default `0` when `min(region) > 0`,
#'   otherwise `-min + 1e-6 * scale`.
#' @param k Cubic B-spline bases per side before dropping the first two
#'   (fda convention).
#' @param region Length-2 numeric domain (required).
#' @param ref_value Reference location: `"left"` (default), `"middle"`,
#'   `"right"`, or a numeric value inside `region`.
#' @param boundary_sd Prior SD for the boundary column (`Inf` = flat).
#' @param n_quad Gauss-Legendre nodes per knot span.
#' @param init,lower,upper,parscale,log SD hyperparameter controls.
#' @return Named list with the [`mgp-class`] term and its
#'   [`mgp_boundary-class`] companion.
#' @export
#' @rdname mgp-class
mgp <- function(
    x,
    a = 2,
    c = NULL,
    k = 30,
    region = NULL,
    ref_value = "left",
    boundary_sd = Inf,
    n_quad = 10,
    init = 0.02,
    lower = 1e-9,
    upper = Inf,
    parscale = 1,
    log = TRUE) {
  if (is.null(region)) {
    stop("region is required for mgp()", call. = FALSE)
  }
  region <- range(as.numeric(region), na.rm = TRUE)
  if (length(region) != 2L || any(!is.finite(region)) || diff(region) <= 0) {
    stop("region must be a finite length-2 interval with positive length", call. = FALSE)
  }
  if (length(a) != 1L || (!is.finite(a) && !isTRUE(all.equal(a, Inf))) || a == 0) {
    stop("a must be a single nonzero value, or Inf", call. = FALSE)
  }
  a <- as.numeric(a)
  k <- as.integer(k)
  n_quad <- as.integer(n_quad)
  if (length(boundary_sd) != 1L || !(is.finite(boundary_sd) || is.infinite(boundary_sd)) ||
    boundary_sd <= 0) {
    stop("boundary_sd must be a single positive value (or Inf)", call. = FALSE)
  }

  x <- adlaplace::strip_term_name(as.character(x))
  c_val <- mgp_shift(c, region)
  mgp_check_domain(region, c_val)
  x0 <- mgp_reference(region, ref_value)
  cs <- x0 + c_val
  if (!is.finite(cs) || cs <= 0) {
    stop("internal shift x0 + c must be positive", call. = FALSE)
  }

  fem <- mgp_fem(
    domain = region, x0 = x0, a = a, cs = cs, k = k, n_quad = n_quad
  )

  the_f <- stats::as.formula(paste0("~ 0 + ", x), env = new.env())
  mgp_name <- paste(c(x, "mgp"), collapse = "_")
  boundary_name <- paste(c(x, "mgp_boundary"), collapse = "_")

  mgp_term <- methods::new(
    "mgp",
    name = x,
    label = mgp_name,
    formula = the_f,
    p.order = 2L,
    ref_value = x0,
    knots = region,
    init = as.numeric(init)[1L],
    lower = as.numeric(lower)[1L],
    upper = as.numeric(upper)[1L],
    parscale = as.numeric(parscale)[1L],
    log = as.logical(log)[1L],
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_mult",
    ad_kind = "random",
    package = "adlaplace",
    a = a,
    c = c_val,
    shift = cs,
    region = region,
    k = k,
    n_quad = n_quad,
    fem = fem,
    two_sided = isTRUE(fem$two_sided)
  )

  boundary_term <- methods::new(
    "mgp_boundary",
    name = x,
    label = boundary_name,
    formula = the_f,
    p.order = 1L,
    ref_value = x0,
    knots = region,
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random",
    package = "adlaplace",
    sd = as.numeric(boundary_sd)[1L],
    a = a,
    c = c_val,
    shift = cs
  )

  result <- list()
  result[[mgp_name]] <- mgp_term
  result[[boundary_name]] <- boundary_term
  result
}

.mgp_basis_names <- function(term, n) {
  pad <- formatC(seq_len(n), width = max(1L, ceiling(log10(n + 1e-6))), flag = "0")
  paste0(term@name, "_mgp_k", pad)
}

#' @describeIn mgp-class Sparse B-spline design at `data[[x]]`.
#' @param term An `mgp` term.
#' @param data Model data frame.
#' @export
setMethod("design", "mgp", function(term, data) {
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  x <- as.numeric(data[[term@name]])
  if (any(!is.finite(x))) {
    stop("mgp covariate '", term@name, "' has non-finite values", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  if (any(x < term@region[[1L]] - tol | x > term@region[[2L]] + tol)) {
    stop(
      "mgp covariate '", term@name, "' has values outside region [",
      term@region[[1L]], ", ", term@region[[2L]], "]",
      call. = FALSE
    )
  }
  B <- mgp_fem_design(x, term@fem)
  colnames(B) <- .mgp_basis_names(term, ncol(B))
  B
})

#' @describeIn mgp-class Fixed `random_mult` precision payload.
#' @export
setMethod("precision", "mgp", function(term, data) {
  Q <- term@fem$Q
  list(
    Q = Q,
    log_det = term@fem$log_det,
    rank = as.integer(term@fem$rank)
  )
})

#' @describeIn mgp-class One SD theta.
#' @export
setMethod("theta_info", "mgp", function(term) {
  data.frame(
    term = term@name,
    model = "mgp",
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

#' @describeIn mgp-class No fixed effects.
#' @export
setMethod("beta_info", "mgp", function(term, data) {
  NULL
})

#' @describeIn mgp-class One row per FEM weight.
#' @export
setMethod("random_info", "mgp", function(term, data) {
  n <- as.integer(term@fem$rank)
  basis <- seq_len(n)
  result <- data.frame(
    term = term@name,
    model = "mgp",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result$gamma_label <- .mgp_basis_names(term, n)
  result
})

#' @describeIn mgp-class No companion log-det density (`random_mult` includes it).
#' @export
setMethod("extra_density", "mgp", function(term) {
  NULL
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "mgp", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  x <- methods::callNextMethod(x, map)
  x@label <- paste(c(x@name, "mgp"), collapse = "_")
  x
})

#' @describeIn mgp-class Predictive SD via Green's-function quadrature.
#' @param h Step size(s).
#' @param x Location(s) on the original scale (default `ref_value`).
#' @param sd SD hyperparameter scale (default 1).
#' @param ... Unused; for S4 generic compatibility.
#' @export
setMethod("psd", "mgp", function(term, h, x = NULL, sd = 1, ...) {
  if (is.null(x)) {
    x <- term@ref_value
  }
  x <- as.numeric(x)
  h <- as.numeric(h)
  n <- max(length(x), length(h))
  x <- rep_len(x, n)
  h <- rep_len(h, n)
  sd <- as.numeric(sd)[1L]
  if (!is.finite(sd) || sd < 0) {
    stop("sd must be a non-negative finite value", call. = FALSE)
  }
  if (any(x + term@c <= 0) || any(x + h + term@c <= 0)) {
    stop(
      "psd() for mgp requires x + c > 0 and x + h + c > 0 on the original scale",
      call. = FALSE
    )
  }
  u <- x - term@ref_value
  abs(sd) * mgp_psd_unit(u, h, a = term@a, cs = term@shift)
})

# --- mgp_boundary methods ----------------------------------------------------

#' @describeIn mgp_boundary-class Normalized Box-Cox design column.
#' @param term An `mgp_boundary` term.
#' @param data Model data frame.
#' @export
setMethod("design", "mgp_boundary", function(term, data) {
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  x <- as.numeric(data[[term@name]])
  u <- x - term@ref_value
  m <- mgp_boxcox(u, a = term@a, cs = term@shift)
  D <- matrix(m, ncol = 1L)
  colnames(D) <- paste0(term@name, "_mgp_boundary")
  D
})

#' @describeIn mgp_boundary-class Diagonal precision, or `NULL` when `sd = Inf`.
#' @export
setMethod("precision", "mgp_boundary", function(term, data) {
  if (is.infinite(term@sd)) {
    return(NULL)
  }
  precision_mat <- Matrix::Diagonal(n = 1L, x = 1 / (term@sd^2))
  dimnames(precision_mat) <- list(
    paste0(term@name, "_mgp_boundary")
  )[c(1, 1)]
  precision_mat
})

#' @describeIn mgp_boundary-class One random coefficient.
#' @export
setMethod("random_info", "mgp_boundary", function(term, data) {
  data.frame(
    term = term@name,
    model = "mgp_boundary",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = NA,
    order = 1L,
    gamma_label = paste0(term@label),
    stringsAsFactors = FALSE
  )
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "mgp_boundary", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  x <- methods::callNextMethod(x, map)
  x@label <- paste(c(x@name, "mgp_boundary"), collapse = "_")
  x
})
