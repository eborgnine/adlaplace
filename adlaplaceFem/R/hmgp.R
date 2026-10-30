#' @include mgp.R
NULL

#' Grouped near-monotone mGP
#'
#' One common mGP curve plus an independent mGP deviation per level of
#' `by`, with one shared deviation SD, and an iid boundary slope per level
#' with its own SD. The deviation curve is the same operator as [mgp()], so
#' [adlaplace::psd()] is inherited.
#'
#' @slot by Grouping variable. One deviation curve and one slope per level.
#' @exportClass hmgp
#' @rdname hmgp-class
setClass(
  "hmgp",
  slots = list(by = "by_group"),
  contains = "mgp",
  prototype = prototype(by = adlaplace::by_group())
)

#' Per-group boundary slope of an [hmgp-class] term
#'
#' One Box-Cox column per group, with a shared SD theta (identity precision,
#' as in [adlaplace::iid()]).
#'
#' @slot by Grouping variable, shared with the parent [hmgp-class] term.
#' @exportClass hmgp_boundary
#' @rdname hmgp_boundary-class
setClass(
  "hmgp_boundary",
  slots = list(by = "by_group"),
  contains = "mgp_boundary",
  prototype = prototype(
    by = adlaplace::by_group(),
    sd = NA_real_
  )
)

#' @param x Variable name.
#' @param a Curvature; see [mgp()].
#' @param knots Knot locations on the original scale. Strictly positive.
#' @param ref_value Reference location, as in [mgp()].
#' @param by Grouping column. Levels share the deviation SD and the slope SD.
#' @param levels Optional group levels, in block order. When omitted, levels
#'   are taken from `data` on the first design or precision call. Pass them
#'   explicitly so a prediction grid uses the same column order as a fitted
#'   model.
#' @param n_quad Gauss-Legendre nodes per knot span.
#' @param init,lower,upper,parscale,log SD controls for the deviation curve.
#'   Reused for the global [mgp()] term when `include_global` is `TRUE`.
#' @param slope_init,slope_lower,slope_upper,slope_parscale,slope_log SD
#'   controls for the per-group boundary slope.
#' @param include_global When `TRUE` (default), also return the common
#'   [mgp()] curve and its flat boundary column.
#' @param boundary_sd Prior SD of the global boundary column (`Inf` = flat).
#'   Used only when `include_global` is `TRUE`.
#' @return Named list. With `include_global = TRUE`: `<x>_mgp`,
#'   `<x>_mgp_boundary`, `<x>_hmgp`, `<x>_hslope_mgp`. Otherwise the last two.
#' @examples
#' knots <- seq(8, 100, by = 4)
#' hmgp("pm", a = 2, knots = knots, by = "region", levels = c("BC", "ON"))
#' @export
#' @rdname hmgp-class
hmgp <- function(
    x,
    a = 2,
    knots,
    ref_value = "left",
    by,
    levels = NULL,
    n_quad = 10,
    init = 0.02,
    lower = 1e-9,
    upper = Inf,
    parscale = 1,
    log = TRUE,
    slope_init = 0.02,
    slope_lower = 1e-9,
    slope_upper = Inf,
    slope_parscale = 1,
    slope_log = TRUE,
    include_global = TRUE,
    boundary_sd = Inf) {
  if (missing(by) || length(by) != 1L || !nzchar(as.character(by))) {
    stop("by is required for hmgp()", call. = FALSE)
  }
  .hmgp_check_scalar(init, "init")
  .hmgp_check_scalar(lower, "lower", positive = FALSE)
  .hmgp_check_scalar(upper, "upper", positive = FALSE, allow_inf = TRUE)
  .hmgp_check_scalar(parscale, "parscale")
  .hmgp_check_scalar(slope_init, "slope_init")
  .hmgp_check_scalar(slope_lower, "slope_lower", positive = FALSE)
  .hmgp_check_scalar(slope_upper, "slope_upper", positive = FALSE, allow_inf = TRUE)
  .hmgp_check_scalar(slope_parscale, "slope_parscale")
  if (length(log) != 1L || length(slope_log) != 1L) {
    stop("log and slope_log must be length 1", call. = FALSE)
  }
  if (length(include_global) != 1L) {
    stop("include_global must be a single logical value", call. = FALSE)
  }

  setup <- .mgp_setup(x, a, knots, ref_value, n_quad)
  by_obj <- adlaplace::by_group(term = as.character(by), levels = levels)
  the_f <- stats::as.formula(paste0("~ 0 + ", setup$x), env = new.env())

  result <- list()
  if (isTRUE(include_global)) {
    result <- c(result, mgp(
      x = setup$x,
      a = setup$a,
      knots = setup$knots,
      ref_value = setup$x0,
      boundary_sd = boundary_sd,
      n_quad = setup$n_quad,
      init = init,
      lower = lower,
      upper = upper,
      parscale = parscale,
      log = log
    ))
  }

  hmgp_name <- paste(c(setup$x, "hmgp"), collapse = "_")
  slope_name <- paste(c(setup$x, "hslope_mgp"), collapse = "_")
  result[[hmgp_name]] <- methods::new(
    "hmgp",
    name = setup$x,
    label = hmgp_name,
    formula = the_f,
    p.order = 2L,
    ref_value = setup$x0,
    knots = setup$knots,
    init = as.numeric(init)[1L],
    lower = as.numeric(lower)[1L],
    upper = as.numeric(upper)[1L],
    parscale = as.numeric(parscale)[1L],
    log = as.logical(log)[1L],
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_mult",
    ad_kind = "random",
    package = "adlaplace",
    a = setup$a,
    shift = setup$cs,
    n_quad = setup$n_quad,
    fem = setup$fem,
    two_sided = isTRUE(setup$fem$two_sided),
    by = by_obj
  )
  result[[slope_name]] <- methods::new(
    "hmgp_boundary",
    name = setup$x,
    label = slope_name,
    formula = the_f,
    p.order = 1L,
    ref_value = setup$x0,
    knots = setup$knots,
    init = as.numeric(slope_init)[1L],
    lower = as.numeric(slope_lower)[1L],
    upper = as.numeric(slope_upper)[1L],
    parscale = as.numeric(slope_parscale)[1L],
    log = as.logical(slope_log)[1L],
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random",
    package = "adlaplace",
    sd = NA_real_,
    a = setup$a,
    shift = setup$cs,
    by = by_obj
  )
  result
}

#' @noRd
.hmgp_check_scalar <- function(x, name, positive = TRUE, allow_inf = FALSE) {
  ok <- length(x) == 1L && !is.na(x) &&
    (is.finite(x) || (allow_inf && is.infinite(x)))
  if (positive) {
    ok <- ok && x > 0
  }
  if (!ok) {
    stop(name, " must be a single ", if (positive) "positive " else "", "value", call. = FALSE)
  }
  invisible(x)
}

#' @noRd
.hmgp_k_labels <- function(term, n) {
  pad <- formatC(
    seq_len(n),
    width = max(1L, ceiling(log10(n + 1e-6))),
    flag = "0"
  )
  paste0(term@label, "_k", pad)
}

#' @describeIn hmgp-class Block design, one FEM basis per group.
#' @param term An `hmgp` term.
#' @param data Model data frame with the covariate and the grouping column.
#' @export
setMethod("design", "hmgp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  x <- as.numeric(data[[term@name]])
  if (any(!is.finite(x))) {
    stop("hmgp covariate '", term@name, "' has non-finite values", call. = FALSE)
  }
  tol <- sqrt(.Machine$double.eps)
  kn <- range(term@knots)
  if (any(x < kn[[1L]] - tol | x > kn[[2L]] + tol)) {
    stop(
      "hmgp covariate '", term@name, "' has values outside knots [",
      kn[[1L]], ", ", kn[[2L]], "]",
      call. = FALSE
    )
  }
  B <- mgp_fem_design(x, term@fem)
  g <- adlaplace::group_index(data, term@by@term[1L], term@by@levels)
  A <- adlaplace::scatter_by_group(B, g, length(term@by@levels))
  colnames(A) <- adlaplace::group_gamma_labels(
    .hmgp_k_labels(term, ncol(B)),
    term@by@labels
  )
  A
})

#' @describeIn hmgp-class Block-diagonal `random_mult` precision.
#' @export
setMethod("precision", "hmgp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  n_levels <- length(term@by@levels)
  if (n_levels < 1L) {
    stop("hmgp needs at least one group level", call. = FALSE)
  }
  Q <- Matrix::kronecker(Matrix::Diagonal(n_levels), term@fem$Q)
  Q <- Matrix::forceSymmetric(as_dgc_matrix(Q), uplo = "U")
  list(
    Q = Q,
    log_det = n_levels * term@fem$log_det,
    rank = n_levels * as.integer(term@fem$rank)
  )
})

#' @describeIn hmgp-class One shared SD theta.
#' @export
setMethod("theta_info", "hmgp", function(term) {
  data.frame(
    term = term@name,
    model = "hmgp",
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

#' @describeIn hmgp-class One row per FEM weight and group.
#' @export
setMethod("random_info", "hmgp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  n <- as.integer(term@fem$rank)
  result <- expand.grid(
    basis = seq_len(n),
    by = term@by@levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result$term <- term@name
  result$model <- "hmgp"
  result$label <- term@label
  result$order <- term@p.order
  result$by_labels <- term@by@labels[match(result$by, term@by@levels)]
  result$gamma_label <- adlaplace::group_gamma_labels(
    .hmgp_k_labels(term, n),
    term@by@labels
  )
  result
})

#' @describeIn hmgp-class No fixed effects.
#' @export
setMethod("beta_info", "hmgp", function(term, data) {
  NULL
})

#' @describeIn hmgp-class No companion log-det density.
#' @export
setMethod("extra_density", "hmgp", function(term) {
  NULL
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "hmgp", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  x <- methods::callNextMethod(x, map)
  if (length(x@by@term) && x@by@term[1L] %in% names(map)) {
    x@by@term <- as.character(map[[x@by@term[1L]]])
  }
  x@label <- paste(c(x@name, "hmgp"), collapse = "_")
  x
})

#' @describeIn hmgp_boundary-class Box-Cox column scattered by group.
#' @param term An `hmgp_boundary` term.
#' @param data Model data frame.
#' @export
setMethod("design", "hmgp_boundary", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  x <- as.numeric(data[[term@name]])
  u <- x - term@ref_value
  m <- mgp_boxcox(u, a = term@a, cs = term@shift)
  B <- matrix(m, ncol = 1L)
  g <- adlaplace::group_index(data, term@by@term[1L], term@by@levels)
  A <- adlaplace::scatter_by_group(B, g, length(term@by@levels))
  colnames(A) <- adlaplace::group_gamma_labels(term@label, term@by@labels)
  A
})

#' @describeIn hmgp_boundary-class Identity precision; the SD is a theta.
#' @export
setMethod("precision", "hmgp_boundary", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  n_levels <- length(term@by@levels)
  if (n_levels < 1L) {
    stop("hmgp_boundary needs at least one group level", call. = FALSE)
  }
  labels <- adlaplace::group_gamma_labels(term@label, term@by@labels)
  precision_mat <- Matrix::Diagonal(n = n_levels, x = 1)
  dimnames(precision_mat) <- list(labels, labels)
  precision_mat
})

#' @describeIn hmgp_boundary-class One SD theta.
#' @export
setMethod("theta_info", "hmgp_boundary", function(term) {
  data.frame(
    term = term@name,
    model = "hmgp_boundary",
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

#' @describeIn hmgp_boundary-class One slope coefficient per group.
#' @export
setMethod("random_info", "hmgp_boundary", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  data.frame(
    term = term@name,
    model = "hmgp_boundary",
    label = term@label,
    by = term@by@levels,
    by_labels = term@by@labels,
    basis = NA,
    order = 1L,
    gamma_label = adlaplace::group_gamma_labels(term@label, term@by@labels),
    stringsAsFactors = FALSE
  )
})

#' @describeIn hmgp_boundary-class No fixed effects.
#' @export
setMethod("beta_info", "hmgp_boundary", function(term, data) {
  NULL
})

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "hmgp_boundary", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  x <- methods::callNextMethod(x, map)
  if (length(x@by@term) && x@by@term[1L] %in% names(map)) {
    x@by@term <- as.character(map[[x@by@term[1L]]])
  }
  x@label <- paste(c(x@name, "hslope_mgp"), collapse = "_")
  x
})
