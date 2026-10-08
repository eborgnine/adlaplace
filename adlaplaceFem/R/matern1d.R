#' One-dimensional linear B-spline FEM (mass and stiffness)
#'
#' @param knots Strictly increasing knot positions, including the domain ends.
#' @param degree B-spline degree. The Matérn-1/2 term uses degree 1.
#' @return List with open `knots`, `degree`, sparse Grams `C` and `G`, and
#'   `n` basis functions.
#' @export
fem_bspline_1d <- function(knots, degree = 1L) {
  degree <- as.integer(degree)
  knots <- sort(unique(as.numeric(knots)))
  if (length(knots) < 2L || any(!is.finite(knots))) {
    stop("knots must have at least two finite positions", call. = FALSE)
  }
  open <- axis_to_open_knots(knots, degree)
  C <- gram_1d(open, degree, 0L, 0L)
  G <- gram_1d(open, degree, 1L, 1L)
  list(knots = open, degree = degree, C = C, G = G, n = nrow(C))
}

#' @noRd
.rsmatern_levels <- function(term, data) {
  if (length(term@by@levels)) {
    return(term@by@levels)
  }
  by_var <- term@by@term[1L]
  if (!length(by_var) || !by_var %in% names(data)) {
    stop("rsmatern grouping variable not found in data", call. = FALSE)
  }
  x <- data[[by_var]]
  if (is.factor(x)) levels(x) else as.character(unique(x))
}

#' @noRd
.rsmatern_gamma_labels <- function(term, levels, n_basis) {
  bpad <- formatC(
    seq_len(n_basis),
    width = max(1L, ceiling(log10(n_basis + 1e-6))),
    flag = "0"
  )
  unlist(lapply(levels, function(lv) {
    paste0(term@label, "_", lv, "_b", bpad)
  }), use.names = FALSE)
}

#' @noRd
.rsmatern_fem <- function(term) {
  if (is.list(term@fem) && !is.null(term@fem$C)) {
    return(term@fem)
  }
  fem_bspline_1d(term@knots, degree = 1L)
}

#' Random-slope Matérn-1/2 (Ornstein-Uhlenbeck) FEM term
#'
#' Stationary one-dimensional Matérn field with smoothness 1/2, the continuous
#' analogue of AR(1), multiplied by `(mult - ref_mult)`. Groups in `by` share
#' one `(range, sd)` and have independent fields. The precision is the linear
#' B-spline FEM
#' `Q = tau^2 (kappa^2 C + G)` with practical range `rho = 2/kappa` and
#' `tau = 1 / (sd * sqrt(2 * kappa))`.
#'
#' @slot mult Exposure column multiplied into the basis.
#' @slot ref_mult Reference subtracted from `mult`.
#' @slot by Grouping variable; one independent field per level.
#' @slot fem Cached 1D Grams from [fem_bspline_1d()].
#' @export
#' @rdname rsmatern-class
setClass(
  "rsmatern",
  slots = c(
    mult = "character",
    ref_mult = "numeric",
    by = "by_group",
    fem = "ANY"
  ),
  contains = "model_term",
  prototype = prototype(
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_fem_ssq_1",
    ad_kind = "random",
    package = "adlaplaceFem",
    p.order = 1L
  )
)

#' @param x Time column.
#' @param mult Exposure column. The basis row is multiplied by
#'   `data[[mult]] - ref_mult`.
#' @param by Grouping column. Levels share `(range, sd)`.
#' @param knots Increasing knot positions on the time axis. Linear B-splines.
#' @param ref_mult Reference exposure (default 0).
#' @param levels Optional group levels, in the order used for the blocks.
#' @param init Length-2 starting values `(range, sd)`.
#' @param lower,upper Length-2 bounds. The default lower range is the knot
#'   spacing.
#' @param log Log-transform the two hyperparameters (default `TRUE`).
#' @param parscale Optional length-2 step scale.
#' @param term An [`rsmatern-class`] term.
#' @param data Model data frame (passed by adlaplace generics).
#' @return An `rsmatern` term.
#' @export
#' @rdname rsmatern-class
rsmatern <- function(
    x,
    mult,
    by,
    knots,
    ref_mult = 0,
    levels = NULL,
    init = c(1, 0.05),
    lower = NULL,
    upper = NULL,
    log = TRUE,
    parscale = NULL) {
  knots <- sort(unique(as.numeric(knots)))
  if (length(knots) < 2L) {
    stop("knots must have at least two positions", call. = FALSE)
  }
  spacing <- min(diff(knots))
  span <- max(knots) - min(knots)
  if (is.null(lower)) {
    lower <- c(spacing, 1e-6)
  }
  if (is.null(upper)) {
    upper <- c(span, 2)
  }
  init <- as.numeric(init)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  if (length(init) != 2L || length(lower) != 2L || length(upper) != 2L) {
    stop("init, lower, and upper must each have length 2 (range, sd)", call. = FALSE)
  }
  if (is.null(parscale)) {
    parscale <- rep(0.01, 2L)
  }
  fem <- fem_bspline_1d(knots, degree = 1L)
  label <- paste(c(x, mult, "rsmatern"), collapse = "_")
  methods::new(
    "rsmatern",
    name = x,
    label = label,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_fem_ssq_1",
    ad_kind = "random",
    package = "adlaplaceFem",
    p.order = 1L,
    init = init,
    lower = lower,
    upper = upper,
    log = rep_len(as.logical(log), 2L),
    parscale = as.numeric(parscale),
    mult = as.character(mult),
    ref_mult = as.numeric(ref_mult),
    by = adlaplace::by_group(term = as.character(by), levels = levels),
    knots = knots,
    fem = fem
  )
}

#' @describeIn rsmatern-class Design: basis at `x`, times `(mult - ref_mult)`,
#'   written into the block for `by`.
#' @export
setMethod("design", "rsmatern", function(term, data) {
  fem <- .rsmatern_fem(term)
  if (!term@name %in% names(data)) {
    stop("column '", term@name, "' not found", call. = FALSE)
  }
  if (!term@mult %in% names(data)) {
    stop("column '", term@mult, "' not found", call. = FALSE)
  }
  levels <- .rsmatern_levels(term, data)
  by_var <- term@by@term[1L]
  reg <- match(as.character(data[[by_var]]), levels)
  if (anyNA(reg)) {
    stop("rsmatern: data has grouping levels that are not in the term", call. = FALSE)
  }
  B <- bspline_eval(fem$knots, data[[term@name]], fem$degree, 0L)
  B <- methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix")
  scale <- as.numeric(data[[term@mult]]) - term@ref_mult
  B <- methods::as(Matrix::Diagonal(x = scale) %*% B, "CsparseMatrix")
  nb <- ncol(B)
  counts <- diff(B@p)
  cols <- rep.int(seq_along(counts) - 1L, counts)
  rows <- B@i
  col_out <- cols + (reg[rows + 1L] - 1L) * nb
  A <- Matrix::sparseMatrix(
    i = rows + 1L,
    j = col_out + 1L,
    x = B@x,
    dims = c(nrow(B), nb * length(levels)),
    repr = "C"
  )
  colnames(A) <- .rsmatern_gamma_labels(term, levels, nb)
  A
})

#' @describeIn rsmatern-class Block-diagonal FEM precision, one block per group.
#' @export
setMethod("precision", "rsmatern", function(term, data) {
  fem <- .rsmatern_fem(term)
  levels <- .rsmatern_levels(term, data)
  nlev <- length(levels)
  if (nlev < 1L) {
    stop("rsmatern needs at least one group level", call. = FALSE)
  }
  eye <- Matrix::sparseMatrix(
    i = seq_len(nlev), j = seq_len(nlev), x = 1,
    dims = c(nlev, nlev)
  )
  fem_precision_payload(
    list(C = Matrix::kronecker(eye, fem$C), G = Matrix::kronecker(eye, fem$G)),
    alpha = 1L
  )
})

#' @describeIn rsmatern-class One metadata row per basis weight and group.
#' @export
setMethod("random_info", "rsmatern", function(term, data) {
  fem <- .rsmatern_fem(term)
  levels <- .rsmatern_levels(term, data)
  basis <- seq_len(fem$n)
  result <- expand.grid(
    basis = basis,
    by = levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  # Block order matches design(): all basis functions of level 1, then level 2.
  result <- result[order(match(result$by, levels), result$basis), , drop = FALSE]
  rownames(result) <- NULL
  result$term <- term@name
  result$model <- "rsmatern"
  result$label <- term@label
  result$by_labels <- result$by
  result$order <- NA
  result$gamma_label <- .rsmatern_gamma_labels(term, levels, fem$n)
  result
})

#' @describeIn rsmatern-class Hyperparameters `(range, sd)`.
#' @export
setMethod("theta_info", "rsmatern", function(term) {
  data.frame(
    term = term@name,
    model = "rsmatern",
    label = paste0(term@label, c("_range", "_sd")),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    model_role = term@model_role,
    log = rep_len(term@log, 2L),
    stringsAsFactors = FALSE
  )
})

#' @describeIn rsmatern-class No fixed effects.
#' @export
setMethod("beta_info", "rsmatern", function(term, data) {
  NULL
})

#' @describeIn rsmatern-class Log-determinant density `random_fem_det_1`.
#' @export
setMethod("extra_density", "rsmatern", function(term) {
  "random_fem_det_1"
})

#' Replace data-column names in an rsmatern term
#'
#' Method for [adlaplace::replace_vars()]. Rewrites the time column, the
#' exposure column (`mult`), and the grouping variable (`by`).
#'
#' @param x An `rsmatern` term.
#' @param map Named character vector: names are current column names, values
#'   are replacement column names. Names absent from `map` are left unchanged.
#' @return The term with column references and the label updated.
#' @seealso [adlaplace::replace_vars()]
#' @name replace_vars-methods
#' @rdname replace_vars-methods
NULL

#' @rdname replace_vars-methods
#' @export
setMethod("replace_vars", "rsmatern", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  x <- methods::callNextMethod(x, map)
  if (length(x@mult) && x@mult %in% names(map)) {
    x@mult <- as.character(map[[x@mult]])
  }
  if (length(x@by@term) && x@by@term[1L] %in% names(map)) {
    x@by@term <- as.character(map[[x@by@term[1L]]])
  }
  x@label <- paste(c(x@name, x@mult, "rsmatern"), collapse = "_")
  x
})
