#' Hierarchical Integrated Wiener Process Model Term
#'
#' @description Creates and manages hierarchical integrated Wiener process (HIWP) model terms.
#' @name hiwp-class
#' @aliases hiwp
#' @docType class
#' @exportClass hiwp
#' @title Hierarchical Integrated Wiener Process Model Term
#' @importFrom adlaplace iwp
#' @importFrom adlaplace rpoly
#' @importFrom adlaplace fpoly
#' @section Methods:
#' The following methods are available for `hiwp` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for HIWP term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for HIWP term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#'   \item{\code{replace_vars(x, map)}}{Rewrites data-column names, including the grouping factor}
#' }
NULL

setClass(
  "hiwp",
  slots = list(by = "by_group"),
  contains = "iwp",
  prototype = prototype(by = adlaplace::by_group())
)

#' @export
#' @rdname hiwp-class
#' @param x Variable name
#' @param p Order of the integrated Wiener process (default: 2)
#' @param ref_value Reference value for the basis
#' @param knots Vector of knot locations
#' @param by Grouping variable for hierarchical structure
#' @param levels Optional group levels, in block order. When omitted, levels
#'   are taken from the data on the first design or precision call.
#' @param init,lower,upper,parscale,log Scalar SD controls for the group
#'   deviations. Reused for the global \code{iwp} when \code{include_global}
#'   is \code{TRUE}.
#' @param slope_init,slope_lower,slope_upper,slope_parscale,slope_log SD
#'   controls for the per-group \code{hrpoly} slopes, recycled to \code{p - 1}.
#' @param boundary_is_random Whether the global boundary should be random
#' @param include_poly Whether to include polynomial terms
#' @param include_global Whether to include the global IWP and its boundary
#' @return A named list. With the defaults: \code{<x>_iwp}, \code{<x>_rpoly}
#'   (or \code{<x>_fpoly}), \code{<x>_hiwp}, then one \code{<x>_hrpoly_<d>}
#'   per slope order.
#' @examples
#' # Example usage:
#' # knots <- seq(0, 1, length.out = 5)
#' # hiwp_term <- hiwp(x = "age", knots = knots, by = "group")
hiwp <- function(
  x, p = 2, ref_value = 0, knots,
  by,
  levels = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  log = TRUE,
  slope_init = .my_theta_init,
  slope_lower = .my_theta_lower,
  slope_upper = .my_theta_upper,
  slope_parscale = .my_theta_parscale,
  slope_log = TRUE,
  boundary_is_random = TRUE,
  include_poly = TRUE,
  include_global = TRUE
) {
  if (missing(by) || length(by) != 1L || !nzchar(as.character(by))) {
    stop("by is required for hiwp()", call. = FALSE)
  }
  if (length(init) != 1L) stop("init must be a single value", call. = FALSE)
  if (length(lower) != 1L) stop("lower must be a single value", call. = FALSE)
  if (length(upper) != 1L) stop("upper must be a single value", call. = FALSE)
  if (length(parscale) != 1L) stop("parscale must be a single value", call. = FALSE)
  if (length(log) != 1L) stop("log must be a single value", call. = FALSE)
  if (length(knots) < 2L) stop("knots must have length >= 2", call. = FALSE)

  x <- adlaplace::strip_term_name(as.character(x))
  ref_value <- adlaplace::ref_align(ref_value, knots)
  by_obj <- adlaplace::by_group(term = as.character(by), levels = levels)
  the_f <- stats::as.formula(paste0("~ 0 + ", x), env = new.env())
  n_slope <- max(0L, as.integer(p) - 1L)
  if (n_slope > 0L) {
    slope_init <- rep_len(slope_init, n_slope)
    slope_lower <- rep_len(slope_lower, n_slope)
    slope_upper <- rep_len(slope_upper, n_slope)
    slope_parscale <- rep_len(slope_parscale, n_slope)
    slope_log <- rep_len(slope_log, n_slope)
  }

  result <- list()
  if (include_global) {
    result <- c(result, adlaplace::iwp(
      x = x,
      p = p,
      ref_value = ref_value,
      knots = knots,
      init = init,
      lower = lower,
      upper = upper,
      parscale = parscale,
      log = log,
      boundary_is_random = boundary_is_random,
      include_poly = include_poly
    ))
  }

  hiwp_name <- paste(c(x, "hiwp"), collapse = "_")
  result[[hiwp_name]] <- methods::new(
    "hiwp",
    name = x,
    label = hiwp_name,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    knots = knots,
    by = by_obj,
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    log = log,
    model_role = factor("random", levels = adlaplace::.model_role_levels),
    density = "random_diagonal",
    ad_kind = "random"
  )
  if (include_poly && n_slope > 0L) {
    for (D_poly in seq_len(n_slope)) {
      hrpoly_name <- paste(c(x, "hrpoly", D_poly), collapse = "_")
      result[[hrpoly_name]] <- hrpoly(
        x = x, p = D_poly, ref_value = ref_value,
        by = by_obj,
        init = slope_init[D_poly],
        lower = slope_lower[D_poly],
        upper = slope_upper[D_poly],
        parscale = slope_parscale[D_poly],
        log = slope_log[D_poly]
      )
    }
  }
  result
}

#' @noRd
.hiwp_k_labels <- function(term, n) {
  width <- max(1L, ceiling(log10(n + 1e-6)))
  pad <- formatC(seq_len(n), width = width, flag = "0")
  paste0(term@label, "_k", pad)
}

#' @rdname hiwp-class
#' @param term An `hiwp` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A design matrix for the HIWP term.
#' @export
setMethod("design", "hiwp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  A0 <- methods::callNextMethod()
  g <- adlaplace::group_index(data, term@by@term[1L], term@by@levels)
  A <- adlaplace::scatter_by_group(A0, g, length(term@by@levels))
  colnames(A) <- adlaplace::group_gamma_labels(
    .hiwp_k_labels(term, ncol(A0)),
    term@by@labels
  )
  A
})

#' @rdname hiwp-class
#' @param term An `hiwp` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A precision matrix for the HIWP term.
#' @export
setMethod("precision", "hiwp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  iwp_precision <- methods::callNextMethod()
  n_levels <- length(term@by@levels)
  result <- Matrix::bdiag(rep(list(iwp_precision), n_levels))
  labels <- adlaplace::group_gamma_labels(
    .hiwp_k_labels(term, ncol(iwp_precision)),
    term@by@labels
  )
  dimnames(result) <- list(labels, labels)
  result
})

#' @rdname hiwp-class
#' @param term An `hiwp` term object.
#' @return A data frame containing theta parameter information for the HIWP term.
#' @export
setMethod("theta_info", "hiwp", function(term) {
  # Global and local parameters
  result <- data.frame(
    term = term@name,
    model = "hiwp",
    label = paste(c(term@name, "hiwp"), collapse = "_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    model_role = term@model_role,
    log = term@log
  )

  result
})

#' @rdname hiwp-class
#' @param term An `hiwp` term object.
#' @param data A data frame containing the variables used in the term.
#' @return NULL (HIWP terms don't have beta parameters).
#' @export
setMethod("beta_info", "hiwp", function(term, data) {
  NULL
})

#' @rdname hiwp-class
#' @param term An `hiwp` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A data frame containing random effects information for the HIWP term.
#' @export
setMethod("random_info", "hiwp", function(term, data) {
  term <- adlaplace::add_by_levels(term, data)
  basis <- seq(1, length.out = length(term@knots) - 1)

  result <- expand.grid(
    basis = basis,
    by = term@by@levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result$term <- term@name
  result$model <- "hiwp"
  result$label <- term@label
  result$order <- term@p.order
  result$by_labels <- term@by@labels[match(result$by, term@by@levels)]
  result$gamma_label <- adlaplace::group_gamma_labels(
    .hiwp_k_labels(term, length(basis)),
    term@by@labels
  )
  result
})

