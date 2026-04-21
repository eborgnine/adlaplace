#' Random Slope Fixed Polynomial Model Term
#'
#' @description Creates a fixed polynomial model term for random slope models.
#' @name rsfpoly-class
#' @docType class
#' @exportClass rsfpoly
#'
#' @section Methods:
#' The following methods are available for `rsfpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for rsfpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for rsfpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("rsfpoly",
  representation = representation(
    mult = "character",
    ref_mult = "numeric"
  ),
  contains = "model",
  prototype = list(
    by = character(0),
    knots = numeric(0),
    type = factor("fixed", levels = adlaplace::.type_factor_levels)
  )
)

#' @rdname rsfpoly-class
#' @export
rsfpoly <- function(
  x, mult, p = 2,
  ref_value = 0, ref_mult = 0,
  init = .my_beta_init,
  lower = .my_beta_lower,
  upper = .my_beta_upper,
  parscale = .my_beta_parscale
) {
  methods::new("rsfpoly",
    term = x,
    mult = mult,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    ref_mult = ref_mult,
    init = rep_len(init, p),
    lower = rep_len(lower, p),
    upper = rep_len(upper, p),
    parscale = rep_len(parscale, p)
  )
}

#' @describeIn rsfpoly-class Creates design matrix for rsfpoly term
#' @param term A rsfpoly term object
#' @param data A data frame containing the term variables
#' @return A design matrix for the random slope fixed polynomial term, or NULL if p.order is 0
#' @export
setMethod("design", "rsfpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }

  mult_vec <- data[[term@mult]] - term@ref_mult

  D <- stats::poly(
    data[[term@term]] - term@ref_value,
    degree = term@p.order,
    raw = TRUE
  )
  D <- D[, 1:ncol(D), drop = F]
  D <- D * mult_vec
  seq_order <- seq.int(1, length.out = term@p.order)

  colnames(D) <- paste(term@term, term@mult, "rsfpoly", seq_order, sep = "_")
  D
})

#' @describeIn rsfpoly-class Creates precision matrix for rsfpoly term
#' @param term A rsfpoly term object
#' @param data A data frame containing the term variables
#' @return NULL (fixed effects don't have precision matrices)
#' @export
setMethod("precision", "rsfpoly", function(term, data) {
  # Fixed effects don't have precision matrices
  NULL
})

#' @describeIn rsfpoly-class Extracts theta parameter information for rsfpoly term
#' @param term A rsfpoly term object
#' @return NULL (fixed effects don't have theta parameters)
#' @export
setMethod("theta_info", "rsfpoly", function(term) {
  # not thetas for fpoly
  NULL
})

#' @describeIn rsfpoly-class Extracts beta parameter information for rsfpoly term
#' @param term A rsfpoly term object
#' @param data A data frame containing the term variables
#' @return A data frame containing beta parameter information for the random slope fixed polynomial term
#' @export
setMethod("beta_info", "rsfpoly", function(term, data) {
  the_label <- paste(term@term, term@mult, "rsfpoly", sep = "_")
  seq_order <- seq.int(1, length.out = term@p.order)

  result <- data.frame(
    term = term@term,
    model = "fpoly",
    label = the_label,
    order = seq_order,
    beta_label = paste(the_label, seq_order, sep = "_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = as.character(term@type)
  )
  return(result)
})

#' @describeIn rsfpoly-class Extracts random effects information for rsfpoly term
#' @param term A rsfpoly term object
#' @param data A data frame containing the term variables
#' @return NULL (fixed effects don't have random effects information)
#' @export
setMethod("random_info", "rsfpoly", function(term, data) {
  # no gammas for fpoly
  NULL
})
