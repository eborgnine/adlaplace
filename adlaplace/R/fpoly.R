#' Fixed Polynomial Model Term
#'
#' @description Creates and manages fixed polynomial model terms.
#' @name fpoly-class
#' @aliases fpoly
#' @docType class
#' @exportClass fpoly
#'
#' @section Methods:
#' The following methods are available for `fpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for fpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for fpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("fpoly",
         slots = list(
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         contains = "model",
         prototype = prototype(
           knots = numeric(0),
           by = character(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           type = factor("fixed", levels = .type_factor_levels)
         )
)

#' @rdname fpoly-class
#' @param x Variable name.
#' @param x Variable name.
#' @param p Polynomial degree (default: 2).
#' @param ref_value Reference value for the polynomial.
#' @param init Initial values for beta parameters.
#' @param lower Lower bounds for beta parameters.
#' @param upper Upper bounds for beta parameters.
#' @param parscale Parameter scales for optimization.
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A `fpoly` term object.
#' @export
#' @examples
#' # Create a fixed polynomial term
#' fpoly_term <- fpoly(x = "temp", p = 3, ref_value = 15)
fpoly <- function(x, p = 2, ref_value = 0,
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale) {
  methods::new("fpoly",
    term = x,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    init = rep_len(init, p),
    lower = rep_len(lower, p),
    upper = rep_len(upper, p),
    parscale = rep_len(parscale, p)
  )
}

#' @describeIn fpoly-class Design method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A design matrix for the fixed polynomial term, or NULL if p.order is 0.
setMethod("design", "fpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  D <- stats::poly(
    data[[term@term]] - term@ref_value,
    degree = term@p.order,
    raw = TRUE
  )
  D <- D[, 1:ncol(D), drop = FALSE]
  seq_order <- seq.int(1, length.out = term@p.order)

  colnames(D) <- paste0(term@term, "_fpoly_", seq_order)
  D
})

#' @describeIn fpoly-class Precision method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return NULL (fixed effects don't have precision matrices).
setMethod("precision", "fpoly", function(term, data) {
  # Fixed effects don't have precision matrices
  NULL
})

#' @describeIn fpoly-class Theta info method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @return NULL (fixed effects don't have random effects parameters).
setMethod("theta_info", "fpoly", function(term) {
  # Fixed effects don't have random effects parameters
  return(NULL)
})

#' @describeIn fpoly-class Beta info method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A data frame containing beta parameter information for the fixed polynomial term.
setMethod("beta_info", "fpoly", function(term, data) {
  the_colnames <- colnames(design(term, data))
  the_label <- paste(term@term, "fpoly", sep = "_")

  result <- data.frame(
    term = term@term,
    model = "fpoly",
    label = the_label,
    order = NA,
    beta_label = the_colnames,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )

  return(result)
})

#' @describeIn fpoly-class Random info method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return NULL (fixed effects don't have random effects information).
setMethod("random_info", "fpoly", function(term, data) {
  # Fixed effects don't have random effects information
  return(NULL)
})
