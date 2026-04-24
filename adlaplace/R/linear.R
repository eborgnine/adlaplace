#' Linear Model Term
#'
#' @description Creates and manages linear model terms for fixed effects.
#' @name linear-class
#' @aliases linear
#' @docType class
#' @title Linear Model Term
#' @exportClass linear
#'
#' @section Methods:
#' The following methods are available for `linear` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for linear term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for linear term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("linear",
         representation = representation(
         ),
         contains = "model",
         prototype = list(
                    ref_value = numeric(0),
          p.order = integer(0),
          knots = numeric(0),
          by = character(0),
    type = factor("fixed", levels = .type_factor_levels)
         )
)

#' @export
#' @rdname linear-class
#' @param x Variable name
#' @param init Initial value for beta parameter (default: 0)
#' @param lower Lower bound for beta parameter (default: -Inf)
#' @param upper Upper bound for beta parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#' @return A linear term object
#' @examples
#' # Create a linear term
#' linear_term <- linear(x = "temperature")
#' @export
linear <- function(x, 
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale) {
  methods::new("linear",
    term = x,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn linear-class Creates design matrix for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "linear", function(term, data){
    res <- Matrix::sparse.model.matrix(term@formula, data, drop.unused.levels = FALSE)
    if (is.factor(data[[term@term]])) {
      res <- res[, -1, drop = FALSE]
    }
    colnames(res) = paste0(term@term, "_linear_", colnames(res))
    res
})

#' @describeIn linear-class Creates precision matrix for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "linear", function(term, data) {
  # Linear terms don't have precision matrices
  NULL
})

#' @describeIn linear-class Extracts theta parameter information for linear term
#' @param term A linear term object
#' @export
setMethod("theta_info", "linear", function(term) {
  # Linear terms don't have random effects parameters
  return(NULL)
})

#' @describeIn linear-class Extracts beta parameter information for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "linear", function(term, data) {
  the_colnames = colnames(design(term, data))
  the_label = paste(term@term, "linear", sep="_")

  result <- data.frame(
    term = term@term,
    model = "linear",
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

#' @describeIn linear-class Extracts random effects information for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "linear", function(term, data) {
  NULL
})
