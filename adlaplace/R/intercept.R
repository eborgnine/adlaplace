#' intercept Model Term
#'
#' @description Creates and manages intercept model terms for fixed effects.
#' @name intercept-class
#' @aliases intercept
#' @docType class
#' @title intercept Model Term
#' @exportClass intercept
#'
#' @section Methods:
#' The following methods are available for `intercept` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for intercept term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for intercept term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL


setClass("intercept",
  representation = representation(),
  contains = "model",
  prototype = list(
    term = "intercept",
    formula = . ~ 1,
    ref_value = numeric(0),
    p.order = integer(0),
    knots = numeric(0),
    by = character(0),
    type = factor("fixed", levels = .type_factor_levels)
  )
)

#' @export
#' @rdname intercept-class
#' @param init Initial value for beta parameter (default: 0)
#' @param lower Lower bound for beta parameter (default: -Inf)
#' @param upper Upper bound for beta parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#' @return A intercept term object
#' @examples
#' # Create a intercept term
#' intercept_term <- intercept(init=2)
#' @export
intercept <- function(init = .my_beta_init,
                      lower = .my_beta_lower,
                      upper = .my_beta_upper,
                      parscale = .my_beta_parscale) {
  methods::new("intercept",
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn intercept-class Creates design matrix for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "intercept", function(term, data) {
  matrix(1, nrow(data), ncol=1,
    dimnames = list(rownames(data),"intercept")
    )
})

#' @describeIn intercept-class Creates precision matrix for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "intercept", function(term, data) {
  # intercept terms don't have precision matrices
  NULL
})

#' @describeIn intercept-class Extracts theta parameter information for intercept term
#' @param term A intercept term object
#' @export
setMethod("theta_info", "intercept", function(term) {
  # intercept terms don't have random effects parameters
  return(NULL)
})

#' @describeIn intercept-class Extracts beta parameter information for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "intercept", function(term, data) {
  
  data.frame(
    term = "intercept",
    model = "intercept",
    label = "intercept",
    order = NA,
    beta_label = "intercept",
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )
})

#' @describeIn intercept-class Extracts random effects information for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "intercept", function(term, data) {
  NULL
})
