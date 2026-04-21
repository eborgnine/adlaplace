#' Overdispersion Class
#'
#' @description A class for representing overdispersion terms in hierarchical models.
#' This class handles the additional variance introduced by overdispersion
#' in count data models.
#' @name overdispersion-class
#' @aliases overdispersion
#' @docType class
#' @title Overdispersion Class
#' @exportClass overdispersion
#'
#' @section Methods:
#' The following methods are available for `overdispersion` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for overdispersion term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for overdispersion term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("overdispersion",
  representation = representation(
  ),
  contains = "model",           # Inherits from base model class
  prototype = prototype(
           term = character(0),
           formula = formula(),
           knots = numeric(0),
           ref_value = numeric(0),
           p.order = integer(0),
           by = character(0)
  )
)

#' Overdispersion Term Constructor
#'
#' @description Creates an overdispersion term for hierarchical models.
#' @rdname overdispersion-class
#' @param x Variable name (not used, included for interface consistency)
#' @param init Initial value for the overdispersion parameter (default: 1e-3)
#' @param lower Lower bound for the overdispersion parameter (default: 0)
#' @param upper Upper bound for the overdispersion parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#' @return An overdispersion object
#' @examples
#' overdisp_term <- overdispersion()
#' print(overdisp_term)
#' @export
overdispersion <- function(
  x = NA,
  init = 1e-3,
  lower = 0,
  upper = Inf,
  parscale = 1
) {
  methods::new("overdispersion",
    term = character(0),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    type = factor("family", levels = .type_factor_levels)
  )
}

# Method implementations for overdispersion class
#' @describeIn overdispersion-class Creates design matrix for overdispersion term
#' @param term An overdispersion term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "overdispersion", function(term, data) {
  NULL
})

#' @describeIn overdispersion-class Creates precision matrix for overdispersion term
#' @param term An overdispersion term object
#' @export
setMethod("precision", "overdispersion", function(term) {
  NULL
})

#' @describeIn overdispersion-class Extracts theta parameter information for overdispersion term
#' @param term An overdispersion term object
#' @export
setMethod("theta_info", "overdispersion", function(term) {
  data.frame(
    term = NA,
    model = "overdispersion",
    label = "overdispersion",
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    stringsAsFactors = FALSE
  )
})

#' @describeIn overdispersion-class Extracts beta parameter information for overdispersion term
#' @param term An overdispersion term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "overdispersion", function(term, data) {
  NULL
})

#' @describeIn overdispersion-class Extracts random effects information for overdispersion term
#' @param term An overdispersion term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "overdispersion", function(term, data) {
  NULL
})