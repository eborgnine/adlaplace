#' Overdispersion Class
#'
#' @description
#' A class for representing overdispersion terms in hierarchical models.
#' This class handles the additional variance introduced by overdispersion
#' in count data models.

# Overdispersion class definition
setClass("overdispersion",
  representation = representation(
  ),
  contains = "model",           # Inherits from base model class
  prototype = prototype(
  )
)

#' @title Create an Overdispersion Term
#'
#' @description
#' Creates an overdispersion term for hierarchical models.
#'
#' @param init Initial value for the overdispersion parameter (default: 1e-3)
#' @param lower Lower bound for the overdispersion parameter (default: 0)
#' @param upper Upper bound for the overdispersion parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#'
#' @return An overdispersion object
#'
#' @examples
#' overdisp_term <- overdispersion()
#' print(overdisp_term)

overdispersion <- function(
  x = NA,
  init = 1e-3,
  lower = 0,
  upper = Inf,
  parscale = 1
) {
  new("overdispersion",
    term = character(0),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    type = factor("family", levels = .type_factor_levels)
  )
}

# Method implementations for overdispersion class

#' @title Design Matrix for Overdispersion
#'
#' @description
#' Creates design matrix for overdispersion term.
#'
#' @param object overdispersion object
#' @param data data frame containing the variables
#' @param ... additional arguments
#'
#' @return NULL (overdispersion doesn't contribute to design matrix)

setMethod("design", "overdispersion", function(term, data) {
  NULL
})

#' @title Precision Matrix for Overdispersion
#'
#' @description
#' Creates precision matrix for overdispersion term.
#'
#' @param object overdispersion object
#' @param ... additional arguments
#'
#' @return NULL (overdispersion doesn't have random effects)

setMethod("precision", "overdispersion", function(term) {
  NULL
})

#' @title Theta Information for Overdispersion
#'
#' @description
#' Extracts theta parameter information for overdispersion term.
#'
#' @param object overdispersion object
#' @param ... additional arguments
#'
#' @return data frame with theta parameter information

setMethod("theta_info", "overdispersion", function(term) {
  data.frame(
    term = NA,
    model = "overdispersion",
    order = NA,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    stringsAsFactors = FALSE
  )
})

#' @title Beta Information for Overdispersion
#'
#' @description
#' Extracts beta parameter information for overdispersion term.
#'
#' @param object overdispersion object
#' @param data data frame
#' @param ... additional arguments
#'
#' @return NULL (overdispersion doesn't have beta parameters)

setMethod("beta_info", "overdispersion", function(term, data) {
  NULL
})

#' @title Random Information for Overdispersion
#'
#' @description
#' Extracts random effects information for overdispersion term.
#'
#' @param object overdispersion object
#' @param data data frame
#' @param ... additional arguments
#'
#' @return NULL (overdispersion doesn't have random effects)

setMethod("random_info", "overdispersion", function(term, data) {
  NULL
})