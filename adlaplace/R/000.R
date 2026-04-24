#' @useDynLib adlaplace, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL



# Define factor levels for model type
#' @export
.type_factor_levels <- c("fixed", "random", "response")

#' Model Term Generics
#'
#' @description Generic functions for extracting model matrices and parameter information from model terms.
#' @param term A model term object.
#' @param data A data frame containing the variables.
#' @return Method-specific return values.
#' @name model-generics
NULL

#' @rdname model-generics
#' @export
setGeneric("design", function(term, data) standardGeneric("design"))
#' @rdname model-generics
#' @export
setGeneric("precision", function(term, data) standardGeneric("precision"))
#' @rdname model-generics
#' @export
setGeneric("theta_info", function(term) standardGeneric("theta_info"))
#' @rdname model-generics
#' @export
setGeneric("beta_info", function(term, data) standardGeneric("beta_info"))
#' @rdname model-generics
#' @export
setGeneric("random_info", function(term, data) standardGeneric("random_info"))

#' Base Model Class
#'
#' @description The base S4 class that all specific model term classes inherit from.
#' @exportClass model
setClass("model",
         slots = list(
           term = "character",
           formula = "formula",
           knots = "numeric",
           ref_value = "numeric",
           p.order = "integer",
           by = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           type = "factor"
         ),
         prototype = prototype(
           term = character(0),
           formula = formula(),
           knots = numeric(0),
           ref_value = numeric(0),
           p.order = as.integer(1),
           by = character(0),
           init = numeric(0),
           lower = numeric(0),  
           upper = numeric(0),
           parscale = numeric(0),
           type = factor("fixed", levels = .type_factor_levels)
         )
)

.my_beta_init <- 0
.my_beta_lower <- -Inf
.my_beta_upper <- Inf
.my_beta_parscale <- 1

.my_theta_init <- 0.02
.my_theta_lower <- 1e-9
.my_theta_upper <- Inf
.my_theta_parscale <- 1
