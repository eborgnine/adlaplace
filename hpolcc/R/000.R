#' @useDynLib hpolcc, .registration = TRUE
#' @importFrom Rcpp evalCpp

# Define factor levels for model type
.type_factor_levels <- c("fixed", "random", "family")

# Define generic functions for model methods
setGeneric("design", function(term, data) standardGeneric("design"))
setGeneric("precision", function(term, data) standardGeneric("precision"))
setGeneric("theta_info", function(term) standardGeneric("theta_info"))
setGeneric("beta_info", function(term, data) standardGeneric("beta_info"))
setGeneric("random_info", function(term, data) standardGeneric("random_info"))

# Define the base model class that all specific model classes inherit from
setClass("model",
         representation(
           term = "character",
           formula = "formula",
           knots = "numeric",
           ref_value = "numeric",
           p.order = "numeric",
           by = "character",
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric",
           type = "factor"
         ),
         prototype = list(
           term = character(0),
           formula = formula(),
           knots = numeric(0),
           ref_value = as.numeric(NA),
           p.order = 1,
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

ref_align <- function(ref_value, knots) {
  knots[which.min(abs(knots - ref_value))]
}
