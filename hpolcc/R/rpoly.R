#' Random Polynomial Model Term
#'
#' @description Creates a random polynomial model term.
#'
#' @param x Variable name
#' @param p Polynomial degree (default: 2)
#' @param ref_value Reference value for the polynomial
#' @param sd Standard deviation for random effects
#' @param prefix Optional prefix for term names
#'
#' @return A rpoly term object

# RPoly class definition
setClass("rpoly",
         representation = representation(
           sd = "numeric"
         ),
         contains = "model",
         prototype = list(
           sd = numeric(0)
         )
)

rpoly <- function(x, p = 2, ref_value, sd = Inf, prefix = NULL) {
  new("rpoly",
    term = x,
    formula = formula(paste0("~ 0 + ", prefix, x)),
    p.order = p,
    ref_value = ref_value,
    sd = rep_len(sd, p),
    type = factor("random", levels = .type_factor_levels)
  )
}

# Design matrix for rpoly terms
setMethod("design", "rpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }

  D <- poly(data[[term@term]]-term@ref_value, raw=TRUE, degree = term@p.order)
  D <- D[,1:ncol(D),drop=FALSE]
  colnames(D) <- paste(term@term, 1:ncol(D), sep="_")
  D
})

# Precision matrix for rpoly terms
setMethod("precision", "rpoly", function(term, data) {
  Matrix::Diagonal(term@p.order, 1)
})

# Theta info for rpoly terms
setMethod("theta_info", "rpoly", function(term) {
  NULL
})

# Beta info for rpoly terms
setMethod("beta_info", "rpoly", function(term) {
  # Rpoly terms don't have beta parameters (random effects only)
  return(NULL)
})

# Gamma info for rpoly terms
setMethod("random_info", "rpoly", function(term, data) {
  order <- seq_len(term@p)
  basis <- NA

  result <- expand.grid(
    term = term@term,
    model = "rpoly",
    name = term@term,
    groups = NA,
    groups_string = NA,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE,
    type = as.character(term@type)
  )
  result$gamma_label <- paste(result$term, result$order, sep = "_")

  result
})
