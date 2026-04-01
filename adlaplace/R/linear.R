#' Linear Model Term
#'
#' @description Creates a linear model term for fixed effects.
#'
#' @param x Variable name
#' @param prefix Optional prefix for term names
#' @param init Initial value for beta parameter (default: 0)
#' @param lower Lower bound for beta parameter (default: -Inf)
#' @param upper Upper bound for beta parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#'
#' @return A linear term object
#'
#' @examples
#' # Create a linear term
#' linear_term <- linear(x = "temperature")

# Linear class definition
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
linear <- function(x, 
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale) {
  new("linear",
    term = x,
    formula = as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

# Design matrix for linear terms
setMethod("design", "linear", function(term, data){
    res <- Matrix::sparse.model.matrix(term@formula, data, drop.unused.levels = FALSE)
    if (is.factor(data[[term@term]])) {
      res <- res[, -1, drop = FALSE]
    }
    colnames(res) = paste0(term@term, "_linear_", colnames(res))
    res
})

# Precision matrix for linear terms
setMethod("precision", "linear", function(term, data) {
  # Linear terms don't have precision matrices
  NULL
})

# Theta info for linear terms
setMethod("theta_info", "linear", function(term) {
  # Linear terms don't have random effects parameters
  return(NULL)
})

# Beta info for linear terms
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

# Gamma info for linear terms
setMethod("random_info", "linear", function(term, data) {
  NULL
})
