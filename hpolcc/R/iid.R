#' Independent and Identically Distributed Random Effects
#'
#' @description Creates an IID (independent and identically distributed) random effects term.
#'
#' @param x Variable name (typically a factor)
#' @param init Initial value for theta parameter
#' @param lower Lower bound for theta parameter
#' @param upper Upper bound for theta parameter
#' @param parscale Parameter scale for optimization
#' @param prefix Optional prefix for term names
#'
#' @return An iid term object

# IID class definition
setClass("iid",
         representation = representation(
         ),
         contains = "model",
         prototype = prototype(
          ref_value = numeric(0),
          p.order = as.integer(0),
          knots = numeric(0),
          by = character(0),
          type = factor("random", levels = .type_factor_levels)
         )
)

iid <- function(x,
                init = .my_theta_init,
                lower = .my_theta_lower,
                upper = .my_theta_upper,
                parscale = .my_theta_parscale) {
  new("iid",
    term = x,
    formula = formula(paste0("~ 0 + ", prefix, x)),
    init = init ,
    lower = lower ,
    upper = upper ,
    parscale = parscale
  )
}

# Design matrix for iid terms
setMethod("design", "iid", function(term, data){
  ff <- paste0("~ 0 + factor(", term@term, ")") |> formula()
  Matrix::sparse.model.matrix(ff, data) |> as("TsparseMatrix")
})

# Precision matrix for iid terms
setMethod("precision", "iid", function(term, data) {
  # Identity matrix for iid terms
  n <- length(unique(data[[term@term]]))
  Matrix::Diagonal(n, 1)
})

# Theta info for iid terms
setMethod("theta_info", "iid", function(term) {
  result <- data.frame(
    term = term@term, model = "iid", label = paste(c("iid", term@term), collapse = "_"),
    init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

# Beta info for iid terms
setMethod("beta_info", "iid", function(term) {
  # IID terms don't have beta parameters
  return(NULL)
})

# Gamma info for iid terms
setMethod("random_info", "iid", function(term, data) {
    

  result <- expand.grid(
    term = term@term,
    model = "iid",
    label = paste("iid", term@term, sep="_"),
    by = NA,
    by_label = NA,
    basis = sort(unique(data[[term@term]])),
    order = NA,
    stringsAsFactors = FALSE
  )
  result$gamma_label = paste(result$label, result$basis, sep="_")

  result
})
