#' Hierarchical Random Polynomial Model Term
#'
#' @description Creates a hierarchical random polynomial model term.
#'
#' @param x Variable name
#' @param p Polynomial degree (default: 1)
#' @param ref_value Reference value for the polynomial
#' @param by Grouping variable for hierarchical structure
#' @param init Initial value for theta parameter
#' @param lower Lower bound for theta parameter
#' @param upper Upper bound for theta parameter
#' @param parscale Parameter scale for optimization
#' @param prefix Optional prefix for term names
#'
#' @return A hrpoly term object

# HRPoly class definition
setClass("hrpoly",
         representation = representation(
           by_levels = "integer",
           by_labels = "character"
         ),
         contains = "model",
         prototype = list(
           by_levels = integer(0),
           by_labels = character(0)
         )
)

hrpoly <- function(
  x,
  p = 1,
  ref_value,
  by,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  prefix = NULL
) {
  new("hrpoly",
    term = x,
    formula = formula(paste0("~ 0 + ", prefix, x)),
    p.order = p,
    ref_value = ref_value,
    by = by,
    by_levels = integer(0),  # Will be set later when data is available
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1],
    type = factor("random", levels = .type_factor_levels)
  )
}

# Design matrix for hrpoly terms
setMethod("design", "hrpoly", function(term, data) {

  id_split <- split(1:nrow(data), 
    factor(data[[term@by]], levels = term@by_levels), 
  drop = FALSE)
  mm <- c(0, sapply(id_split, length)) |> cumsum()

  A0 <- drop(poly(
    data[[term@term]] - term@ref_value, 
    raw = TRUE, degree = term@p.order
  )[,term@p.order])
    

  Afinal <- Matrix::Matrix(0, nrow = nrow(data), ncol = length(term@groups)) |> as("TsparseMatrix")

  for (k in seq_along(id_split)) {
    Afinal[id_split[[k]], k] <- A0[id_split[[k]]]
  }
  colnames(Afinal) <- paste(
    term@term,
    term@by_labels,
    sep = "_"
  )

  Afinal
})

# Precision matrix for hrpoly terms
setMethod("precision", "hrpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  unique_values <- term@by_levels
  if (!length(unique_values)) {
    unique_values <- unique(data[[term@term]])
  }
  Matrix::Diagonal(length(unique_values), 1)
})

# Theta info for hrpoly terms
setMethod("theta_info", "hrpoly", function(term) {
  result <- data.frame(
    term = term@term, model = "hrpoly", 
    label = paste(c("hrpoly", term@term), collapse = "_"),
    global = FALSE, 
    order = term@p.order,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )
  return(result)
})

# Beta info for hrpoly terms
setMethod("beta_info", "hrpoly", function(term) {
  # HRpoly terms don't have beta parameters (random effects only)
  return(NULL)
})

# Gamma info for hrpoly terms
setMethod("random_info", "hrpoly", function(term, data) {
  basis <- NA
  if(!length(term@by_levels)) {
    term@by_levels <- unique(data[[term@by]])
  }

  result <- expand.grid(
    term = term@term,
    model = "hrpoly",
    name = term@term,
    by = term@by_levels,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result$by_labels <- term@by_labels[match(result$by, term@by_levels)]
  result$gamma_label <- paste(result$term, result$by_labels, sep = "_")
  result
})
