#' Integrated Wiener Process Model Term
#'
#' @description Creates an integrated Wiener process (IWP) model term for use in hierarchical models.
#'
#' @param x Variable name
#' @param p Order of the integrated Wiener process (default: 2)
#' @param ref_value Reference value for the basis (default: mean of knots)
#' @param knots Vector of knot locations
#' @param range Range of the data (optional)
#' @param init Initial value for theta parameter
#' @param lower Lower bound for theta parameter
#' @param upper Upper bound for theta parameter
#' @param parscale Parameter scale for optimization
#' @param boundary_is_random Whether boundary should be treated as random
#' @param include_poly Whether to include polynomial terms
#' @param prefix Optional prefix for term names
#'
#' @return A list containing iwp term object and optionally polynomial terms
#'
#' @examples
#' # Create an IWP term
#' iwp_term <- iwp(x = "pm25", knots = seq(0, 30, by = 5), ref_value = 10)

# IWP class definition
setClass("iwp",
  representation = representation(
  ),
  contains = "model",
  prototype = list(
  )
)

iwp <- function(
  x, p = 2,
  ref_value, knots, range = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE,
  prefix = ""
) {
  the_f <- formula(paste0("~ 0 + ", x))
  ref_value <- ref_align(ref_value, knots)
  result <- list()
  iwp_name <- paste(c(x, "iwp"), collapse = "_")
  result[[iwp_name]] <- new("iwp",
    term = x,
    formula = the_f,
    p.order = p,
    ref_value = ref_value,
    knots = knots,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1],
    type = factor("random", levels = .type_factor_levels)
  )
  if (include_poly) {
    poly_name <- paste(c(prefix, x, "poly"), collapse = "_")
    if (boundary_is_random) {
      result[[poly_name]] <- rpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    } else {
      result[[poly_name]] <- fpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    }
  }
  result
}

# Design matrix for iwp terms
setMethod("design", "iwp", function(term, data) {
  res <- local_poly(
    knots = term@knots - term@ref_value,
    refined_x = data[[term@term]] - term@ref_value, p = term@p.order
  )
  res <- methods::as(res, "TsparseMatrix")

  bnumPad <- formatC(1:ncol(res),
    width = max(ceiling(c(1, log10(ncol(res)))), na.rm = TRUE),
    flag = "0"
  )
  dimnames(res) <- list(
    rownames(data),
    paste(term@term, bnumPad, sep = "_")
  )
  res
})

# Precision matrix for iwp terms
setMethod("precision", "iwp", function(term, data) {
  result = as(compute_weights_precision(knots = term@knots), "TsparseMatrix")
  dimnames(result) = list(paste(
    term@term,
    1:ncol(result),
    sep="_"
  ))[c(1,1)]
  result
})

# Theta info for iwp terms
setMethod("theta_info", "iwp", function(term) {
  result <- data.frame(
    term = term@term, model = "iwp", 
    label = paste(c("iwp", term@term), collapse = "_"),
    order = NA, init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

# Beta info for iwp terms
setMethod("beta_info", "iwp", function(term) {
  # IWP terms don't have beta parameters
  return(NULL)
})

# Gamma info for iwp terms
setMethod("random_info", "iwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)
  order <- NA

 bnumPad <- formatC(basis,
    width = max(ceiling(c(1, log10(max(basis)))), na.rm = TRUE),
    flag = "0"
  )
 
  result <- expand.grid(
    term = term@term,
    model = "iwp",
    label = paste("model", term@term,sep="_"),
    by = NA,
    basis = basis,
    gamma_label <- paste(result$label, bnumPad, sep = "_"),
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result
})
