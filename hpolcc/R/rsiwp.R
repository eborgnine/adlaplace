

# IWP class definition
setClass("rsiwp",
  representation = representation(
    mult = "numeric"
  ),
  contains = "model",
  prototype = prototype(
    by = character(0),
    type = factor("random", levels = adlaplace::.type_factor_levels)
  )
)

#' @title Random Slope Integrated Wiener Process Term
#'
#' @description
#' Creates a Random Slope Integrated Wiener Process term.
#'
#' @param x Variable name
#' @param mult variable to multiply the IWP by
#' @param p Order of the integrated Wiener process (default: 2)
#' @param ref_value Reference value for the basis
#' @param ref_mult Reference value for the covariate
#' @param knots Vector of knot locations
#' @param range Range of the data (optional)
#' @param init Initial values for theta parameters
#' @param lower Lower bounds for theta parameters
#' @param upper Upper bounds for theta parameters
#' @param parscale Parameter scales for optimization
#' @param boundary_is_random Whether boundary should be treated as random
#' @param include_poly Whether to include polynomial terms
#'
#' @return A list containing iwp term object and optionally polynomial terms

#' @export
rsiwp <- function(
  x, 
  mult, 
  p = 2,
  ref_value = 0, 
  ref_mult= 0, 
  knots, 
  range = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE
) {
  # Check all arguments other than knots are length 1. knots must be length > 1
  if (length(x) != 1) stop("x must be a single variable name")
  if (length(p) != 1) stop("p must be a single value")
  if (length(ref_value) != 1) stop("ref_value must be a single value")
  if (length(ref_mult) != 1) stop("ref_value must be a single value")
  if (!is.null(range) && length(range) != 2) stop("range must be a vector of length 2")
  if (length(init) != 1) stop("init must be a single value")
  if (length(lower) != 1) stop("lower must be a single value")
  if (length(upper) != 1) stop("upper must be a single value")
  if (length(parscale) != 1) stop("parscale must be a single value")
  if (length(knots) < 2) stop("knots must have length >= 2")
  if (length(boundary_is_random) != 1) stop("boundary_is_random must be a single value")
  if (length(include_poly) != 1) stop("include_poly must be a single value")

  the_f <- as.formula(paste0("~ 0 + ", x), env=new.env())
  result <- list()
  iwp_name <- paste("rsiwp", x, sep = "_")

  ref_value = adlaplace::ref_align(ref_value, knots)


  result[[iwp_name]] <- new("rsiwp",
    term = x,
    mult = mult,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    ref_mult = ref_mult,
    knots = knots,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
    # type is already set in prototype, no need to repeat
  )
  if (include_poly) {
    poly_name <- paste(c(x, "rspoly"), collapse = "_")
    if (boundary_is_random) {
      result[[poly_name]] <- rsrpoly(
        x = x, mult=mult,
        p = p - 1,
        ref_value = ref_value, ref_mult = ref_mult
      )
    } else {
      result[[poly_name]] <- rsfpoly(
        x = x, mult=mult,
        p = p - 1,
        ref_value = ref_value,
        ref_mult = ref_mult
      )
    }
  }
  result
}

# Design matrix for iwp terms
setMethod("design", "rsiwp", function(term, data) {
  refined_x <- data[[term@term]] - term@ref_value
  mult_vec <- data[[term@mult]] - term@ref_mult


  basis <- adlaplace::local_poly(term@knots, refined_x, term@p.order)
  result <- basis[,1:ncol(basis),drop=FALSE]

  result = result * mult_vec

  knots_string = formatC(seq.int(ncol(result)), 
    width = ceiling(log10(ncol(result))), flag = "0")

  colnames(result) <- paste0(term@term, "_", term@mult, "_rsiwp_k", knots_string)
  result
})

# Precision matrix for iwp terms
setMethod("precision", "rsiwp", function(term, data) {
  result = Matrix::Matrix(adlaplace::compute_weights_precision(term@knots))
  
  knots_string = formatC(seq.int(nrow(result)), 
    width = ceiling(log10(nrow(result))), flag = "0")

  dimnames(result) = list(
    paste0(term@term, "_", term@mult, "_rsiwp_k", knots_string)
  )[c(1,1)]
  result
})

# Theta info for iwp terms
setMethod("theta_info", "rsiwp", function(term) {
  result <- data.frame(
    term = term@term, model = "rsiwp", 
    label = paste(c(term@term, term@mult, "rsiwp"), collapse = "_"),
    init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

# Beta info for iwp terms
setMethod("beta_info", "rsiwp", function(term) {
  # IWP terms don't have beta parameters
  return(NULL)
})

# Gamma info for iwp terms
setMethod("random_info", "rsiwp", function(term, data) {

  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "rsiwp",
    label = paste(c(term@term, term@mult, "rsiwp"), collapse = "_"),
    by = NA,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )

  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$by_labels <- NA  # iwp doesn't have by_labels
  result$gamma_label <- paste0(result$label, "_k", bnumPad)

  result
})
