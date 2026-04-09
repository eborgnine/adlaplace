#' Random Polynomial Model for Random Slope
#'
#' @description Creates a random polynomial model term for use in a random slope model.
#'
#' @param x Variable name
#' @param mult covariate to be mulitplied
#' @param p Polynomial degree (default: 2)
#' @param ref_value Reference value for the polynomial
#' @param ref_mult Reference value for the covariate
#' @param sd Standard deviation for random effects
#'
#' @return A rsrpoly term object

setClass("rsrpoly",
  representation = representation(
    mult = "numeric",
    ref_mult = "numeric",
    sd = "numeric"
  ),
  contains = "model",
  prototype = prototype(
    knots = numeric(0),
    by = character(0),
    sd = numeric(0),
    type = factor("random", levels = adlaplace::.type_factor_levels)
  )
)

#' @export
rsrpoly <- function(x, mult, p = 2, ref_value = 0, ref_mult = 0, sd = Inf) {
  # check sd is positive and length 1.  check p length 2 integer.  check ref value length 1
  if (!missing(sd) && (length(sd) != 1 || sd <= 0)) {
    stop("sd must be a single positive numeric value")
  }
  if (!missing(p) && (length(p) != 1 || (round(p) != p) || p < 0)) {
    stop("p must be a single non-negative integer")
  }
  if (!missing(ref_value) && length(ref_value) != 1) {
    stop("ref_value must be a single value")
  }
  if (!missing(ref_mult) && length(ref_mult) != 1) {
    stop("ref_mult must be a single value")
  }

  new("rsrpoly",
    term = x,
    mult = mult,
    formula = as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    ref_mult = ref_mult,
    sd = rep_len(sd, p)
  )
}

setMethod("design", "rsrpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }

  mult_vec <- data[[term@mult]] - term@ref_mult

  a_matrix <- poly(
    data[[term@term]] - term@ref_value,
    raw = TRUE, degree = term@p.order
  )
  a_matrix <- a_matrix[, 1:ncol(a_matrix), drop = FALSE]
  a_matrix <- a_matrix * mult_vec
  colnames(a_matrix) <- paste(
    term@term, term@mult, 
    "rsrpoly", 1:ncol(a_matrix), sep = "_"
    )
  a_matrix
})

# Precision matrix for rpoly terms
setMethod("precision", "rsrpoly", function(term, data) {
  the_sd <- term@sd
  names(the_sd) <- paste(term@term, "rsrpoly", seq.int(from = 1, length.out = term@p.order), sep = "_")
  Matrix::Diagonal(length(the_sd), the_sd^(-2), names = TRUE)
})

# Theta info for rpoly terms
setMethod("theta_info", "rsrpoly", function(term) {
  NULL
})

# Beta info for rpoly terms
setMethod("beta_info", "rsrpoly", function(term) {
  # Rpoly terms don't have beta parameters (random effects only)
  return(NULL)
})

# Gamma info for rpoly terms
setMethod("random_info", "rsrpoly", function(term, data) {
  order <- seq_len(term@p.order)

  result <- expand.grid(
    term = term@term,
    model = "rsrpoly",
    label = paste(c(term@term, term@mult, "rsrpoly"), collapse = "_"),
    by = NA,
    basis = NA,
    order = order,
    by_labels = NA,
    stringsAsFactors = FALSE
  )
  result$gamma_label <- paste(result$label, result$order, sep = "_")

  result
})
