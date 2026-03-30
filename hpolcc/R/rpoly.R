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
         prototype = prototype(
           knots = numeric(0),
           by = character(0),
           sd = numeric(0),
           type = factor("random", levels = .type_factor_levels)
         )
)

rpoly <- function(x, p = 2, ref_value = 0, sd = Inf) {
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

  new("rpoly",
    term = x,
    formula = formula(paste0("~ 0 + ", x)),
    p.order = as.integer(p),
    ref_value = ref_value,
    sd = rep_len(sd, p)
    # type is already set in prototype
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
  if(term@sd == Inf) {
    return(NULL)
  }
  Matrix::Diagonal(term@p.order, term@sd^(-2))
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
  order <- seq_len(term@p.order)
  basis <- NA

  result <- expand.grid(
    term = term@term,
    model = "rpoly",
    label = paste(c("rpoly", term@term), collapse = "_"),
    by = NA,
    basis = basis,
    order = order,
    by_labels = NA,
    stringsAsFactors = FALSE
  )
  result$gamma_label <- paste(result$label, result$order, sep = "_")

  result
})
