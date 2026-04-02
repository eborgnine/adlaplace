#' Fixed Polynomial Model Term
#'
#' @description Creates a fixed polynomial model term.
#'
#' @param x Variable name
#' @param p Polynomial degree (default: 2)
#' @param ref_value Reference value for the polynomial
#' @param init Initial values for beta parameters
#' @param lower Lower bounds for beta parameters
#' @param upper Upper bounds for beta parameters
#' @param parscale Parameter scales for optimization
#'
#' @return A fpoly term object
#'
#' @examples
#' # Create a fixed polynomial term
#' fpoly_term <- fpoly(x = "temp", p = 3, ref_value = 15)

# FPoly class definition
setClass("fpoly",
         representation = representation(
         ),
         contains = "model",
         prototype = list(
                     by = character(0),
                     knots = numeric(0),
    type = factor("fixed", levels = .type_factor_levels)
         )
)

#' @export
#' @export
fpoly <- function(x, p = 2, ref_value = 0,
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale
                  ) {
  new("fpoly",
    term = x,
    formula = as.formula(paste0("~ 0 + ", x), env=new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    init = rep_len(init, p),
    lower = rep_len(lower, p),
    upper = rep_len(upper, p),
    parscale = rep_len(parscale, p)  )
}

# Design matrix for fpoly terms
setMethod("design", "fpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  D <- poly(
    data[[term@term]]-term@ref_value, 
    degree = term@p.order,
    raw = TRUE
    )
  D <- D[,1:ncol(D),drop=F]
  seq_order = seq.int(1, len=term@p.order)

  colnames(D) <- paste0(term@term, "_fpoly_", seq_order)
  D
})

# Precision matrix for fpoly terms
setMethod("precision", "fpoly", function(term, data) {
  # Fixed effects don't have precision matrices
  NULL
})

# Theta info for fpoly terms
setMethod("theta_info", "fpoly", function(term) {
  # not thetas for fpoly
  NULL
  })

# Beta info for fpoly terms
setMethod("beta_info", "fpoly", function(term, data) {
  the_label = paste(term@term, "fpoly", sep="_")
    seq_order = seq.int(1, len=term@p.order)

  result <- data.frame(
    term = term@term,
    model = "fpoly",
    label = the_label,
    order = seq_order,
    beta_label = paste(the_label, seq_order, sep="_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = as.character(term@type)
  )
  return(result)
})

# Gamma info for fpoly terms
setMethod("random_info", "fpoly", function(term, data) {
  # no gammas for fpoly
  NULL
})
