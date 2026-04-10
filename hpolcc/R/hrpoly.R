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
#'
#' @return A hrpoly term object

#' @importFrom adlaplace fpoly rpoly
# HRPoly class definition
setClass("hrpoly",
         representation = representation(
           by_levels = "integer",
           by_labels = "character"
         ),
         contains = "model",
         prototype = list(
           by_levels = integer(0),
           by_labels = character(0),
           knots = numeric(0),
           type = factor("random", levels = adlaplace::.type_factor_levels)
         )
)

#' @export
hrpoly <- function(
  x,
  p = 1,
  ref_value = 0,
  by,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale
) {
  # Check all arguments are of length 1
  if (length(x) != 1) stop("x must be a single variable name")
  if (length(p) != 1) stop("p must be a single value")
  if (length(ref_value) != 1) stop("ref_value must be a single value")
  if (length(by) != 1) stop("by must be a single variable name")
  if (length(init) != 1) stop("init must be a single value")
  if (length(lower) != 1) stop("lower must be a single value")
  if (length(upper) != 1) stop("upper must be a single value")
  if (length(parscale) != 1) stop("parscale must be a single value")
  
  # Check that numeric arguments are valid
  if (p <= 0) stop("p must be positive")
  if (any(lower >= upper)) stop("lower bounds must be less than upper bounds")
  if (any(parscale <= 0)) stop("parscale must be positive")

  new("hrpoly",
    term = x,
    formula = as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    by = by,
    by_levels = integer(0),  # Will be set later when data is available
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
    # type is already set in prototype
  )
}


# Design matrix for hrpoly terms
setMethod("design", "hrpoly", function(term, data) {

  term = get_by_levels(term, data)

  a_base = adlaplace::design(as(term, "rpoly"), data)[,term@p.order]

  a_split = mapply(
    function(x, a_base, id) {
      the_i = which(x==id)
      data.frame(
        i = the_i,
        j_orig = id,
        x = a_base[the_i]
      )
    },
    id = term@by_levels,
    MoreArgs = list(x=data[[term@by]], a_base = a_base),
    SIMPLIFY=FALSE
  )
  a_df = do.call(rbind, a_split)
  a_df$j = match(a_df$j_orig, term@by_levels)

  result = Matrix::sparseMatrix(i=a_df$i, j=a_df$j, 
    x = a_df$x, dims = c(length(a_base), length(term@by_levels)),
    dimnames = list(NULL, paste0(
    term@term,
    "_hrpoly_",
    term@p.order,
    "_g",
    term@by_labels
  )))

  result
})

# Precision matrix for hrpoly terms
setMethod("precision", "hrpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  term = get_by_levels(term, data)

  result = Matrix::Diagonal(length(term@by_levels), 1)
  dimnames(result) = list(
    paste0(term@term, "_hrpoly_", term@p.order, "_g", term@by_labels)
  )[c(1,1)]
  result

})

# Theta info for hrpoly terms
setMethod("theta_info", "hrpoly", function(term) {
  result <- data.frame(
    term = term@term, model = "hrpoly", 
    label = paste(c(term@term, "hrpoly", term@p.order), collapse = "_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type
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
    label = paste(c(term@term, "hrpoly", term@p.order), collapse = "_"),
    by = term@by_levels,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result$by_labels <- term@by_labels[match(result$by, term@by_levels)]
  result$gamma_label <- paste0(result$label,  "_g", result$by_labels)
  result
})

methods::setAs("hrpoly", "rpoly", function(from) {
  rpoly(
    x = from@term,
    p = from@p.order,
    ref_value = from@ref_value,
    sd = 1
  )
})
