#' Independent and Identically Distributed Random Effects
#'
#' @description Creates an IID (independent and identically distributed) random effects term.
#'
#' @param x Variable name (typically a factor)
#' @param init Initial value for theta parameter
#' @param lower Lower bound for theta parameter
#' @param upper Upper bound for theta parameter
#' @param parscale Parameter scale for optimization
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

#' @export
iid <- function(x,
                init = .my_theta_init,
                lower = .my_theta_lower,
                upper = .my_theta_upper,
                parscale = .my_theta_parscale) {
  result = list(new("iid",
    term = x,
    formula = as.formula(paste0("~ 0 + ", x), env=new.env()),
    init = init ,
    lower = lower ,
    upper = upper ,
    parscale = parscale
  ))
  names(result) = result[[1]]@term
  result
}

# Design matrix for iid terms
setMethod("design", "iid", function(term, data){
  if(is.numeric(data[[term@term]])) {
    data[[term@term]] = factor(data[[term@term]])
  }
  result = as(Matrix::sparse.model.matrix(term@formula, data), "TsparseMatrix")
  colnames(result) = gsub(paste0("^", term@term), paste0(term@term, "_iid_"), colnames(result))
  result
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
    term = term@term, model = "iid", 
    label = paste(c(term@term,"iid"), 
    collapse = "_"),
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
    label = paste(term@term,"iid", sep="_"),
    by = NA,
    by_labels = NA,
    basis = sort(unique(data[[term@term]])),
    order = NA,
    stringsAsFactors = FALSE
  )
  result$gamma_label = paste(result$label, result$basis, sep="_")

  result
})
