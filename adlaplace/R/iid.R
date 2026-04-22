#' IID Random Effects Term
#'
#' @description Creates and manages IID (independent and identically distributed) random effects terms.
#' @name iid-class
#' @aliases iid
#' @exportClass iid
#' @export
setClass("iid",
         slots = list(),
         contains = "model",
         prototype = prototype(
          ref_value = numeric(0),
          p.order = as.integer(0),
          knots = numeric(0),
          by = character(0),
          type = factor("random", levels = .type_factor_levels)
         )
)

#' @param x Variable name (typically a factor).
#' @param init Initial value for theta parameter.
#' @param lower Lower bound for theta parameter.
#' @param upper Upper bound for theta parameter.
#' @param parscale Parameter scale for optimization.
#' @return An `iid` term object (in a named list).
#' @rdname iid-class
#' @export
iid <- function(x,
                init = .my_theta_init,
                lower = .my_theta_lower,
                upper = .my_theta_upper,
                parscale = .my_theta_parscale) {
  result = list(methods::new("iid",
    term = x,
    formula = stats::as.formula(paste0("~ 0 + ", x), env=new.env()),
    init = init ,
    lower = lower ,
    upper = upper ,
    parscale = parscale
  ))
  names(result) = result[[1]]@term
  result
}

#' @describeIn iid-class Creates design matrix for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "iid", function(term, data){
  if(is.numeric(data[[term@term]])) {
    data[[term@term]] = factor(data[[term@term]])
  }
  result = methods::as(Matrix::sparse.model.matrix(term@formula, data), "TsparseMatrix")
  colnames(result) = gsub(paste0("^", term@term), paste0(term@term, "_iid_"), colnames(result))
  result
})

#' @describeIn iid-class Creates precision matrix for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "iid", function(term, data) {
  # Identity matrix for iid terms
  term_here = data[[term@term]]
  if(is.factor(term_here)) {
    n = nlevels(term_here)
    labels_here = levels(term_here)
  } else {
    labels_here = as.character(unique(term_here))
    n <- length(labels_here)
  }
  result = Matrix::Diagonal(n, 1)
  dimnames(result) = list(paste0(term@term , "_iid_", labels_here))[c(1,1)]
  result
})

#' @describeIn iid-class Extracts random effects information for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "iid", function(term, data) {
  
  result <- data.frame(
    term = term@term,
    model = "iid",
    label = paste(term@term, "iid", sep = "_"),
    by = NA,
    by_labels = NA,
    basis = as.character(sort(unique(data[[term@term]]))),
    order = NA
  )
  result$gamma_label = paste(result$label, result$basis, sep="_")

  result
})

#' @describeIn iid-class Extracts theta parameter information for iid term
#' @param term An iid term object
#' @export
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

#' @describeIn iid-class Extracts beta parameter information for iid term
#' @param term An iid term object
#' @export
setMethod("beta_info", "iid", function(term) {
  # IID terms don't have beta parameters
  return(NULL)
})
