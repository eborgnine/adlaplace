#' independent Random Slope
#'
#' @description Creates an iid random slope model term
#' @name rsiid-class
#' @docType class
#' @exportClass rsiid
#'
#' @section Methods:
#' The following methods are available for `rsiid` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for rsiid term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for rsiid term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("rsiid",
  representation = representation(
    mult = "character",
    ref_mult = "numeric"
  ),
  contains = "model",
  prototype = prototype(
    type = factor("random", levels = adlaplace::.type_factor_levels)
  )
)

# Register the coercion method properly
methods::setAs(
  "rsiid", "iid",
  function(from) {
    methods::new("iid",
      term = from@term,
      formula = from@formula,
      init = from@init,
      lower = from@lower,
      upper = from@upper,
      parscale = from@parscale
    )
  }
)

#' @rdname rsiid-class
#' @param x Variable name.
#' @param mult Variable to multiply the polynomial by.
#' @param ref_mult Reference value for the covariate.
#' @param sd Standard deviation for random effects.
#' @return A `rsiid` term object.
#' @export
rsiid <- function(
  x, mult, ref_mult = 0,
  init = .my_theta_init, lower = .my_theta_lower,
  upper = .my_theta_upper, parscale = .my_theta_parscale
) {
  if (!missing(ref_mult) && length(ref_mult) != 1) {
    stop("ref_mult must be a single value")
  }
  result = list(methods::new("rsiid",
    term = x,
    mult = mult,
    ref_mult = ref_mult,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init, lower = lower, upper = upper, parscale = parscale
    ))
  names(result) = result[[1]]@term
  result
}

#' @describeIn rsiid-class Creates design matrix for rsiid term
#' @param term A rsiid term object
#' @param data A data frame containing the term variables
#' @return A design matrix for the random slope polynomial term, or NULL if p.order is 0
#' @export
setMethod("design", "rsiid", function(term, data) {
  mult_vec <- data[[term@mult]] - term@ref_mult

  if (any(is.na(mult_vec))) {
    warning("missing values in", term@mult)
  }

  term_iid = methods::as(term, "iid")
  a_matrix <- adlaplace::design(term_iid, data)
  a_matrix <- a_matrix * mult_vec

  colnames(a_matrix) <- gsub("_iid_", "_rsiid_", colnames(a_matrix))
  a_matrix
})

# Precision matrix for rpoly terms
#' @describeIn rsiid-class Creates precision matrix for rsiid term
#' @param term A rsiid term object
#' @param data A data frame containing the term variables
#' @return A precision matrix for the random slope polynomial term
#' @export
setMethod("precision", "rsiid", function(term, data) {
  result <- adlaplace::precision(methods::as(term, "iid"), data)
  rownames(result) <- colnames(result) <- 
    gsub("_iid_", "_rsiid_", colnames(result))
  result
})

# Theta info for rpoly terms
#' @describeIn rsiid-class Extracts theta parameter information for rsiid term
#' @param term A rsiid term object
#' @return A data frame containing theta parameter information for the random slope term
#' @export
setMethod("theta_info", "rsiid", function(term) {
  result <- adlaplace::theta_info(methods::as(term, "iid"))
  result$model = "rsiid"
  result$label = gsub("_iid$", "_rsiid", result$label)
  result
})

# Beta info for rpoly terms
#' @describeIn rsiid-class Extracts beta parameter information for rsiid term
#' @param term A rsiid term object
#' @return NULL (random slope polynomial terms don't have beta parameters)
#' @export
setMethod("beta_info", "rsiid", function(term) {

  return(NULL)
})

# Gamma info for rpoly terms
#' @describeIn rsiid-class Extracts random effects information for rsiid term
#' @param term A rsiid term object
#' @param data A data frame containing the term variables
#' @return A data frame containing random effects information for the random slope polynomial term
#' @export
setMethod("random_info", "rsiid", function(term, data) {

  result = adlaplace::random_info(methods::as(term, "iid"), data)
  result$model = "rsiid"
  result$label = gsub("_iid$", "_rsiid", result$label)
  result$gamma_label = gsub("_iid_", "_rsiid_", result$gamma_label)
  result
})
