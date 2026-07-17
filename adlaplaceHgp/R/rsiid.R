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
    type = factor("random", levels = adlaplace::.type_factor_levels),
    ad_fun = "random_diagonal",
    ad_kind = "random"
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
#' @param init Initial values for beta parameters.
#' @param lower Lower bounds for beta parameters.
#' @param upper Upper bounds for beta parameters.
#' @param parscale Parameter scales for optimization.
#' @return A `rsiid` term object.
#' @export
rsiid <- function(
  x, mult, ref_mult = 0,
  init = NULL, lower = NULL,
  upper = NULL, parscale = NULL
) {
  if (!missing(ref_mult) && length(ref_mult) != 1) {
    stop("ref_mult must be a single value")
  }
  if (is.null(init)) init <- .my_theta_init
  if (is.null(lower)) lower <- .my_theta_lower
  if (is.null(upper)) upper <- .my_theta_upper
  if (is.null(parscale)) parscale <- .my_theta_parscale
  rsiid_label <- paste(c(x, mult, "rsiid"), collapse = "_")

  result <- list(methods::new("rsiid",
    term = x,
    label = rsiid_label,
    mult = mult,
    ref_mult = ref_mult,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init, lower = lower, upper = upper, parscale = parscale
  ))
  names(result) <- result[[1]]@term
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
    warning("missing values in ", term@mult)
  }

  term_iid = methods::as(term, "iid")
  a_matrix <- adlaplace::design(term_iid, data)
  a_matrix <- a_matrix * mult_vec

  colnames(a_matrix) <- gsub(
    "_iid_", paste0("_rsiid_", term@mult, "_"), colnames(a_matrix)
  )
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
  rownames(result) <- colnames(result) <- gsub(
    "_iid_", paste0("_rsiid_", term@mult, "_"), colnames(result)
  )
  result
})

# Theta info for rpoly terms
#' @describeIn rsiid-class Extracts theta parameter information for rsiid term
#' @param term A rsiid term object
#' @return A data frame containing theta parameter information for the random slope term
#' @export
setMethod("theta_info", "rsiid", function(term) {
  result <- data.frame(
    term = term@term,
    model = "rsiid",
    label = term@label,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

# Beta info for rpoly terms
#' @describeIn rsiid-class Extracts beta parameter information for rsiid term
#' @param term A rsiid term object
#' @param data A data frame containing the term variables (unused)
#' @return NULL (random slope polynomial terms don't have beta parameters)
#' @export
setMethod("beta_info", "rsiid", function(term, data) {
  NULL
})

# Gamma info for rpoly terms
#' @describeIn rsiid-class Extracts random effects information for rsiid term
#' @param term A rsiid term object
#' @param data A data frame containing the term variables
#' @return A data frame containing random effects information for the random slope polynomial term
#' @export
setMethod("random_info", "rsiid", function(term, data) {
  term_values <- data[[term@term]]
  basis_labels <- if (is.factor(term_values)) {
    levels(term_values)
  } else {
    as.character(sort(unique(term_values)))
  }

  result <- data.frame(
    term = term@term,
    model = "rsiid",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = basis_labels,
    order = NA
  )
  result$gamma_label <- paste0(
    term@term, "_rsiid_", term@mult, "_", result$basis
  )

  result
})
