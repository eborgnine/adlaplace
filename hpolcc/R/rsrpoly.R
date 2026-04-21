#' Random Polynomial Model for Random Slope
#'
#' @description Creates a random polynomial model term for use in a random slope model.
#' @name rsrpoly-class
#' @docType class
#' @exportClass rsrpoly
#'
#' @section Methods:
#' The following methods are available for `rsrpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for rsrpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for rsrpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("rsrpoly",
  representation = representation(
    mult = "character",
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

#' @rdname rsrpoly-class
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

  methods::new("rsrpoly",
    term = x,
    mult = mult,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    ref_mult = ref_mult,
    sd = rep_len(sd, p)
  )
}

#' @describeIn rsrpoly-class Creates design matrix for rsrpoly term
#' @param term A rsrpoly term object
#' @param data A data frame containing the term variables
#' @return A design matrix for the random slope polynomial term, or NULL if p.order is 0
#' @export
setMethod("design", "rsrpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }

  mult_vec <- data[[term@mult]] - term@ref_mult

  if(any(is.na(mult_vec))) {
    warning("missing values in", term@mult)
  }

  a_matrix <- stats::poly(
    data[[term@term]] - term@ref_value,
    raw = TRUE, degree = term@p.order
  )
  a_matrix <- a_matrix[, 1:ncol(a_matrix), drop = FALSE]
  a_matrix <- a_matrix * mult_vec
  colnames(a_matrix) <- paste(
    term@term, term@mult,
    "rsrpoly", 1:ncol(a_matrix),
    sep = "_"
  )
  a_matrix
})

# Precision matrix for rpoly terms
#' @describeIn rsrpoly-class Creates precision matrix for rsrpoly term
#' @param term A rsrpoly term object
#' @param data A data frame containing the term variables
#' @return A precision matrix for the random slope polynomial term
#' @export
setMethod("precision", "rsrpoly", function(term, data) {
  the_sd <- term@sd
  names(the_sd) <- paste(
    term@term, term@mult, "rsrpoly",
    seq.int(from = 1, length.out = term@p.order),
    sep = "_"
  )
  Matrix::Diagonal(length(the_sd), the_sd^(-2), names = TRUE)
})

# Theta info for rpoly terms
#' @describeIn rsrpoly-class Extracts theta parameter information for rsrpoly term
#' @param term A rsrpoly term object
#' @return NULL (random slope polynomial terms don't have theta parameters)
#' @export
setMethod("theta_info", "rsrpoly", function(term) {
  NULL
})

# Beta info for rpoly terms
#' @describeIn rsrpoly-class Extracts beta parameter information for rsrpoly term
#' @param term A rsrpoly term object
#' @return NULL (random slope polynomial terms don't have beta parameters)
#' @export
setMethod("beta_info", "rsrpoly", function(term) {
  # Rpoly terms don't have beta parameters (random effects only)
  return(NULL)
})

# Gamma info for rpoly terms
#' @describeIn rsrpoly-class Extracts random effects information for rsrpoly term
#' @param term A rsrpoly term object
#' @param data A data frame containing the term variables
#' @return A data frame containing random effects information for the random slope polynomial term
#' @export
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
