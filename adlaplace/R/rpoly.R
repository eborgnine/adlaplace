#' @include 000.R
#' Random Polynomial Model Term
#'
#' @description Creates and manages random polynomial model terms.
#' @name rpoly-class
#' @aliases rpoly
#' @docType class
#' @title Random Polynomial Model Term
#' @exportClass rpoly
#'
#' @section Methods:
#' The following methods are available for `rpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for rpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for rpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("rpoly",
         slots = list(
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

#' @rdname rpoly-class
#' @param x Variable name.
#' @param p Polynomial degree (default: 2).
#' @param ref_value Reference value for the polynomial.
#' @param sd Standard deviation for random effects.
#' @param term A `rpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A `rpoly` term object.
#' @export
rpoly <- function(x, p = 2, ref_value = 0, sd = Inf) {
  if (!missing(sd) && (length(sd) != 1 || sd <= 0)) {
    stop("sd must be a single positive numeric value")
  }
  if (!missing(p) && (length(p) != 1 || (round(p) != p) || p < 0)) {
    stop("p must be a single non-negative integer")
  }
  if (!missing(ref_value) && length(ref_value) != 1) {
    stop("ref_value must be a single value")
  }

  methods::new("rpoly",
    term = x,
    formula = stats::as.formula(paste0("~ 0 + ", x), env=new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    sd = rep_len(sd, p)
    # type is already set in prototype
  )
}

#' @describeIn rpoly-class Theta info method for rpoly objects
#' @export
#' @param term A `rpoly` term object.
#' @return NULL (random polynomial terms don't have theta info).
setMethod("theta_info", "rpoly", function(term) {
  # Random polynomial terms don't have theta info
  return(NULL)
})

#' @describeIn rpoly-class Design method for rpoly objects
#' @export
#' @param term A `rpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A design matrix for the random polynomial term.
setMethod("design", "rpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  D <- stats::poly(
    data[[term@term]] - term@ref_value,
    degree = term@p.order,
    raw = TRUE
  )
  D <- D[, 1:ncol(D), drop = FALSE]
  
  colnames(D) <- paste0(term@term, "_rpoly_", seq.int(1, length.out = term@p.order))
  D
})

#' @describeIn rpoly-class Precision method for rpoly objects
#' @export
#' @param term A `rpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A precision matrix for the random polynomial term.
setMethod("precision", "rpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  
  # Create precision matrix based on standard deviation
  p <- term@p.order
  sd_values <- rep_len(term@sd, p)
  
  # For random effects, we typically use identity matrix scaled by 1/sd^2
  precision_mat <- diag(1 / (sd_values^2), nrow = p, ncol = p)
  
  return(precision_mat)
})

#' @describeIn rpoly-class Random info method for rpoly objects
#' @export
#' @param term A `rpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A data frame containing random effects information for the random polynomial term.
setMethod("random_info", "rpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  
  the_colnames <- colnames(design(term, data))
  the_label <- paste(term@term, "rpoly", sep = "_")
  
  result <- data.frame(
    term = term@term,
    model = "rpoly",
    label = the_label,
    order = NA,
    random_label = the_colnames,
    sd = rep_len(term@sd, term@p.order)
  )
  
  return(result)
})

#' @describeIn rpoly-class Beta info method for rpoly objects
#' @export
#' @param term A `rpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return NULL (random polynomial terms don't have beta info).
setMethod("beta_info", "rpoly", function(term, data) {
  # Random polynomial terms don't have beta info
  return(NULL)
})

