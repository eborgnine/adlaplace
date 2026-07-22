#' intercept Model Term
#'
#' @description Creates and manages intercept model terms for fixed effects.
#' @name intercept-class
#' @aliases intercept
#' @docType class
#' @title intercept Model Term
#' @exportClass intercept
#'
#' @section Methods:
#' The following methods are available for `intercept` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for intercept term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for intercept term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL


setClass("intercept",
  representation = representation(),
  contains = "model_term",
  prototype = list(
    name = "intercept",
    label = "intercept",
    formula = . ~ 1,
    ref_value = numeric(0),
    p.order = integer(0),
    knots = numeric(0),
    model_role = factor("fixed", levels = .model_role_levels),
    density = NA_character_,
    ad_kind = NA_character_
  )
)

#' @export
#' @rdname intercept-class
#' @param init Initial value for beta parameter (default: 0)
#' @param lower Lower bound for beta parameter (default: -Inf)
#' @param upper Upper bound for beta parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#' @return A intercept term object
#' @examples
#' # Create a intercept term
#' intercept_term <- intercept(init=2)
#' @export
intercept <- function(init = .my_beta_init,
                      lower = .my_beta_lower,
                      upper = .my_beta_upper,
                      parscale = .my_beta_parscale) {
  methods::new("intercept",
    label = "intercept",
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn intercept-class Creates design matrix for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "intercept", function(term, data) {
  matrix(1, nrow(data), ncol=1,
    dimnames = list(rownames(data),"intercept")
    )
})

#' @describeIn intercept-class Extracts beta parameter information for intercept term
#' @param term A intercept term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "intercept", function(term, data) {
  
  data.frame(
    term = "intercept",
    model = "intercept",
    label = term@label,
    order = NA,
    beta_label = "intercept",
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )
})

#' Linear Model Term
#'
#' @description Creates and manages linear model terms for fixed effects.
#' @name linear-class
#' @aliases linear
#' @docType class
#' @title Linear Model Term
#' @exportClass linear
#'
#' @section Methods:
#' The following methods are available for `linear` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for linear term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for linear term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("linear",
  representation = representation(),
  contains = "model_term",
  prototype = list(
    ref_value = numeric(0),
    p.order = integer(0),
    knots = numeric(0),
    model_role = factor("fixed", levels = .model_role_levels),
    density = NA_character_,
    ad_kind = NA_character_
  )
)

#' @export
#' @rdname linear-class
#' @param x Variable name (character) or symbol in a formula (e.g. \code{linear(x)} in
#'   \code{model_data()}).
#' @param init Initial value for beta parameter (default: 0)
#' @param lower Lower bound for beta parameter (default: -Inf)
#' @param upper Upper bound for beta parameter (default: Inf)
#' @param parscale Parameter scale for optimization (default: 1)
#' @return A linear term object
#' @examples
#' # Create a linear term
#' linear_term <- linear(x = "temperature")
#' @export
linear <- function(x,
                   init = .my_beta_init,
                   lower = .my_beta_lower,
                   upper = .my_beta_upper,
                   parscale = .my_beta_parscale) {
  if (is.symbol(x) || is.name(x)) {
    x <- as.character(x)
  } else if (!is.character(x)) {
    x <- as.character(x)
  }
  if (length(x) != 1L) {
    stop("x must be a single variable name", call. = FALSE)
  }
  methods::new("linear",
    name = x,
    label = paste(x, "linear", sep = "_"),
    formula = stats::as.formula(paste0("~ 0 + ", x)),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn linear-class Creates design matrix for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "linear", function(term, data) {
  res <- Matrix::sparse.model.matrix(term@formula, data, drop.unused.levels = FALSE)
  if (is.factor(data[[term@name]])) {
    res <- res[, -1, drop = FALSE]
  }
  colnames(res) <- paste0(term@name, "_linear_", colnames(res))
  res
})

#' @describeIn linear-class Extracts beta parameter information for linear term
#' @param term A linear term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "linear", function(term, data) {
  the_colnames <- colnames(design(term, data))
  the_label <- term@label

  result <- data.frame(
    term = term@name,
    model = "linear",
    label = the_label,
    order = NA,
    beta_label = the_colnames,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )

  return(result)
})

#' Fixed Polynomial Model Term
#'
#' @description Creates and manages fixed polynomial model terms.
#' @name fpoly-class
#' @aliases fpoly
#' @docType class
#' @exportClass fpoly
#'
#' @section Methods:
#' The following methods are available for `fpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for fpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for fpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("fpoly",
         slots = list(
           init = "numeric",
           lower = "numeric",
           upper = "numeric",
           parscale = "numeric"
         ),
         contains = "model_term",
         prototype = prototype(
           knots = numeric(0),
           init = numeric(0),
           lower = numeric(0),
           upper = numeric(0),
           parscale = numeric(0),
           model_role = factor("fixed", levels = .model_role_levels),
           density = NA_character_,
           ad_kind = NA_character_
         )
)

#' @rdname fpoly-class
#' @param x Variable name.
#' @param x Variable name.
#' @param p Polynomial degree (default: 2).
#' @param ref_value Reference value for the polynomial.
#' @param init Initial values for beta parameters.
#' @param lower Lower bounds for beta parameters.
#' @param upper Upper bounds for beta parameters.
#' @param parscale Parameter scales for optimization.
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A `fpoly` term object.
#' @export
#' @examples
#' # Create a fixed polynomial term
#' fpoly_term <- fpoly(x = "temp", p = 3, ref_value = 15)
fpoly <- function(x, p = 2, ref_value = 0,
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale) {
  methods::new("fpoly",
    name = x,
    label = paste(x, "fpoly", sep = "_"),
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    init = rep_len(init, p),
    lower = rep_len(lower, p),
    upper = rep_len(upper, p),
    parscale = rep_len(parscale, p)
  )
}

#' @describeIn fpoly-class Design method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A design matrix for the fixed polynomial term, or NULL if p.order is 0.
setMethod("design", "fpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }
  D <- stats::poly(
    data[[term@name]] - term@ref_value,
    degree = term@p.order,
    raw = TRUE
  )
  D <- D[, 1:ncol(D), drop = FALSE]
  seq_order <- seq.int(1, length.out = term@p.order)

  colnames(D) <- paste0(term@name, "_fpoly_", seq_order)
  D
})

#' @describeIn fpoly-class Beta info method for fpoly objects
#' @export
#' @param term A `fpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A data frame containing beta parameter information for the fixed polynomial term.
setMethod("beta_info", "fpoly", function(term, data) {
  the_colnames <- colnames(design(term, data))
  the_label <- term@label

  result <- data.frame(
    term = term@name,
    model = "fpoly",
    label = the_label,
    order = NA,
    beta_label = the_colnames,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )

  return(result)
})

