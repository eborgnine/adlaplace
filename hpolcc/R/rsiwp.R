#' Random Slope Integrated Wiener Process Term
#'
#' @description Creates and manages random slope integrated Wiener process (RSIWP) model terms.
#' @name rsiwp-class
#' @docType class
#' @exportClass rsiwp
#'
#' @section Methods:
#' The following methods are available for `rsiwp` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for RSIWP term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for RSIWP term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("rsiwp",
  representation = representation(
    mult = "character",
    ref_mult = "numeric"
  ),
  contains = "model",
  prototype = prototype(
    by = character(0),
    type = factor("random", levels = adlaplace::.type_factor_levels)
  )
)

# Register the coercion method properly
methods::setAs(
  "rsiwp", "iwp",
  function(from) {
    # Create a new iwp object with the same basic properties
    methods::new("iwp",
      term = from@term,
      formula = from@formula,
      knots = from@knots,
      ref_value = from@ref_value,
      p.order = from@p.order,
      lower = from@lower,
      upper = from@upper,
      parscale = from@parscale
    )
  }
)


#' Random Slope Integrated Wiener Process Term Constructor
#'
#' @description Creates a random slope integrated Wiener process (RSIWP) model term.
#' @name rsiwp
#' @param x Variable name.
#' @param mult Variable to multiply the IWP by.
#' @param p Order of the integrated Wiener process (default: 2).
#' @param ref_value Reference value for the basis.
#' @param ref_mult Reference value for the covariate.
#' @param knots Vector of knot locations.
#' @param range Range of the data (optional).
#' @param init Initial values for theta parameters.
#' @param lower Lower bounds for theta parameters.
#' @param upper Upper bounds for theta parameters.
#' @param parscale Parameter scales for optimization.
#' @param boundary_is_random Whether boundary should be treated as random.
#' @param include_poly Whether to include polynomial terms.
#' @return A list containing the `rsiwp` term object and optionally polynomial terms.
#' @export
rsiwp <- function(
  x,
  mult,
  p = 2,
  ref_value = 0,
  ref_mult = 0,
  knots,
  range = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE
) {
  # Check all arguments other than knots are length 1. knots must be length > 1
  if (length(x) != 1) stop("x must be a single variable name")
  if (length(p) != 1) stop("p must be a single value")
  if (length(ref_value) != 1) stop("ref_value must be a single value")
  if (length(ref_mult) != 1) stop("ref_value must be a single value")
  if (!is.null(range) && length(range) != 2) stop("range must be a vector of length 2")
  if (length(init) != 1) stop("init must be a single value")
  if (length(lower) != 1) stop("lower must be a single value")
  if (length(upper) != 1) stop("upper must be a single value")
  if (length(parscale) != 1) stop("parscale must be a single value")
  if (length(knots) < 2) stop("knots must have length >= 2")
  if (length(boundary_is_random) != 1) stop("boundary_is_random must be a single value")
  if (length(include_poly) != 1) stop("include_poly must be a single value")

  the_f <- stats::as.formula(paste0("~ 0 + ", x), env = new.env())
  result <- list()
  iwp_name <- paste("rsiwp", x, sep = "_")

  ref_value <- adlaplace::ref_align(ref_value, knots)
  knots_ref <- knots - ref_value


  result[[iwp_name]] <- methods::new("rsiwp",
    term = x,
    mult = mult,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    ref_mult = ref_mult,
    knots = knots_ref,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
    # type is already set in prototype, no need to repeat
  )
  if(p==1) {
    include_poly = FALSE
  }
  if (include_poly) {
    if (boundary_is_random) {
    poly_name <- paste(c(x, "rsrpoly"), collapse = "_")
      result[[poly_name]] <- rsrpoly(
        x = x, mult = mult,
        p = p - 1,
        ref_value = ref_value, ref_mult = ref_mult
      )
    result[[paste0(x, "_linear")]] = adlaplace::rpoly(
      x, p=1, ref_value = ref_value
    )      
    } else {
    poly_name <- paste(c(x, "rsfpoly"), collapse = "_")
      result[[poly_name]] <- rsfpoly(
        x = x, mult = mult,
        p = p - 1,
        ref_value = ref_value,
        ref_mult = ref_mult
      )
    result[[paste0(x, "_linear")]] = adlaplace::fpoly(
      x, p=1, ref_value = ref_value
    )      
    }
  }
  result
}
#' @describeIn rsiwp-class Creates design matrix for RSIWP term
#' @param term An rsiwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "rsiwp", function(term, data) {
  design_iwp <- adlaplace::design(methods::as(term, "iwp"), data)

  mult_vec <- data[[term@mult]] - term@ref_mult
  result <- design_iwp * mult_vec

  knots_string <- gsub(".*_k", "", colnames(design_iwp))

  colnames(result) <- paste0(
    term@term, "_", term@mult,
    "_rsiwp_k", knots_string
  )
  result
})

#' @describeIn rsiwp-class Creates precision matrix for RSIWP term
#' @param term An rsiwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "rsiwp", function(term, data) {
  term_iwp = methods::as(term, "iwp")
  result <- Matrix::Matrix(adlaplace::precision(term_iwp))

  knots_string <- formatC(
    seq.int(nrow(result)),
    width = ceiling(log10(nrow(result))),
    flag = "0"
  )

  dimnames(result) <- list(
    paste0(term@term, "_", term@mult, "_rsiwp_k", knots_string)
  )[c(1, 1)]
  result
})

#' @describeIn rsiwp-class Extracts theta parameter information for RSIWP term
#' @param term An rsiwp term object
#' @export
setMethod("theta_info", "rsiwp", function(term) {
  result <- data.frame(
    term = term@term, model = "rsiwp",
    label = paste(c(term@term, term@mult, "rsiwp"), collapse = "_"),
    init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

#' @describeIn rsiwp-class Extracts beta parameter information for RSIWP term
#' @param term An rsiwp term object
#' @export
setMethod("beta_info", "rsiwp", function(term) {
  # IWP terms don't have beta parameters
  return(NULL)
})

#' @describeIn rsiwp-class Extracts random effects information for RSIWP term
#' @param term An rsiwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "rsiwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "rsiwp",
    label = paste(c(term@term, term@mult, "rsiwp"), collapse = "_"),
    by = NA,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )

  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$by_labels <- NA # iwp doesn't have by_labels
  result$gamma_label <- paste0(result$label, "_k", bnumPad)

  result
})
