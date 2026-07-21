#' @include 000.R
#' Skew-normal observation model term
#'
#' @description Model term for the observation-level skew-normal log density
#' registered as \code{skewnormal_obs}.
#' @name skewnormal-class
#' @docType class
#' @title Skew-normal observation term
#' @param term A \code{skewnormal} object.
#' @param data A data frame (unused for this observation term).
#' @exportClass skewnormal
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"skewnormal_obs"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#'   \item{\code{package}}{Character scalar \code{"adlaplaceExample"}.}
#' }
NULL

setClass(
  "skewnormal",
  contains = "model",
  prototype = prototype(
    knots = numeric(0),
    ref_value = numeric(0),
    p.order = as.integer(1),
    init = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    parscale = numeric(0),
    type = factor("response", levels = adlaplace:::.type_factor_levels),
    ad_fun = "skewnormal_obs",
    ad_kind = "observations",
    package = "adlaplaceExample"
  )
)

#' Skew-normal observation term constructor
#'
#' @description Creates a response term wired to the \code{skewnormal_obs} AD shard.
#' @param x Outcome variable name.
#' @param init Initial values for \code{omega} and \code{alpha} (log scale for omega).
#' @param lower Lower bounds for \code{omega} and \code{alpha}.
#' @param upper Upper bounds for \code{omega} and \code{alpha}.
#' @param parscale Parameter scales for optimization.
#' @param log Whether each theta is optimized on the log scale (length 2:
#'   \code{omega}, \code{alpha}).
#' @return A \code{skewnormal} object.
#' @export
skewnormal <- function(x,
                       init = c(0.1, 0.1),
                       lower = c(1e-9, -Inf),
                       upper = c(Inf, Inf),
                       parscale = c(1, 1),
                       log = c(TRUE, FALSE)) {
  x <- adlaplace::strip_term_name(as.character(x))
  log <- rep_len(log, 2L)
  methods::new(
    "skewnormal",
    term = x,
    label = paste(x, "skewnormal", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    log = log
  )
}

#' @describeIn skewnormal-class Design method (not used for observation density).
#' @export
setMethod("design", "skewnormal", function(term, data) {
  NULL
})

#' @describeIn skewnormal-class Precision method (not used).
#' @export
setMethod("precision", "skewnormal", function(term, data) {
  NULL
})

#' @describeIn skewnormal-class Theta info for skew-normal parameters.
#'
#' Scale \code{omega} uses \code{log = TRUE} (log scale when
#' \code{config$transform_theta} is not \code{FALSE}); shape \code{alpha} is untransformed.
#' @export
setMethod("theta_info", "skewnormal", function(term) {
  n <- length(term@init)
  data.frame(
    term = term@term,
    model = "skewnormal",
    label = paste0(term@label, c("_omega", "_alpha")),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    log = rep_len(term@log, n),
    stringsAsFactors = FALSE
  )
})

#' @describeIn skewnormal-class Beta info method (not used).
#' @export
setMethod("beta_info", "skewnormal", function(term, data) {
  NULL
})

#' @describeIn skewnormal-class Random info method (not used).
#' @export
setMethod("random_info", "skewnormal", function(term, data) {
  NULL
})
