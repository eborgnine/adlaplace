#' @include 000.R
#' Negative binomial observation model term
#'
#' @description Model term for the observation-level negative binomial log density
#' registered as \code{nbinom_obs}.
#' @name nbinom-class
#' @aliases nbinom
#' @docType class
#' @title Negative binomial observation term
#' @exportClass nbinom
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"nbinom_obs"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#' }
NULL

setClass(
  "nbinom",
  contains = "model",
  prototype = prototype(
    knots = numeric(0),
    ref_value = numeric(0),
    p.order = as.integer(1),
    init = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    parscale = numeric(0),
    type = factor("response", levels = .type_factor_levels),
    ad_fun = "nbinom_obs",
    ad_kind = "observations"
  )
)

#' Negative binomial observation term constructor
#'
#' @description Creates a response term wired to the \code{nbinom_obs} AD shard.
#' @rdname nbinom-class
#' @param x Outcome variable name.
#' @param init Initial value for the overdispersion parameter (log-SD scale).
#' @param lower Lower bound for the overdispersion parameter.
#' @param upper Upper bound for the overdispersion parameter.
#' @param parscale Parameter scale for optimization.
#' @return A \code{nbinom} object.
#' @export
nbinom <- function(x,
                   init = .my_theta_init,
                   lower = .my_theta_lower,
                   upper = .my_theta_upper,
                   parscale = .my_theta_parscale) {
  x <- strip_term_name(as.character(x))
  methods::new(
    "nbinom",
    term = x,
    label = paste(x, "nbinom_sd", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn nbinom-class Design method (not used for observation density).
#' @export
setMethod("design", "nbinom", function(term, data) {
  NULL
})

#' @describeIn nbinom-class Precision method (not used).
#' @export
setMethod("precision", "nbinom", function(term, data) {
  NULL
})

#' @describeIn nbinom-class Theta info for the overdispersion parameter.
#' @export
setMethod("theta_info", "nbinom", function(term) {
  data.frame(
    term = term@term,
    model = "nbinom",
    label = term@label,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    stringsAsFactors = FALSE
  )
})

#' @describeIn nbinom-class Beta info method (not used).
#' @export
setMethod("beta_info", "nbinom", function(term, data) {
  NULL
})

#' @describeIn nbinom-class Random info method (not used).
#' @export
setMethod("random_info", "nbinom", function(term, data) {
  NULL
})
