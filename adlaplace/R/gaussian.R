#' @include 000.R
#' Gaussian observation model term
#'
#' @description Model term for the observation-level Gaussian log density
#' registered as \code{gaussian_obs}. This is the default observation term:
#' a bare response symbol on the left-hand side of a formula (e.g.
#' \code{y ~ x}) is coerced to \code{gaussian(y)} by
#' \code{\link{collect_terms}}.
#' @name gaussian-class
#' @aliases gaussian
#' @docType class
#' @title Gaussian observation term
#' @exportClass gaussian
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"gaussian_obs"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#' }
NULL

setClass(
  "gaussian",
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
    ad_fun = "gaussian_obs",
    ad_kind = "observations"
  )
)

#' Gaussian observation term constructor
#'
#' @description Creates a response term wired to the \code{gaussian_obs} AD
#' shard, with a single residual standard deviation parameter (log scale
#' during optimization).
#'
#' When called with no arguments this function falls back to
#' \code{stats::gaussian()}, so \code{glm(..., family = gaussian)} keeps
#' working with \pkg{adlaplace} attached.
#'
#' @rdname gaussian-class
#' @param x Outcome variable name.
#' @param init Initial value for the residual standard deviation.
#' @param lower Lower bound for the residual standard deviation.
#' @param upper Upper bound for the residual standard deviation.
#' @param parscale Parameter scale for optimization.
#' @return A \code{gaussian} object (or a \code{stats::family} when called
#'   with no arguments).
#' @export
gaussian <- function(x,
                     init = 1,
                     lower = .my_theta_lower,
                     upper = .my_theta_upper,
                     parscale = .my_theta_parscale) {
  if (missing(x)) {
    return(stats::gaussian())
  }
  x <- strip_term_name(as.character(x))
  methods::new(
    "gaussian",
    term = x,
    label = paste(x, "gaussian_sd", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
  )
}

#' @describeIn gaussian-class Design method (not used for observation density).
#' @export
setMethod("design", "gaussian", function(term, data) {
  NULL
})

#' @describeIn gaussian-class Precision method (not used).
#' @export
setMethod("precision", "gaussian", function(term, data) {
  NULL
})

#' @describeIn gaussian-class Theta info for the residual standard deviation.
#' @export
setMethod("theta_info", "gaussian", function(term) {
  data.frame(
    term = term@term,
    model = "gaussian",
    label = term@label,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    stringsAsFactors = FALSE
  )
})

#' @describeIn gaussian-class Beta info method (not used).
#' @export
setMethod("beta_info", "gaussian", function(term, data) {
  NULL
})

#' @describeIn gaussian-class Random info method (not used).
#' @export
setMethod("random_info", "gaussian", function(term, data) {
  NULL
})
