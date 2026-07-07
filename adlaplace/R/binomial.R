#' @include 000.R
#' Bernoulli (binomial logit) observation model term
#'
#' @description Model term for the observation-level Bernoulli log density with
#' a logit link, registered as \code{binomial_obs}. For each observation with
#' linear predictor \code{eta = X beta + A gamma} the log density is
#' \code{y * eta - log(1 + exp(eta))}. The response \code{y} must be coded as
#' 0/1. There are no observation-level hyperparameters.
#' @name binomial-class
#' @aliases binomial
#' @docType class
#' @title Bernoulli observation term
#' @exportClass binomial
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"binomial_obs"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#' }
NULL

setClass(
  "binomial",
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
    ad_fun = "binomial_obs",
    ad_kind = "observations"
  )
)

#' Bernoulli observation term constructor
#'
#' @description Creates a response term wired to the \code{binomial_obs} AD
#' shard. A Bernoulli likelihood has no observation-level hyperparameters.
#'
#' When called with no response variable this function falls back to
#' \code{stats::binomial()}, so \code{glm(..., family = binomial)} and
#' \code{MASS::glmmPQL(..., family = binomial)} keep working with
#' \pkg{adlaplace} attached.
#'
#' @rdname binomial-class
#' @param x Outcome variable name (0/1 coded).
#' @param link Link passed to \code{stats::binomial()} in the fallback case.
#' @return A \code{binomial} object (or a \code{stats::family} when called
#'   with no response variable).
#' @export
binomial <- function(x, link = "logit") {
  if (missing(x)) {
    return(stats::binomial(link = link))
  }
  x <- strip_term_name(as.character(x))
  methods::new(
    "binomial",
    term = x,
    label = paste(x, "binomial", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env())
  )
}

#' @describeIn binomial-class Design method (not used for observation density).
#' @export
setMethod("design", "binomial", function(term, data) {
  NULL
})

#' @describeIn binomial-class Precision method (not used).
#' @export
setMethod("precision", "binomial", function(term, data) {
  NULL
})

#' @describeIn binomial-class Theta info (Bernoulli has no hyperparameters).
#' @export
setMethod("theta_info", "binomial", function(term) {
  NULL
})

#' @describeIn binomial-class Beta info method (not used).
#' @export
setMethod("beta_info", "binomial", function(term, data) {
  NULL
})

#' @describeIn binomial-class Random info method (not used).
#' @export
setMethod("random_info", "binomial", function(term, data) {
  NULL
})
