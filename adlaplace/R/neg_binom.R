#' @include 000.R
#' Negative binomial observation model term
#'
#' @description Model term for the observation-level negative binomial log density
#' registered as \code{neg_binom_obs}.
#' @name neg_binom-class
#' @aliases neg_binom
#' @docType class
#' @title Negative binomial observation term
#' @exportClass neg_binom
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"neg_binom_obs"}.}
#'   \item{\code{as_kind}}{Character scalar \code{"observations"}.}
#' }
NULL

setClass(
  "neg_binom",
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
    ad_fun = "neg_binom_obs",
    as_kind = "observations"
  )
)

#' Negative binomial observation term constructor
#'
#' @description Creates a response term wired to the \code{neg_binom_obs} AD shard.
#' @rdname neg_binom-class
#' @param x Outcome variable name.
#' @return A \code{neg_binom} object.
#' @export
neg_binom <- function(x) {
  methods::new(
    "neg_binom",
    term = as.character(x),
    formula = stats::as.formula(paste(x, "~."), env = new.env())
  )
}

#' @describeIn neg_binom-class Design method (not used for observation density).
#' @export
setMethod("design", "neg_binom", function(term, data) {
  NULL
})

#' @describeIn neg_binom-class Precision method (not used).
#' @export
setMethod("precision", "neg_binom", function(term, data) {
  NULL
})

#' @describeIn neg_binom-class Theta info method (not used).
#' @export
setMethod("theta_info", "neg_binom", function(term) {
  NULL
})

#' @describeIn neg_binom-class Beta info method (not used).
#' @export
setMethod("beta_info", "neg_binom", function(term, data) {
  NULL
})

#' @describeIn neg_binom-class Random info method (not used).
#' @export
setMethod("random_info", "neg_binom", function(term, data) {
  NULL
})

#' Negative binomial extra (single-shard) model term
#'
#' @description Model term for the single-density negative binomial piece registered
#' as \code{neg_binom_extra}.
#' @name neg_binom_extra-class
#' @aliases neg_binom_extra
#' @docType class
#' @title Negative binomial extra term
#' @exportClass neg_binom_extra
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"neg_binom_extra"}.}
#'   \item{\code{as_kind}}{Character scalar \code{"parameters"}.}
#' }
NULL

setClass(
  "neg_binom_extra",
  contains = "model",
  prototype = prototype(
    term = character(0),
    formula = formula(),
    knots = numeric(0),
    ref_value = numeric(0),
    p.order = integer(0),
    type = factor("response", levels = .type_factor_levels),
    ad_fun = "neg_binom_extra",
    as_kind = "parameters"
  )
)

#' Negative binomial extra term constructor
#'
#' @description Creates a term wired to the \code{neg_binom_extra} AD shard.
#' @rdname neg_binom_extra-class
#' @return A \code{neg_binom_extra} object.
#' @export
neg_binom_extra <- function() {
  methods::new("neg_binom_extra")
}

#' @describeIn neg_binom_extra-class Design method (not used).
#' @export
setMethod("design", "neg_binom_extra", function(term, data) {
  NULL
})

#' @describeIn neg_binom_extra-class Precision method (not used).
#' @export
setMethod("precision", "neg_binom_extra", function(term, data) {
  NULL
})

#' @describeIn neg_binom_extra-class Theta info method (not used).
#' @export
setMethod("theta_info", "neg_binom_extra", function(term) {
  NULL
})

#' @describeIn neg_binom_extra-class Beta info method (not used).
#' @export
setMethod("beta_info", "neg_binom_extra", function(term, data) {
  NULL
})

#' @describeIn neg_binom_extra-class Random info method (not used).
#' @export
setMethod("random_info", "neg_binom_extra", function(term, data) {
  NULL
})
