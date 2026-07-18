#' Binomial observation model term
#'
#' @description Model term for the observation-level binomial log density with
#' a logit link, registered as \code{binomial_obs}. For each observation with
#' linear predictor \code{eta = X beta + A gamma}, trial count \code{N}, and
#' success count \code{y}, the log density contribution is
#' \code{y * eta - N * log(1 + exp(eta))} (up to a data-only constant). When
#' \code{size} is omitted, \code{N = 1} (Bernoulli). There are no
#' observation-level hyperparameters.
#' @name binomial-class
#' @aliases binomial
#' @docType class
#' @title Binomial observation term
#' @exportClass binomial
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"binomial_obs"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#' }
#' @section Slots:
#' \describe{
#'   \item{\code{size}}{Optional column name for per-observation trial counts.}
#' }
NULL

setClass(
  "binomial",
  slots = c(size = "character"),
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
    ad_kind = "observations",
    size = character(0)
  )
)

#' Binomial observation term constructor
#'
#' @description Creates a response term wired to the \code{binomial_obs} AD
#' shard. A Bernoulli likelihood (\code{size} omitted) has no observation-level
#' hyperparameters. With \code{size}, trial counts are read from that column
#' of \code{data} and stored as observation weights.
#'
#' When called with no response variable this function falls back to
#' \code{stats::binomial()}, so \code{glm(..., family = binomial)} and
#' \code{MASS::glmmPQL(..., family = binomial)} keep working with
#' \pkg{adlaplace} attached.
#'
#' @rdname binomial-class
#' @param x Outcome variable name (counts of successes; 0/1 when Bernoulli).
#' @param size Optional column name for the number of trials per observation.
#' @param link Link passed to \code{stats::binomial()} in the fallback case.
#' @return A \code{binomial} object (or a \code{stats::family} when called
#'   with no response variable).
#' @export
binomial <- function(x, size = NULL, link = "logit") {
  if (missing(x)) {
    return(stats::binomial(link = link))
  }
  x <- strip_term_name(as.character(x))
  size_col <- if (is.null(size)) {
    character(0)
  } else {
    strip_term_name(as.character(size))
  }
  methods::new(
    "binomial",
    term = x,
    label = paste(x, "binomial", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env()),
    size = size_col
  )
}

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
#' @param term A model term object (S4 methods).
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
#' @param term A model term object (S4 methods).
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

