#' @include 000.R
#' Dirichlet-multinomial case-crossover response term
#'
#' @description Model term for the case-crossover Dirichlet-multinomial log density
#' registered as \code{dirichlet_multinomial}.
#' @name dirichlet_multinom-class
#' @aliases dirichlet_multinom
#' @docType class
#' @title Dirichlet-multinomial response term
#' @exportClass dirichlet_multinom
#'
#' @slot by Character vector of stratification variables for the case-crossover design.
#'
#' @section Slots (inherited from \code{model}):
#' \describe{
#'   \item{\code{ad_fun}}{Character scalar \code{"dirichlet_multinomial"}.}
#'   \item{\code{ad_kind}}{Character scalar \code{"observations"}.}
#'   \item{\code{package}}{Character scalar \code{"hpolcc"}.}
#' }
NULL

setClass(
  "dirichlet_multinom",
  slots = c(by = "character"),
  contains = "model",
  prototype = prototype(
    knots = numeric(0),
    ref_value = numeric(0),
    p.order = as.integer(1),
    init = 1e-3,
    lower = 0,
    upper = Inf,
    parscale = 1,
    type = factor("response", levels = adlaplace:::.type_factor_levels),
    ad_fun = "dirichlet_multinomial",
    ad_kind = "observations",
    package = "hpolcc",
    by = character(0)
  )
)

#' Dirichlet-multinomial case-crossover response term
#'
#' @param x Outcome variable name.
#' @param by Character vector of stratification variables (e.g. \code{c("year", "month", "dow")}).
#' @param init Initial value for overdispersion SD \code{tau}.
#' @param lower Lower bound for \code{tau}.
#' @param upper Upper bound for \code{tau}.
#' @param parscale Parameter scale for optimization.
#' @return A \code{dirichlet_multinom} object.
#' @export
dirichlet_multinom <- function(x,
                               by,
                               init = 1e-3,
                               lower = 0,
                               upper = Inf,
                               parscale = 1) {
  x <- adlaplace:::strip_term_name(as.character(x))
  by <- as.character(by)
  methods::new(
    "dirichlet_multinom",
    term = x,
    label = paste(x, "dirichlet_multinom", sep = "_"),
    formula = stats::as.formula(paste(x, "~."), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    by = by
  )
}

#' Case-crossover stratum map for AD shards
#'
#' @param term A \code{dirichlet_multinom} term object.
#' @param data Analysis-ready data frame.
#' @return Sparse \code{ngCMatrix} mapping observations to strata columns.
#' @describeIn dirichlet_multinom-class Case-crossover stratum map for AD shards.
#' @export
setMethod("elgm_matrix", "dirichlet_multinom", function(term, data) {
  if (!length(term@by)) {
    stop("dirichlet_multinom(..., by = ...) is required", call. = FALSE)
  }
  cc_design <- ccDesign(strat_vars = term@by)
  setStrata(cc_design, data.table::as.data.table(data), outcome = term@term)
})

#' @describeIn dirichlet_multinom-class Design method (not used).
#' @export
setMethod("design", "dirichlet_multinom", function(term, data) {
  NULL
})

#' @describeIn dirichlet_multinom-class Precision method (not used).
#' @export
setMethod("precision", "dirichlet_multinom", function(term, data) {
  NULL
})

#' @describeIn dirichlet_multinom-class Theta info for overdispersion \code{tau}.
#' @export
setMethod("theta_info", "dirichlet_multinom", function(term) {
  data.frame(
    term = term@term,
    model = "dirichlet_multinom",
    label = paste0(term@label, "_tau"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    transform = TRUE,
    stringsAsFactors = FALSE
  )
})

#' @describeIn dirichlet_multinom-class Beta info method (not used).
#' @export
setMethod("beta_info", "dirichlet_multinom", function(term, data) {
  NULL
})

#' @describeIn dirichlet_multinom-class Random info method (not used).
#' @export
setMethod("random_info", "dirichlet_multinom", function(term, data) {
  NULL
})
