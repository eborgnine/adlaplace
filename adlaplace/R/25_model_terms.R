#' By-Group Class
#'
#' @description A class for grouping information in hierarchical model terms.
#' Contains the grouping variable name and its levels/labels.
#' @slot term Character vector of length 1, the grouping variable name
#' @slot levels Character vector of unique levels for the grouping variable
#' @slot labels Character vector of formatted labels for the levels
#' @exportClass by_group
setClass(
  "by_group",
  slots = list(
    term = "character",
    levels = "character",
    labels = "character"
  ),
  prototype = prototype(
    term = character(0),
    levels = character(0),
    labels = character(0)
  )
)

#' Base Model Class
#'
#' @description The base S4 class that all specific model term classes inherit from.
#' Subclasses that need hierarchical grouping should include a \code{by} slot
#' (either \code{character} or \code{by_group} type) and optionally \code{by_levels}
#' and \code{by_labels} slots.
#' @slot term Character vector of length 1, the term name
#' @slot label Character scalar term label reused across metadata builders.
#' @slot formula Formula object for the term
#' @slot knots Numeric vector of knot locations (for spline terms)
#' @slot ref_value Numeric reference value for centering
#' @slot p.order Integer polynomial order
#' @slot init Numeric vector of initial values
#' @slot lower Numeric vector of lower bounds
#' @slot upper Numeric vector of upper bounds
#' @slot parscale Numeric vector of parameter scales for optimization
#' @slot type Factor indicating term type ("fixed", "random", or "response")
#' @slot ad_fun Registered AD density name for \code{ad_fun_ptr()}, or \code{NA}
#'   when the term does not map to a density shard.
#' @slot ad_kind Shard kind passed to \code{ad_fun_ptr()} (\code{"observations"},
#'   \code{"parameters"}, \code{"random"}, or \code{NA}).
#' @slot package Package whose \code{.so} records AD tapes for this term
#'   (default \code{"adlaplace"} for built-in densities).
#' @exportClass model
setClass(
  "model",
  slots = list(
    term = "character",
    label = "character",
    formula = "formula",
    knots = "numeric",
    ref_value = "numeric",
    p.order = "integer",
    init = "numeric",
    lower = "numeric",
    upper = "numeric",
    parscale = "numeric",
    type = "factor",
    ad_fun = "character",
    ad_kind = "character",
    package = "character"
  ),
  prototype = prototype(
    term = character(0),
    label = character(0),
    formula = formula(),
    knots = numeric(0),
    ref_value = numeric(0),
    p.order = as.integer(1),
    init = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    parscale = numeric(0),
    type = factor("fixed", levels = .type_factor_levels),
    ad_fun = NA_character_,
    ad_kind = NA_character_,
    package = "adlaplace"
  )
)

#' Model Term Generics
#'
#' @description Generic functions for extracting model matrices and parameter information from model terms.
#' @param term A model term object.
#' @param data A data frame containing the variables.
#' @return Method-specific return values.
#' @name model-generics
NULL

#' @rdname model-generics
#' @export
setGeneric("design", function(term, data) standardGeneric("design"))
#' @rdname model-generics
#' @export
setGeneric("precision", function(term, data) standardGeneric("precision"))
#' @rdname model-generics
#' @export
setGeneric("theta_info", function(term) standardGeneric("theta_info"))
#' @rdname model-generics
#' @export
setGeneric("beta_info", function(term, data) standardGeneric("beta_info"))
#' @rdname model-generics
#' @export
setGeneric("random_info", function(term, data) standardGeneric("random_info"))
