#' Build raw AD handle for one density shard
#'
#' Constructs CppAD tapes for a single shard. Merge shards with \code{c()} before
#' calling \code{\link{ad_fun}}.
#'
#' @param config Model configuration list (\code{beta}, \code{gamma}, \code{theta}, etc.).
#' @param kind \code{"observations"}, \code{"parameters"}, or \code{"random"}.
#' @param name Registered density name.
#' @param model An \code{ad_model} object.
#' @param precision For \code{kind = "random"}, list with \code{Q}. For
#'   \code{random_diagonal}, \code{Q} defaults to \code{rep(1, ncol(gamma_map))}
#'   when omitted.
#' @param package Backend package (currently only \code{"adlaplace"}).
#' @return External pointer of class \code{ad_fun_ptr}.
#' @export
ad_fun_ptr <- function(config,
                       kind,
                       name,
                       model,
                       precision = NULL,
                       package = "adlaplace") {
  if (missing(name) || is.null(name) || !nzchar(name)) {
    stop("`name` is required")
  }
  if (missing(model) || is.null(model)) {
    stop("`model` is required")
  }
  if (!is(model, "ad_model")) {
    stop("`model` must be an ad_model object")
  }
  if (!identical(package, "adlaplace")) {
    stop("ad_fun_ptr currently supports package = 'adlaplace' only")
  }

  validate_ad_model_maps(model, kind)

  if (identical(kind, "observations")) {
    validate_config_layout(model, config, kind)
    get_ad_fun_raw_obs(model, config, name)
  } else if (identical(kind, "parameters")) {
    validate_config_layout(model, config, kind)
    get_ad_fun_raw_parameters(model, config, name)
  } else if (identical(kind, "random")) {
    if (is(precision, "ad_model")) {
      stop("`precision` for kind = 'random' must be a list (e.g. list(Q = ...))")
    }
    if (is.null(precision)) {
      precision <- list()
    }
    validate_config_layout(model, config, kind)
    if (identical(name, "random_diagonal")) {
      if (is.null(precision$Q)) {
        precision$Q <- rep(1, ncol(model@gamma_map))
      }
      n_active <- Matrix::nnzero(model@gamma_map)
      if (length(precision$Q) != n_active) {
        stop(
          "length(precision$Q) (", length(precision$Q),
          ") must match active gamma_map entries (", n_active, ")"
        )
      }
    }
    get_ad_fun_raw_random(model, precision, config, name)
  } else {
    stop(
      "unknown `kind`: ", kind,
      " (expected 'observations', 'parameters', or 'random')"
    )
  }
}

#' @export
#' @method c ad_fun_ptr
c.ad_fun_ptr <- function(x, ...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return(x)
  }
  shards <- c(list(x), dots)
  if (!all(vapply(shards, function(s) is(s, "ad_fun_ptr"), logical(1L)))) {
    stop("all arguments must be ad_fun_ptr")
  }
  if (length(shards) == 1L) {
    return(shards[[1L]])
  }
  c_ad_fun_ptr(shards)
}
