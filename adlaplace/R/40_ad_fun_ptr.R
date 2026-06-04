#' Build raw AD handle for one density shard
#'
#' Constructs CppAD tapes for a single shard. The density kind and name come
#' from \code{data@ad_kind} and \code{data@ad_fun}. For \code{ad_kind = "random"},
#' \code{data@precision} is passed straight to the backend. For
#' \code{random_diagonal} it must be a numeric vector of diagonal precision
#' weights with \code{length == ncol(gamma_map)}. Missing precision means the
#' shard is not built (e.g. diffuse \code{rpoly} with \code{sd = Inf} via
#' \code{model_data()}); calling \code{ad_fun_ptr()} without it is an error.
#'
#' Merge handles for multiple shards with \code{c()} before calling
#' \code{\link{ad_fun}}.
#'
#' @param data An \code{ad_data} object with \code{ad_kind} and \code{ad_fun}
#'   slots set.
#' @param config Model configuration list (\code{beta}, \code{theta}, etc.).
#'   \code{gamma} is optional; when missing, zeros of length
#'   \code{nrow(gamma_map)} are used as AD tape seeds.
#' @return External pointer of class \code{ad_fun_ptr}.
#' @export
ad_fun_ptr <- function(data, config) {
  if (missing(data) || !is(data, "ad_data")) {
    stop("`data` must be an ad_data object", call. = FALSE)
  }
  kind <- data@ad_kind
  name <- data@ad_fun
  if (length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("data@ad_fun is required (density name)", call. = FALSE)
  }
  if (length(kind) != 1L || is.na(kind) || !nzchar(kind)) {
    stop("data@ad_kind is required (observations/parameters/random)", call. = FALSE)
  }

  validate_ad_data_maps(data, kind)
  config <- normalize_config_for_ptr(config, data, kind)
  validate_config_layout(data, config, kind)

  switch(
    kind,
    observations = get_ad_fun_raw_obs(data, config, name),
    parameters = get_ad_fun_raw_parameters(data, config, name),
    random = get_ad_fun_raw_random(data, data@precision, config, name),
    stop(
      "unknown ad_kind `", kind,
      "`; expected observations, parameters, or random",
      call. = FALSE
    )
  )
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

#' Deep copy of an \code{ad_fun_ptr} handle
#'
#' Returns a new \code{ad_fun_ptr} with independent C++ state (CppAD tapes and
#' sparsity patterns). Hessian templates are not copied; call \code{\link{ad_fun}}
#' on the clone to attach maps. The source handle is left unchanged, unlike
#' \code{\link{c.ad_fun_ptr}} which moves shards and clears the inputs.
#'
#' @param x An \code{ad_fun_ptr} object.
#' @return A new \code{ad_fun_ptr}.
#' @export
clone_ad_fun_ptr <- function(x) {
  if (!inherits(x, "ad_fun_ptr")) {
    stop("`x` must be an ad_fun_ptr object")
  }
  clone_ad_fun_ptr_(x)
}
