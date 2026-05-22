#' Build raw AD handle for one density shard
#'
#' Constructs CppAD tapes for a single shard (observation groups or one
#' single-density term). Merge several partial handles with \code{\link{combine}}.
#' For the full model with \code{hessian_map}, use \code{\link{get_ad_fun}}.
#'
#' @param data Model data list for this shard. For \code{kind = "obs"} (and
#'   obs-style singles like \code{neg_binom_extra}), use observation data
#'   (\code{y}, \code{ATp}, \code{XTp}, \code{theta_map}). For \code{kind = "single"}
#'   with \code{name = "random_diagonal"}, use random data (\code{Q},
#'   \code{theta_map}, \code{gamma_map}).
#' @param config Model configuration list passed to the backend builder.
#' @param kind Character scalar: \code{"obs"} or \code{"single"}.
#' @param name Registered density name (required; e.g. \code{"neg_binom_obs"},
#'   \code{"neg_binom_extra"}, \code{"random_diagonal"}).
#' @param package Backend package that registered the density; currently only
#'   \code{"adlaplace"} is supported.
#'
#' @return External pointer of class \code{adlaplace_handle_ptr}.
#'
#' @seealso \code{\link{combine}}, \code{\link{get_ad_fun}}
#' @export
get_ad_fun_raw <- function(data, config,
                           kind,
                           name,
                           package = "adlaplace") {
  if (missing(name) || is.null(name) || !nzchar(name)) {
    stop("`name` is required")
  }
  if (!identical(package, "adlaplace")) {
    stop("get_ad_fun_raw currently supports package = 'adlaplace' only")
  }
  if (kind == "obs") {
    get_ad_fun_raw_obs(data, config, name)
  } else if (kind == "single") {
    get_ad_fun_raw_single(data, config, name)
  } else {
    stop("unknown `kind`: ", kind, " (expected 'obs' or 'single')")
  }
}
