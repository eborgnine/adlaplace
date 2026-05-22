#' Build backend AD function handle
#'
#' Safe frontend for constructing AD handles. This is the only public
#' \code{get_ad_fun()} entry point; backend packages should expose
#' \code{getAdFun_r()}.
#'
#' @param data Model data list (unused for handle-only attachment; retained for
#'   API compatibility and future \code{model_setup} integration).
#' @param config Model configuration list passed to the backend builder.
#'   Must contain \code{beta}, \code{gamma}, and \code{theta} (vectors of starting
#'   values), and may include \code{verbose}, \code{num_threads}, and other backend-
#'   specific options. The \code{config$package} field is ignored here — use the
#'   \code{package} argument instead.
#' @param ad_fun Combined raw handle from \code{\link{combine}}. Required when
#'   \code{package = "adlaplace"}.
#' @param package Character scalar naming the backend package to use for
#'   \code{getAdFun_r()}. Defaults to \code{"adlaplace"}. Supported backends must
#'   export a \code{getAdFun_r(data, config)} function.
#'
#' @return Backend object returned by \code{<package>::getAdFun_r()}. For the
#'   default \pkg{adlaplace} backend, this is a list containing \code{ad_fun}
#'   (external pointer handle), \code{sparsity}, and \code{hessians}. The object
#'   carries an \code{"adlaplace.backend"} attribute indicating which package built it,
#'   enabling runtime validation in functions like \code{log_lik_laplace()}.
#'
#' @seealso
#' \code{\link{inner_opt}}, \code{\link{log_lik_laplace}}, \code{\link{hessian_map}},
#' \code{\link{combine}}, \code{\link{get_ad_fun_raw}}.
#' Backend packages should export \code{getAdFun_r()}.
#'
#' @export
get_ad_fun <- function(data, config, ad_fun,
                       package = c(config$package, "adlaplace")[1]) {
  if (!is.character(package) || length(package) < 1 || is.na(package[[1]]) || !nzchar(package[[1]])) {
    stop("`package` must be a non-empty character scalar")
  }
  package <- package[[1]]

  if (identical(package, "adlaplace")) {
    if (missing(ad_fun) || is.null(ad_fun)) {
      stop("for package = 'adlaplace', supply `ad_fun` as a combined handle from combine()")
    }
    if (!is_adlaplace_handle(ad_fun)) {
      stop("`ad_fun` must be an adlaplace_handle_ptr from combine()")
    }
    return(attach_ad_fun_artifacts(ad_fun, config))
  }

  if (!missing(ad_fun)) {
    warning("`ad_fun` is ignored when package is not 'adlaplace'")
  }

  builder <- tryCatch(
    getExportedValue(package, "getAdFun_r"),
    error = function(e) NULL
  )
  if (is.null(builder)) {
    stop(
      "Backend package '", package,
      "' does not export `getAdFun_r(data, config)`. ",
      "Use adlaplace::get_ad_fun(..., package='adlaplace') or update backend package exports."
    )
  }

  out <- builder(data, config)
  if (isTRUE(config$verbose)) {
    cat("got AFun\n")
  }
  attr(out, "adlaplace.backend") <- package
  out
}
