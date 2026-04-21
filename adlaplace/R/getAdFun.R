#' Build backend AD function handle
#'
#' Safe frontend for constructing AD handles. This is the only public
#' \code{getAdFun()} entry point; backend packages should expose
#' \code{getAdFun_r()}.
#'
#' @param data Model data list passed to the backend builder. Must contain at
#'   minimum components required by the selected backend (e.g., \code{y}, \code{X},
#'   \code{Z} for linear mixed models). See backend package documentation.
#' @param config Model configuration list passed to the backend builder.
#'   Must contain \code{beta}, \code{gamma}, and \code{theta} (vectors of starting
#'   values), and may include \code{verbose}, \code{num_threads}, and other backend-
#'   specific options. The \code{config$package} field is ignored here — use the
#'   \code{package} argument instead.
#' @param package Character scalar naming the backend package to use for
#'   \code{getAdFun_r()}. Defaults to \code{"adlaplace"}. Supported backends must
#'   export a \code{getAdFun_r(data, config)} function.
#'
#' @return Backend object returned by \code{<package>::getAdFun_r()}. For the
#'   default \pkg{adlaplace} backend, this is a list containing \code{adFun}
#'   (external pointer handle), \code{sparsity}, and \code{hessians}. The object
#'   carries an \code{"adlaplace.backend"} attribute indicating which package built it,
#'   enabling runtime validation in functions like \code{logLikLaplace()}.
#'
#' @examples
#' \dontrun{
#' # Minimal example — requires valid data and config
#' data <- list(y = rnorm(10), X = matrix(1, 10, 1), Z = matrix(1, 10, 1))
#' config <- list(beta = 0, gamma = rep(0, 1), theta = 1, verbose = FALSE)
#' adFun <- getAdFun(data, config, package = "adlaplace")
#' }
#'
#' @seealso
#' \code{\link[adlaplace]{inner_opt}}, \code{\link[adlaplace]{logLikLaplace}},
#' \code{\link[adlaplace]{getAdFun_r}}
#'
#' @export
getAdFun <- function(data, config, package = c(config$package, "adlaplace")[1]) {
  if (!is.character(package) || length(package) < 1 || is.na(package[[1]]) || !nzchar(package[[1]])) {
    stop("`package` must be a non-empty character scalar")
  }
  package <- package[[1]]

  if (identical(package, "adlaplace")) {
    out <- getAdFun_r(data, config)
    attr(out, "adlaplace.backend") <- package
    return(out)
  }

  builder <- tryCatch(
    getExportedValue(package, "getAdFun_r"),
    error = function(e) NULL
  )
  if (is.null(builder)) {
    stop(
      "Backend package '", package,
      "' does not export `getAdFun_r(data, config)`. ",
      "Use adlaplace::getAdFun(..., package='adlaplace') or update backend package exports."
    )
  }

  out <- builder(data, config)
  if(config$verbose) {
    cat("got AFun\n")
  }
  attr(out, "adlaplace.backend") <- package
  out
}
