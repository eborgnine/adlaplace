#' Build backend AD function handle
#'
#' Safe frontend for constructing AD handles. This is the only public
#' \code{getAdFun()} entry point; backend packages should expose
#' \code{getAdFun_r()}.
#'
#' @param data Model data list passed to the backend builder.
#' @param config Model configuration list passed to the backend builder.
#' @param package Backend package name. Defaults to
#'   \code{c(config$package, "adlaplace")[1]}.
#'
#' @return Backend AD handle returned by \code{<package>::getAdFun_r()}.
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
