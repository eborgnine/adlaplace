#' Inner optimization over gamma
#'
#' Runs the inner optimization (typically over \code{gamma}) using trust-region
#' sparse CG. See \code{\link{innerOpt}} for sign conventions and return values.
#'
#' @param parameters Numeric vector of fixed outer parameters
#'   (\code{beta}, \code{theta}).
#' @param gamma Numeric starting values for inner parameters.
#' @param config Configuration list (dimensions, groups, sparsity).
#' @param control Trust-region control list (see \pkg{trustOptim}).
#' @param ad_fun Backend handle from \code{\link{getAdFun}}.
#' @param deriv Logical: if \code{TRUE}, return full outer gradient and Hessian at
#'   the inner solution; if \code{FALSE}, inner quantities only. When \code{NULL},
#'   uses \code{config$deriv} or legacy \code{config$inner_only}.
#' @param adFun Alias for \code{ad_fun}.
#'
#' @return See \code{\link{innerOpt}}.
#'
#' @rdname innerOpt
#' @export
inner_opt <- function(
    parameters,
    gamma,
    config,
    control = list(),
    ad_fun = NULL,
    deriv = NULL,
    adFun = NULL
) {
  if (!is.null(adFun)) {
    if (!is.null(ad_fun)) {
      stop("supply at most one of ad_fun and adFun")
    }
    ad_fun <- adFun
  }
  if (is.null(ad_fun)) {
    stop("inner_opt requires ad_fun (or adFun)")
  }

  .Call(
    `_adlaplace_inner_opt_cpp`,
    parameters,
    gamma,
    config,
    control,
    ad_fun,
    deriv
  )
}
