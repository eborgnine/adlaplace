#' Extract inner inverse Cholesky factor from Laplace output
#'
#' Returns \eqn{H^{-1/2}} for the inner (random-effect) Hessian from the output
#' of \code{\link{log_lik_laplace}}. Accepts the flat \code{half_H_inv} matrix
#' when \code{deriv = TRUE}, or reconstructs it from \code{hessian$chol_inner}.
#'
#' @param laplace Output of \code{log_lik_laplace(..., deriv = TRUE)} (or
#'   \code{inner_opt(..., deriv = TRUE)} with the same \code{extra$hessian} layout).
#'
#' @return A sparse or dense matrix \eqn{H^{-1/2}} with \code{nrow = Ngamma}.
#'
#' @seealso \code{\link{log_lik_laplace}}, \code{\link{rmvnldl}}
#' @export
laplace_half_H_inv <- function(laplace) {
  if (!is.list(laplace)) {
    stop("laplace must be a list (output of log_lik_laplace)", call. = FALSE)
  }

  half_h <- laplace$extra$half_H_inv
  if (is.null(half_h) && is.list(laplace$extra)) {
    half_h <- laplace$extra$extra$half_H_inv
  }
  if (!is.null(half_h)) {
    return(half_h)
  }

  chol_inner <- NULL
  if (is.list(laplace$extra) && is.list(laplace$extra$hessian)) {
    chol_inner <- laplace$extra$hessian$chol_inner
  }
  if (is.null(chol_inner)) {
    stop(
      "laplace must contain extra$half_H_inv (from log_lik_laplace(..., deriv = TRUE)) ",
      "or extra$hessian$chol_inner",
      call. = FALSE
    )
  }

  ldl <- as_ldl_list(chol_inner)
  half_H_inv_from_ldl(ldl)
}
