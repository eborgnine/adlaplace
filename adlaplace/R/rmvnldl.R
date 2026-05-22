#' @noRd
as_ldl_list <- function(chol_prec) {
  if (inherits(chol_prec, "CHMfactor") || inherits(chol_prec, "Cholesky")) {
    ex <- Matrix::expand2(chol_prec)
    return(list(
      L1 = ex$L1,
      D = as.numeric(ex$D@x),
      perm = as.integer(chol_prec@perm),
      perm_inv = as.integer(order(chol_prec@perm) - 1L)
    ))
  }
  if (is.list(chol_prec) &&
    all(c("L1", "D", "perm", "perm_inv") %in% names(chol_prec))
  ) {
    out <- list(
      L1 = chol_prec$L1,
      D = as.numeric(chol_prec$D),
      perm = as.integer(chol_prec$perm),
      perm_inv = as.integer(chol_prec$perm_inv)
    )
    return(out)
  }
  stop(
    "chol_prec must be a CHMfactor, Cholesky (LDL), or a list with L1, D, perm and perm_inv ",
    "(as in log_lik_laplace()$hessian$chol_inner)"
  )
}

#' @noRd
half_H_inv_from_ldl <- function(ldl) {
  p <- length(ldl$D)
  Linv <- Matrix::solve(ldl$L1)
  halfDinv <- Matrix::Diagonal(p, ldl$D^(-0.5))

  half_H_inv_perm <- Matrix::crossprod(Linv, halfDinv)
  half_H_inv <- half_H_inv_perm[ldl$perm_inv + 1L, ]

  # H_inv = tcrossprod(half_H_inv)

  half_H_inv
}

#' Simulate from a multivariate normal using LDL of the precision matrix
#'
#' Draws \eqn{n} samples from \eqn{N(\mu, H^{-1})} where \eqn{H} is the
#' precision matrix, given its LDL factorization (unit lower \code{L1},
#' diagonal \code{D}, and permutation \code{perm}).
#'
#' @param n Number of samples.
#' @param mean Mean vector \eqn{\mu} (length \code{length(D)}).
#' @param chol_prec Either a \code{CHMfactor}, \code{Cholesky}, or \code{dCHMsimpl}
#'   from \code{Matrix::Cholesky(..., LDL = TRUE)}, or a list like
#'   \code{res$hessian$chol_inner} with components \code{L1}, \code{D},
#'   and \code{perm}.
#'
#' @return An \code{n}-by-\code{p} matrix; each row is one draw.
#'
#' @details
#' Uses \eqn{H^{-1/2} = P^\top L^{-1} D^{-1/2}} (same as \code{reformat_chol()}
#' in \code{log_lik_deriv}), and \code{rnorm} only:
#' \eqn{x = \mu + H^{-1/2} z}, \eqn{z \sim N(0, I)}.
#'
#' @seealso \code{\link[Matrix]{Cholesky}}, \code{\link{log_lik_laplace}}
#' @export
rmvnldl <- function(n, mean = 0, chol_prec) {
  ldl <- as_ldl_list(chol_prec)
  half_H_inv <- half_H_inv_from_ldl(ldl)

  p <- length(ldl$D)
  n <- as.integer(n)
  mu <- rep_len(as.numeric(mean), p)

  z <- matrix(stats::rnorm(p * n), nrow = p, ncol = n)
  draws <- as.matrix(mu + half_H_inv %*% z)

  if (n == 1L) {
    as.vector(draws)
  } else {
    t(draws)
  }
}
