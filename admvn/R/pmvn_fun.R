#' Create a reusable MVN CDF tape
#'
#' Builds a CppAD tape for the multivariate normal CDF over \code{upper},
#' \code{mean}, marginal scales, and the Genz Cholesky factor. The tape is built
#' once; subsequent evaluations pass new parameter values without re-taping.
#' Derivatives with respect to \code{upper} use sparse inner patterns (as in
#' adlaplace).
#'
#' @param lower,mean,sigma Seed integration parameters used to fix the Genz
#'   permutation and build the tape (see [pmvn()]).
#' @param n_points,n_shifts,seed QMC settings passed to [pmvn()].
#'
#' @return An object of class \code{"admvn_tape"} with components \code{perm}
#'   (frozen Genz permutation) and \code{eval()}.
#'
#' @details
#' \code{eval(upper, mean = NULL, sigma = NULL, inner = TRUE)} evaluates the
#' taped function. When \code{mean} or \code{sigma} are omitted, build-time
#' values are used. \code{inner = TRUE} (default) returns the gradient and
#' Hessian with respect to \code{upper} only; \code{inner = FALSE} also returns
#' \code{gradient_mean} and includes mean, scale, and Cholesky blocks in the
#' sparse derivatives.
#'
#' @examples
#' sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' f <- pmvn_fun(lower = c(-Inf, -Inf), mean = c(0, 0), sigma = sigma)
#' f$eval(c(0.5, 0.5))
#'
#' @export
pmvn_fun <- function(lower = -Inf,
                     mean = 0,
                     sigma,
                     n_points = 1021L,
                     n_shifts = 8L,
                     seed = 1L) {
  sigma <- as.matrix(sigma)
  n <- nrow(sigma)
  if (length(lower) == 1L) {
    lower <- rep(lower, n)
  } else {
    lower <- as.numeric(lower)
  }
  if (length(mean) == 1L) {
    mean <- rep(mean, n)
  } else {
    mean <- as.numeric(mean)
  }
  if (length(lower) != n || length(mean) != n) {
    stop("lower and mean must have length ", n, " or 1")
  }
  storage.mode(sigma) <- "double"
  chol(sigma)
  ptr <- pmvn_fun_create_cpp(
    lower = lower,
    mean = mean,
    sigma = sigma,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed)
  )
  seed_mean <- mean
  seed_sigma <- sigma
  structure(
    list(
      ptr = ptr,
      perm = pmvn_fun_perm_cpp(ptr),
      eval = function(upper, mean = NULL, sigma = NULL, inner = TRUE) {
        pmvn_fun_eval_cpp(
          ptr,
          upper = as.numeric(upper),
          mean = mean,
          sigma = sigma,
          inner = isTRUE(inner)
        )
      }
    ),
    class = c("admvn_tape", "list")
  )
}
