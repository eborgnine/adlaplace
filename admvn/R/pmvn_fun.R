#' Create a reusable MVN CDF tape
#'
#' Builds a CppAD tape for the multivariate normal CDF that depends on
#' \code{lower}, \code{mean}, and \code{sigma} but can be evaluated at many
#' different \code{upper} limits without re-taping.
#'
#' @param lower,mean,sigma Fixed integration parameters (see \code{\link{pmvn}}).
#' @param n_points,n_shifts,seed QMC settings passed to \code{\link{pmvn}}.
#'
#' @return An external pointer object of class \code{"admvn_tape_ptr"} with
#'   an \code{eval()} method.
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
  structure(
    list(
      ptr = ptr,
      eval = function(upper) {
        pmvn_fun_eval_cpp(ptr, as.numeric(upper))
      }
    ),
    class = c("admvn_tape", "list")
  )
}
