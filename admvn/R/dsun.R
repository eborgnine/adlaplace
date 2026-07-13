#' SUN(3,3) density with derivatives
#'
#' Evaluates the log-density (or density) of a SUN(3,3) distribution using the
#' block-Cholesky parameterization from \code{make_sun_params()}. For a matrix
#' of observations the returned value is the summed log-likelihood. Derivatives
#' are with respect to the 21-component parameter vector only; the data are
#' treated as fixed.
#'
#' @param x Numeric vector of length 3 or an \eqn{n \times 3} matrix of
#'   observations.
#' @param par Numeric vector of length 21; see [make_sun_params()].
#' @param log If \code{TRUE} (default), return the log-density / log-likelihood.
#' @param deriv Integer derivative order: \code{0} (default) value only,
#'   \code{1} also gradient, \code{2} also Hessian, all with respect to
#'   \code{par}.
#' @param n_points,n_shifts,seed Quasi-Monte Carlo settings for the orthant
#'   probabilities in the SUN density.
#' @param n_threads Number of OpenMP threads for parallel shard evaluation.
#'   Defaults to \code{1}. Use \code{0} in \code{$eval()} to reuse the value
#'   stored when the tape was built.
#'
#' @return A list with \code{value} and \code{error}. If \code{deriv >= 1},
#'   also \code{gradient} (length 21). If \code{deriv >= 2}, also
#'   \code{hessian} (\eqn{21 \times 21}).
#'
#' @seealso [make_sun_params()], [dsun_fun()]
#' @export
dsun <- function(x,
                 par,
                 log = TRUE,
                 deriv = 0L,
                 n_points = 1021L,
                 n_shifts = 8L,
                 seed = 1L,
                 n_threads = 1L) {
  par <- as.numeric(par)
  if (length(par) != 21L) {
    stop("par must have length 21; see make_sun_params()")
  }
  deriv <- as.integer(deriv)
  if (!deriv %in% 0L:2L) {
    stop("deriv must be 0, 1, or 2")
  }
  x <- as.matrix(x)
  if (ncol(x) != 3L) {
    stop("x must be a vector of length 3 or a matrix with 3 columns")
  }
  if (nrow(x) == 0L) {
    out <- list(value = if (log) 0 else 1, error = 0)
    if (deriv >= 1L) {
      out$gradient <- numeric(21L)
    }
    if (deriv >= 2L) {
      out$hessian <- matrix(0, 21L, 21L)
    }
    return(out)
  }
  storage.mode(x) <- "double"
  dsun_cpp(
    x = x,
    par = par,
    log_scale = isTRUE(log),
    deriv = deriv,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed),
    n_threads = as.integer(n_threads)
  )
}

#' Create a reusable SUN(3,3) log-likelihood tape
#'
#' Builds per-observation CppAD shard tapes for the summed SUN log-likelihood
#' with the data \code{x} fixed and the independent variables equal to the
#' 21-parameter vector. Use this for MLE: build once, then call
#' \code{$eval(par)} at new parameter values.
#'
#' @param x Numeric vector of length 3 or \eqn{n \times 3} matrix of
#'   observations fixed in the tape.
#' @param par_seed Seed parameter vector used when building the tape.
#' @param n_points,n_shifts,seed QMC settings passed to [dsun()].
#' @param n_threads Default number of OpenMP threads for \code{$eval()}. When
#'   \code{NULL}, uses the maximum available OpenMP thread count.
#'
#' @return An object of class \code{"admvn_sun_tape"} with components
#'   \code{eval(par, log = TRUE, deriv = 0, n_threads = 0)} and
#'   \code{optim_fns()} for cached fn/gr used by \code{optim()}.
#' @export
dsun_fun <- function(x,
                     par_seed,
                     n_points = 1021L,
                     n_shifts = 8L,
                     seed = 1L,
                     n_threads = NULL) {
  par_seed <- as.numeric(par_seed)
  if (length(par_seed) != 21L) {
    stop("par_seed must have length 21")
  }
  x <- as.matrix(x)
  if (ncol(x) != 3L || nrow(x) < 1L) {
    stop("x must be a vector of length 3 or a matrix with 3 columns")
  }
  storage.mode(x) <- "double"
  if (is.null(n_threads)) {
    n_threads <- dsun_n_threads_default_cpp()
  }
  ptr <- dsun_fun_create_cpp(
    x = x,
    par_seed = par_seed,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed),
    n_threads = as.integer(n_threads)
  )

  same_par <- function(a, b) {
    identical(a, b) || (length(a) == length(b) && all(a == b))
  }

  eval_fn <- function(par = NULL, log = TRUE, deriv = 0L, n_threads = 0L) {
    deriv <- as.integer(deriv)
    if (!deriv %in% 0L:2L) {
      stop("deriv must be 0, 1, or 2")
    }
    dsun_fun_eval_cpp(
      ptr,
      par = par,
      log_scale = isTRUE(log),
      deriv = deriv,
      n_threads = as.integer(n_threads)
    )
  }

  optim_fns <- function(log = TRUE, n_threads = 0L) {
    cache <- new.env(parent = emptyenv())
    cache$par <- NULL
    cache$value <- NULL
    cache$gradient <- NULL

    ensure <- function(p) {
      if (!is.null(cache$par) && same_par(p, cache$par)) {
        return(invisible(NULL))
      }
      out <- eval_fn(p, log = log, deriv = 1L, n_threads = n_threads)
      cache$par <- p
      cache$value <- out$value
      cache$gradient <- out$gradient
      invisible(NULL)
    }

    list(
      fn = function(p) {
        ensure(p)
        cache$value
      },
      gr = function(p) {
        ensure(p)
        cache$gradient
      },
      clear = function() {
        cache$par <- NULL
        cache$value <- NULL
        cache$gradient <- NULL
        invisible(NULL)
      }
    )
  }

  structure(
    list(
      ptr = ptr,
      n_obs = nrow(x),
      n_threads = as.integer(n_threads),
      eval = eval_fn,
      optim_fns = optim_fns
    ),
    class = c("admvn_sun_tape", "list")
  )
}
