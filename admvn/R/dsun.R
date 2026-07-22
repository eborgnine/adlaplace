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
#'   \code{optim_fns()} for \code{optim()}: \code{fn}/\code{gr} share one
#'   \code{deriv = 1} evaluation cached by parameter value.
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

    # One deriv=1 eval per new par (analytic reverse makes this barely
    # more expensive than value-only); fn and gr share the cache.
    # Sanitize for L-BFGS-B: non-finite / huge grads from wild Omega steps
    # otherwise abort optim().
    sanitize <- function(value, gradient) {
      if (!is.finite(value)) {
        value <- if (isTRUE(log)) -1e300 else 0
      }
      gradient <- as.numeric(gradient)
      gradient[!is.finite(gradient)] <- 0
      gradient <- pmin(pmax(gradient, -1e8), 1e8)
      list(value = value, gradient = gradient)
    }

    ensure <- function(p) {
      if (!is.null(cache$par) && same_par(p, cache$par) &&
            !is.null(cache$value) && !is.null(cache$gradient)) {
        return(invisible(NULL))
      }
      out <- eval_fn(p, log = log, deriv = 1L, n_threads = n_threads)
      safe <- sanitize(out$value, out$gradient)
      cache$par <- p
      cache$value <- safe$value
      cache$gradient <- safe$gradient
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

#' SUN(3,3) MLE via C++ trust-region optimization
#'
#' Maximizes the taped SUN log-likelihood using the quasi-Newton trust-region
#' solver from \pkg{trustOptim} (SR1 by default). The `nu` components are
#' optimized on the log scale so that they remain positive. The entire
#' optimization runs in C++ with one OpenMP-parallel value+gradient evaluation
#' per iteration.
#'
#' @param tape An object from [dsun_fun()].
#' @param start Numeric start vector of length 21 (natural parameterization).
#' @param control List of trust-region controls (see \pkg{trustOptim}), e.g.
#'   \code{maxit}, \code{grad.tol}, \code{report.freq},
#'   \code{quasi.newton.method} (1 = SR1, 2 = BFGS).
#' @param n_threads Threads for shard evaluation; \code{0} uses the tape default.
#'
#' @return A list with \code{par}, \code{value} (log-likelihood),
#'   \code{gradient}, \code{iterations}, \code{status}, and related fields.
#' @export
sun_mle <- function(tape, start, control = list(), n_threads = 0L) {
  if (!inherits(tape, "admvn_sun_tape")) {
    stop("tape must be created by dsun_fun()")
  }
  start <- as.numeric(start)
  if (length(start) != 21L) {
    stop("start must have length 21")
  }
  if (any(start[4:6] <= 0)) {
    stop("nu1, nu2, nu3 in start must be positive")
  }
  sun_mle_cpp(
    tape$ptr,
    start = start,
    control = control,
    n_threads = as.integer(n_threads)
  )
}
