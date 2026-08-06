#' SUN(4,4) density with derivatives
#'
#' Evaluates the log-density (or density) of a SUN(4,4) distribution using the
#' block-Cholesky parameterization from \code{make_sun44_params()}. For a matrix
#' of observations the returned value is the summed log-likelihood, or a
#' weighted sum if \code{weights} is supplied. Derivatives are with respect to
#' the 36-component parameter vector only.
#'
#' Orthant probabilities use Genz QMC with CppAD (dimension 4), so prefer a
#' modest \code{n_points}/\code{n_shifts} budget for speed.
#'
#' @param x Numeric vector of length 4 or an \eqn{n \times 4} matrix of
#'   observations.
#' @param par Numeric vector of length 36; see [make_sun44_params()].
#' @param log If \code{TRUE} (default), return the log-density / log-likelihood.
#' @param deriv Integer derivative order: \code{0} (default) value only,
#'   \code{1} also gradient, \code{2} also Hessian, all with respect to
#'   \code{par}.
#' @param n_points,n_shifts,seed Quasi-Monte Carlo settings for orthant CDFs.
#' @param n_threads Number of OpenMP threads for parallel shard evaluation.
#' @param weights Optional numeric vector of length \code{nrow(x)}.
#'
#' @return A list with \code{value} and \code{error}. If \code{deriv >= 1},
#'   also \code{gradient} (length 36). If \code{deriv >= 2}, also
#'   \code{hessian} (\eqn{36 \times 36}).
#'
#' @seealso [make_sun44_params()], [dsun44_fun()], [dsun()]
#' @export
dsun44 <- function(x,
                   par,
                   log = TRUE,
                   deriv = 0L,
                   n_points = 1021L,
                   n_shifts = 8L,
                   seed = 1L,
                   n_threads = 1L,
                   weights = NULL) {
  par <- as.numeric(par)
  if (length(par) != 36L) {
    stop("par must have length 36; see make_sun44_params()")
  }
  deriv <- as.integer(deriv)
  if (!deriv %in% 0L:2L) {
    stop("deriv must be 0, 1, or 2")
  }
  x <- as.matrix(x)
  if (ncol(x) != 4L) {
    stop("x must be a vector of length 4 or a matrix with 4 columns")
  }
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) {
      stop("weights must have length nrow(x)")
    }
  }
  if (nrow(x) == 0L) {
    out <- list(value = if (log) 0 else 1, error = 0)
    if (deriv >= 1L) {
      out$gradient <- numeric(36L)
    }
    if (deriv >= 2L) {
      out$hessian <- matrix(0, 36L, 36L)
    }
    return(out)
  }
  storage.mode(x) <- "double"
  dsun44_cpp(
    x = x,
    par = par,
    log_scale = isTRUE(log),
    deriv = deriv,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed),
    n_threads = as.integer(n_threads),
    weights = weights
  )
}

#' Create a reusable SUN(4,4) log-likelihood tape
#'
#' Builds per-observation CppAD shard tapes for the (weighted) SUN(4,4)
#' log-likelihood with data \code{x} fixed. Use for KL / MLE: build once, then
#' call \code{$eval(par)} or \code{$optim_fns()}.
#'
#' @inheritParams dsun44
#' @param par_seed Seed parameter vector used when building the tape.
#' @param n_threads Default OpenMP threads for \code{$eval()}. When
#'   \code{NULL}, uses the maximum available OpenMP thread count.
#'
#' @return An object of class \code{"admvn_sun44_tape"} with
#'   \code{eval()} and \code{optim_fns()} (cached \code{fn}/\code{gr}).
#' @export
dsun44_fun <- function(x,
                       par_seed,
                       n_points = 1021L,
                       n_shifts = 8L,
                       seed = 1L,
                       n_threads = NULL,
                       weights = NULL) {
  par_seed <- as.numeric(par_seed)
  if (length(par_seed) != 36L) {
    stop("par_seed must have length 36")
  }
  x <- as.matrix(x)
  if (ncol(x) != 4L || nrow(x) < 1L) {
    stop("x must be a vector of length 4 or a matrix with 4 columns")
  }
  storage.mode(x) <- "double"
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) {
      stop("weights must have length nrow(x)")
    }
  }
  if (is.null(n_threads)) {
    n_threads <- dsun_n_threads_default_cpp()
  }
  ptr <- dsun44_fun_create_cpp(
    x = x,
    par_seed = par_seed,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed),
    n_threads = as.integer(n_threads),
    weights = weights
  )

  same_par <- function(a, b) {
    identical(a, b) || (length(a) == length(b) && all(a == b))
  }

  eval_fn <- function(par = NULL, log = TRUE, deriv = 0L, n_threads = 0L) {
    deriv <- as.integer(deriv)
    if (!deriv %in% 0L:2L) {
      stop("deriv must be 0, 1, or 2")
    }
    dsun44_fun_eval_cpp(
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
      weights = weights,
      eval = eval_fn,
      optim_fns = optim_fns
    ),
    class = c("admvn_sun44_tape", "list")
  )
}
