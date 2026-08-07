#' SUN(3,3) density with joint hyperspherical parameterization
#'
#' Like [dsun()], but the length-21 parameter vector uses the joint
#' hyperspherical map of [make_sun_hs_params()] instead of the block-Cholesky
#' layout of [make_sun_params()].
#'
#' @inheritParams dsun
#' @param par Numeric vector of length 21; see [make_sun_hs_params()].
#'
#' @return A list with \code{value} and \code{error}. If \code{deriv >= 1},
#'   also \code{gradient}. If \code{deriv >= 2}, also \code{hessian}.
#'
#' @seealso [make_sun_hs_params()], [dsun_hs_fun()], [dsun()]
#' @export
dsun_hs <- function(x,
                    par,
                    log = TRUE,
                    deriv = 0L,
                    n_points = 1021L,
                    n_shifts = 8L,
                    seed = 1L,
                    n_threads = 1L,
                    weights = NULL) {
  par <- as.numeric(par)
  if (length(par) != 21L) {
    stop("par must have length 21; see make_sun_hs_params()")
  }
  deriv <- as.integer(deriv)
  if (!deriv %in% 0L:2L) {
    stop("deriv must be 0, 1, or 2")
  }
  x <- as.matrix(x)
  if (ncol(x) != 3L) {
    stop("x must be a vector of length 3 or a matrix with 3 columns")
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
      out$gradient <- numeric(21L)
    }
    if (deriv >= 2L) {
      out$hessian <- matrix(0, 21L, 21L)
    }
    return(out)
  }
  storage.mode(x) <- "double"
  dsun_hs_cpp(
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

#' Create a reusable SUN(3,3) hyperspherical log-likelihood tape
#'
#' Like [dsun_fun()], for the joint hyperspherical parameterization of
#' [make_sun_hs_params()].
#'
#' @inheritParams dsun_fun
#' @param par_seed Seed parameter vector (length 21) for [make_sun_hs_params()].
#'
#' @return An object of class \code{"admvn_sun_tape"} with \code{eval()} and
#'   \code{optim_fns()} (same interface as [dsun_fun()]).
#' @seealso [dsun_hs()], [sun_mle()]
#' @export
dsun_hs_fun <- function(x,
                        par_seed,
                        n_points = 1021L,
                        n_shifts = 8L,
                        seed = 1L,
                        n_threads = NULL,
                        weights = NULL) {
  par_seed <- as.numeric(par_seed)
  if (length(par_seed) != 21L) {
    stop("par_seed must have length 21")
  }
  x <- as.matrix(x)
  if (ncol(x) != 3L || nrow(x) < 1L) {
    stop("x must be a vector of length 3 or a matrix with 3 columns")
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
  ptr <- dsun_hs_fun_create_cpp(
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
    class = c("admvn_sun_tape", "list")
  )
}
