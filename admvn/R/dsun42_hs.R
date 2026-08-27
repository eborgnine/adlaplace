#' SUN(4,2) density with joint hyperspherical parameterization
#'
#' @param x Numeric matrix with 4 columns (observations in rows).
#' @param par Numeric vector of length 23; see [make_sun42_hs_params()].
#' @param log Logical; return log-density / log-likelihood if \code{TRUE}.
#' @param deriv Integer 0, 1, or 2 for value / gradient / Hessian.
#' @param n_points,n_shifts Orthant QMC settings.
#' @param seed QMC seed.
#' @param n_threads OpenMP threads.
#' @param weights Optional observation weights (length \code{nrow(x)}).
#' @export
dsun42_hs <- function(x, par, log = TRUE, deriv = 0L, n_points = 1021L,
                      n_shifts = 8L, seed = 1L, n_threads = 1L,
                      weights = NULL) {
  par <- as.numeric(par)
  if (length(par) != 23L) stop("par must have length 23 for SUN(4,2) hs")
  deriv <- as.integer(deriv)
  if (!deriv %in% 0L:2L) stop("deriv must be 0, 1, or 2")
  x <- as.matrix(x)
  if (ncol(x) != 4L) stop("x must have 4 columns for SUN(4,2)")
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) stop("weights must have length nrow(x)")
  }
  if (!nrow(x)) {
    out <- list(value = if (log) 0 else 1, error = 0)
    if (deriv >= 1L) out$gradient <- numeric(23L)
    if (deriv >= 2L) out$hessian <- matrix(0, 23L, 23L)
    return(out)
  }
  storage.mode(x) <- "double"
  dsun42_hs_cpp(
    x, par, isTRUE(log), deriv, as.integer(n_points), as.integer(n_shifts),
    as.integer(seed), as.integer(n_threads), weights)
}

#' Create a reusable SUN(4,2) hyperspherical likelihood tape
#'
#' @param x Numeric matrix with 4 columns.
#' @param par_seed Seed vector for [make_sun42_hs_params()].
#' @param n_points,n_shifts Orthant QMC settings.
#' @param seed QMC seed.
#' @param n_threads OpenMP threads (default from package helper).
#' @param weights Optional observation weights.
#' @export
dsun42_hs_fun <- function(x, par_seed, n_points = 1021L, n_shifts = 8L,
                          seed = 1L, n_threads = NULL, weights = NULL) {
  par_seed <- as.numeric(par_seed)
  if (length(par_seed) != 23L) stop("par_seed must have length 23 for SUN(4,2) hs")
  x <- as.matrix(x)
  if (ncol(x) != 4L || nrow(x) < 1L) stop("x has invalid dimensions")
  storage.mode(x) <- "double"
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) stop("weights must have length nrow(x)")
  }
  if (is.null(n_threads)) n_threads <- dsun_n_threads_default_cpp()
  ptr <- dsun42_hs_fun_create_cpp(
    x, par_seed, as.integer(n_points), as.integer(n_shifts),
    as.integer(seed), as.integer(n_threads), weights)
  eval_fn <- function(par = NULL, log = TRUE, deriv = 0L, n_threads = 0L) {
    deriv <- as.integer(deriv)
    if (!deriv %in% 0L:2L) stop("deriv must be 0, 1, or 2")
    dsun42_fun_eval_cpp(ptr, par, isTRUE(log), deriv, as.integer(n_threads))
  }
  optim_fns <- function(log = TRUE, n_threads = 0L) {
    cache <- new.env(parent = emptyenv())
    cache$par <- cache$value <- cache$gradient <- NULL
    ensure <- function(pv) {
      if (!is.null(cache$par) && identical(pv, cache$par)) return(invisible())
      out <- eval_fn(pv, log, 1L, n_threads)
      cache$par <- pv
      cache$value <- if (is.finite(out$value)) out$value else if (log) -1e300 else 0
      cache$gradient <- pmin(pmax(replace(
        as.numeric(out$gradient), !is.finite(out$gradient), 0), -1e8), 1e8)
      invisible()
    }
    list(
      fn = function(pv) { ensure(pv); cache$value },
      gr = function(pv) { ensure(pv); cache$gradient },
      clear = function() {
        cache$par <- cache$value <- cache$gradient <- NULL
        invisible()
      }
    )
  }
  structure(
    list(
      ptr = ptr, n_obs = nrow(x), n_threads = as.integer(n_threads),
      weights = weights, eval = eval_fn, optim_fns = optim_fns
    ),
    class = c("admvn_sun42_tape", "list")
  )
}
