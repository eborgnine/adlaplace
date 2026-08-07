.dsun_square_hs <- function(x, par, d, cpp, log, deriv, n_points,
                            n_shifts, seed, n_threads, weights) {
  par <- as.numeric(par)
  p <- d * 2L + choose(2L * d, 2L)
  if (length(par) != p) stop("par has the wrong length")
  deriv <- as.integer(deriv)
  if (!deriv %in% 0L:2L) stop("deriv must be 0, 1, or 2")
  x <- as.matrix(x)
  if (ncol(x) != d) stop("x has the wrong number of columns")
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) stop("weights must have length nrow(x)")
  }
  if (!nrow(x)) {
    out <- list(value = if (log) 0 else 1, error = 0)
    if (deriv >= 1L) out$gradient <- numeric(p)
    if (deriv >= 2L) out$hessian <- matrix(0, p, p)
    return(out)
  }
  storage.mode(x) <- "double"
  cpp(x, par, isTRUE(log), deriv, as.integer(n_points), as.integer(n_shifts),
      as.integer(seed), as.integer(n_threads), weights)
}

.dsun_square_hs_fun <- function(x, par_seed, d, create_cpp, eval_cpp,
                                n_points, n_shifts, seed, n_threads, weights,
                                classes) {
  par_seed <- as.numeric(par_seed)
  p <- d * 2L + choose(2L * d, 2L)
  if (length(par_seed) != p) stop("par_seed has the wrong length")
  x <- as.matrix(x)
  if (ncol(x) != d || nrow(x) < 1L) stop("x has invalid dimensions")
  storage.mode(x) <- "double"
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != nrow(x)) stop("weights must have length nrow(x)")
  }
  if (is.null(n_threads)) n_threads <- dsun_n_threads_default_cpp()
  ptr <- create_cpp(x, par_seed, as.integer(n_points), as.integer(n_shifts),
                    as.integer(seed), as.integer(n_threads), weights)
  eval_fn <- function(par = NULL, log = TRUE, deriv = 0L, n_threads = 0L) {
    deriv <- as.integer(deriv)
    if (!deriv %in% 0L:2L) stop("deriv must be 0, 1, or 2")
    eval_cpp(ptr, par, isTRUE(log), deriv, as.integer(n_threads))
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
  structure(list(ptr = ptr, n_obs = nrow(x), n_threads = as.integer(n_threads),
                 weights = weights, eval = eval_fn, optim_fns = optim_fns),
            class = classes)
}

#' SUN(2,2) density with joint hyperspherical parameterization
#' @inheritParams dsun22
#' @param par Numeric vector of length 10; see [make_sun22_hs_params()].
#' @export
dsun22_hs <- function(x, par, log = TRUE, deriv = 0L, n_points = 1021L,
                      n_shifts = 8L, seed = 1L, n_threads = 1L,
                      weights = NULL) {
  .dsun_square_hs(x, par, 2L, dsun22_hs_cpp, log, deriv, n_points, n_shifts,
                  seed, n_threads, weights)
}

#' Create a reusable SUN(2,2) hyperspherical likelihood tape
#' @inheritParams dsun22_fun
#' @param par_seed Seed vector for [make_sun22_hs_params()].
#' @export
dsun22_hs_fun <- function(x, par_seed, n_points = 1021L, n_shifts = 8L,
                          seed = 1L, n_threads = NULL, weights = NULL) {
  .dsun_square_hs_fun(
    x, par_seed, 2L, dsun22_hs_fun_create_cpp, dsun22_fun_eval_cpp,
    n_points, n_shifts, seed, n_threads, weights,
    c("admvn_sun22_tape", "list"))
}
