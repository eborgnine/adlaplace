#' SUN(4,4) density with joint hyperspherical parameterization
#' @inheritParams dsun44
#' @param par Numeric vector of length 36; see [make_sun44_hs_params()].
#' @export
dsun44_hs <- function(x, par, log = TRUE, deriv = 0L, n_points = 1021L,
                      n_shifts = 8L, seed = 1L, n_threads = 1L,
                      weights = NULL) {
  .dsun_square_hs(x, par, 4L, dsun44_hs_cpp, log, deriv, n_points, n_shifts,
                  seed, n_threads, weights)
}

#' Create a reusable SUN(4,4) hyperspherical likelihood tape
#' @inheritParams dsun44_fun
#' @param par_seed Seed vector for [make_sun44_hs_params()].
#' @export
dsun44_hs_fun <- function(x, par_seed, n_points = 1021L, n_shifts = 8L,
                          seed = 1L, n_threads = NULL, weights = NULL) {
  .dsun_square_hs_fun(
    x, par_seed, 4L, dsun44_hs_fun_create_cpp, dsun44_fun_eval_cpp,
    n_points, n_shifts, seed, n_threads, weights,
    c("admvn_sun44_tape", "list"))
}
