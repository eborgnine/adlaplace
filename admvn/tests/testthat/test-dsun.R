test_that("make_sun_params matches expected structure", {
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  expect_length(dp$xi, 3L)
  expect_equal(dim(dp$Omega), c(3L, 3L))
  expect_equal(dim(dp$Delta), c(3L, 3L))
  expect_equal(dp$tau, c(0, 0, 0))
  expect_equal(dim(dp$Gamma), c(3L, 3L))
})

test_that("dsun log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  set.seed(1)
  x <- sn::rsun(8, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  adv <- admvn::dsun(x, par, log = TRUE, deriv = 0L, n_points = 2500L, n_shifts = 12L)
  expect_equal(adv$value, sn_val, tolerance = 0.05)
  expect_null(adv$gradient)
})

test_that("dsun gradient matches numDeriv on summed loglik", {
  skip_if_not_installed("sn")
  skip_if_not_installed("numDeriv")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  set.seed(2)
  x <- sn::rsun(4, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1500L, n_shifts = 10L, n_threads = 1L)
  f <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
  g_num <- numDeriv::grad(f, par, method = "Richardson")
  g_ad <- tape$eval(par, log = TRUE, deriv = 1L)$gradient
  expect_equal(g_ad, as.numeric(g_num), tolerance = 0.35)
})

test_that("dsun_fun reuses tape with fixed data", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  x <- sn::rsun(3, dp = dp)
  f <- dsun_fun(x, par, n_points = 1500L, n_shifts = 10L, seed = 7L, n_threads = 1L)
  direct <- admvn::dsun(
    x, par, log = TRUE, deriv = 1L,
    n_points = 1500L, n_shifts = 10L, seed = 7L, n_threads = 1L
  )
  via_tape <- f$eval(par, log = TRUE, deriv = 1L)
  expect_equal(via_tape$value, direct$value, tolerance = 1e-10)
  expect_equal(via_tape$gradient, direct$gradient, tolerance = 1e-8)
})

test_that("serial and parallel shard eval agree", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  set.seed(3)
  x <- sn::rsun(6, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1200L, n_shifts = 8L, seed = 11L, n_threads = 1L)
  serial <- tape$eval(par, log = TRUE, deriv = 1L, n_threads = 1L)
  parallel <- tape$eval(par, log = TRUE, deriv = 1L, n_threads = 4L)
  expect_equal(parallel$value, serial$value, tolerance = 1e-10)
  expect_equal(parallel$gradient, serial$gradient, tolerance = 1e-8)
})

test_that("optim_fns caches fn and gr for the same par", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
    L11 = 1, L12 = 0.5, L13 = 1,
    L21 = 0.5, L22 = 1, L23 = 1.5,
    L31 = 0, L32 = 0.5, L33 = 1,
    a = 0.3, b = 0.2, c = 0.3
  )
  dp <- make_sun_params(par)
  x <- sn::rsun(3, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1200L, n_shifts = 8L, n_threads = 1L)
  obj <- tape$optim_fns()
  direct <- tape$eval(par, log = TRUE, deriv = 1L)
  expect_equal(obj$fn(par), direct$value, tolerance = 1e-10)
  expect_equal(obj$gr(par), direct$gradient, tolerance = 1e-8)
  expect_equal(obj$fn(par), direct$value, tolerance = 1e-10)
  expect_equal(obj$gr(par), direct$gradient, tolerance = 1e-8)
})
