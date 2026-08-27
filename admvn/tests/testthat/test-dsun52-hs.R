test_that("SUN(5,2) hyperspherical map has expected structure", {
  par <- make_sun52_hs_start_from_normal(rep(0, 5), diag(5), 0.3)
  dp <- make_sun52_hs_params(par)
  expect_length(par, 31L)
  expect_equal(dim(dp$Omega), c(5L, 5L))
  expect_equal(dim(dp$Delta), c(5L, 2L))
  expect_equal(dp$tau, c(0, 0))
  expect_equal(diag(dp$Gamma), rep(1, 2), tolerance = 1e-10)
  expect_equal(dp$Delta[1:2, 1:2], diag(0.3, 2), tolerance = 1e-10)
  expect_equal(length(sun52_hs_bounds()$lower), 31L)
})

test_that("zero-z SUN(5,2) hyperspherical density is Gaussian", {
  skip_if_not_installed("mvtnorm")
  par <- c(rep(0, 5), rep(1, 5), rep(0, 21))
  x <- matrix(c(-0.5, 0.2, 0.1, 0.3, -0.1, 1, -1, 0.3, -0.2, 0.4), 2, 5, byrow = TRUE)
  got <- dsun52_hs(x, par, n_points = 128L, n_shifts = 2L)$value
  ref <- sum(mvtnorm::dmvnorm(x, sigma = diag(5), log = TRUE))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("dsun52_hs_fun tape matches dsun52_hs", {
  skip_if_not_installed("sn")
  par <- make_sun52_hs_start_from_normal(
    c(0.1, -0.2, 0.05, 0, 0.15), diag(c(1.1, 0.9, 1.0, 1.2, 0.85)), 0.2)
  set.seed(52)
  x <- sn::rsun(4, dp = make_sun52_hs_params(par))
  direct <- dsun52_hs(
    x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L, deriv = 1L)
  tape <- dsun52_hs_fun(x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L)
  ev <- tape$eval(par, log = TRUE, deriv = 1L)
  expect_equal(ev$value, direct$value, tolerance = 1e-8)
  expect_equal(ev$gradient, direct$gradient, tolerance = 1e-8)
})

test_that("dsun52_hs log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- make_sun52_hs_start_from_normal(rep(0, 5), diag(5), 0.25)
  dp <- make_sun52_hs_params(par)
  set.seed(52)
  x <- sn::rsun(4, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  adv <- dsun52_hs(x, par, n_points = 256L, n_shifts = 4L, n_threads = 1L)
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
})
