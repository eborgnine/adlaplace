test_that("SUN(4,2) hyperspherical map has expected structure", {
  par <- make_sun42_hs_start_from_normal(c(0, 0, 0, 0), diag(4), 0.3)
  dp <- make_sun42_hs_params(par)
  expect_length(par, 23L)
  expect_equal(dim(dp$Omega), c(4L, 4L))
  expect_equal(dim(dp$Delta), c(4L, 2L))
  expect_equal(dp$tau, c(0, 0))
  expect_equal(diag(dp$Gamma), rep(1, 2), tolerance = 1e-10)
  expect_equal(dp$Delta[1:2, 1:2], diag(0.3, 2), tolerance = 1e-10)
  expect_equal(length(sun42_hs_bounds()$lower), 23L)
})

test_that("zero-z SUN(4,2) hyperspherical density is Gaussian", {
  skip_if_not_installed("mvtnorm")
  par <- c(0, 0, 0, 0, 1, 1, 1, 1, rep(0, 15))
  x <- matrix(c(-0.5, 0.2, 0.1, 0.3, 1, -1, 0.3, -0.2), 2, 4, byrow = TRUE)
  got <- dsun42_hs(x, par, n_points = 128L, n_shifts = 2L)$value
  ref <- sum(mvtnorm::dmvnorm(x, sigma = diag(4), log = TRUE))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("dsun42_hs_fun tape matches dsun42_hs", {
  skip_if_not_installed("sn")
  par <- make_sun42_hs_start_from_normal(
    c(0.1, -0.2, 0.05, 0), diag(c(1.1, 0.9, 1.0, 1.2)), 0.2)
  set.seed(42)
  x <- sn::rsun(4, dp = make_sun42_hs_params(par))
  direct <- dsun42_hs(
    x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L, deriv = 1L)
  tape <- dsun42_hs_fun(x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L)
  ev <- tape$eval(par, log = TRUE, deriv = 1L)
  expect_equal(ev$value, direct$value, tolerance = 1e-8)
  expect_equal(ev$gradient, direct$gradient, tolerance = 1e-8)
})

test_that("dsun42_hs log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- make_sun42_hs_start_from_normal(c(0, 0, 0, 0), diag(4), 0.25)
  dp <- make_sun42_hs_params(par)
  set.seed(42)
  x <- sn::rsun(4, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  adv <- dsun42_hs(x, par, n_points = 256L, n_shifts = 4L, n_threads = 1L)
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
})
