test_that("SUN(3,2) hyperspherical map has expected structure", {
  par <- make_sun32_hs_start_from_normal(c(0, 0, 0), diag(3), 0.3)
  dp <- make_sun32_hs_params(par)
  expect_length(par, 16L)
  expect_equal(dim(dp$Omega), c(3L, 3L))
  expect_equal(dim(dp$Delta), c(3L, 2L))
  expect_equal(dp$tau, c(0, 0))
  expect_equal(diag(dp$Gamma), rep(1, 2), tolerance = 1e-10)
  expect_equal(dp$Delta[1:2, 1:2], diag(0.3, 2), tolerance = 1e-10)
  expect_equal(length(sun32_hs_bounds()$lower), 16L)
})

test_that("zero-z SUN(3,2) hyperspherical density is Gaussian", {
  skip_if_not_installed("mvtnorm")
  par <- c(0, 0, 0, 1, 1, 1, rep(0, 10))
  x <- matrix(c(-0.5, 0.2, 0.1, 1, -1, 0.3), 2, 3, byrow = TRUE)
  got <- dsun32_hs(x, par, n_points = 128L, n_shifts = 2L)$value
  ref <- sum(mvtnorm::dmvnorm(x, sigma = diag(3), log = TRUE))
  expect_equal(got, ref, tolerance = 1e-8)
})

test_that("dsun32_hs log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- make_sun32_hs_start_from_normal(c(0, 0, 0), diag(3), 0.25)
  dp <- make_sun32_hs_params(par)
  set.seed(32)
  x <- sn::rsun(4, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  adv <- dsun32_hs(x, par, n_points = 256L, n_shifts = 4L, n_threads = 1L)
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
})
