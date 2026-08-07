test_that("SUN(2,2) hyperspherical map has expected structure", {
  par <- make_sun22_hs_start_from_normal(c(0, 0), diag(2), 0.3)
  dp <- make_sun22_hs_params(par)
  expect_length(par, 10L)
  expect_equal(dim(dp$Omega), c(2L, 2L))
  expect_equal(diag(dp$Gamma), rep(1, 2), tolerance = 1e-10)
  expect_equal(dp$Delta, diag(0.3, 2), tolerance = 1e-10)
  expect_equal(length(sun22_hs_bounds()$lower), 10L)
})

test_that("zero-z SUN(2,2) hyperspherical density is Gaussian", {
  skip_if_not_installed("mvtnorm")
  par <- c(0, 0, 1, 1, rep(0, 6))
  x <- matrix(c(-0.5, 0.2, 1, -1), 2, 2, byrow = TRUE)
  got <- dsun22_hs(x, par, n_points = 128L, n_shifts = 2L)$value
  ref <- sum(mvtnorm::dmvnorm(x, sigma = diag(2), log = TRUE))
  expect_equal(got, ref, tolerance = 1e-8)
})
