test_that("SUN(4,4) hyperspherical map has expected structure", {
  par <- make_sun44_hs_start_from_normal(rep(0, 4), diag(4), 0.2)
  dp <- make_sun44_hs_params(par)
  expect_length(par, 36L)
  expect_equal(dim(dp$Omega), c(4L, 4L))
  expect_equal(diag(dp$Gamma), rep(1, 4), tolerance = 1e-10)
  expect_equal(dp$Delta, diag(0.2, 4), tolerance = 1e-10)
  expect_equal(length(sun44_hs_bounds()$lower), 36L)
})

test_that("zero-z SUN(4,4) hyperspherical density is Gaussian", {
  skip_if_not_installed("mvtnorm")
  par <- c(rep(0, 4), rep(1, 4), rep(0, 28))
  x <- matrix(c(-0.5, 0.2, 0.1, 0.3), 1, 4)
  got <- dsun44_hs(
    x, par, n_points = 64L, n_shifts = 2L, n_threads = 1L)$value
  ref <- mvtnorm::dmvnorm(x, sigma = diag(4), log = TRUE)
  expect_equal(got, unname(ref), tolerance = 1e-8)
})
