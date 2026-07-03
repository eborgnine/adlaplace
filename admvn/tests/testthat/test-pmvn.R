test_that("pmvn matches mvtnorm for values", {
  skip_if_not_installed("mvtnorm")
  set.seed(1)
  for (n in 1:5) {
    sigma <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n)
    mean <- rnorm(n)
    lower <- rep(-Inf, n)
    upper <- rnorm(n)
    a <- pmvn(
      upper,
      lower = lower,
      mean = mean,
      sigma = sigma,
      n_points = 3000L,
      n_shifts = 12L
    )
    b <- mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean = mean,
      sigma = sigma,
      algorithm = mvtnorm::GenzBretz(maxpts = 25000, releps = 1e-4)
    )
    expect_equal(a$value, as.numeric(b[1]), tolerance = 0.05)
  }
})

test_that("pmvn gradient matches conditional formula", {
  sigma <- matrix(c(1, 0.4, 0.4, 1.5), 2, 2)
  mean <- c(0.1, -0.2)
  lower <- c(-Inf, -Inf)
  upper <- c(0.5, 0.3)

  out <- pmvn(
    upper,
    lower = lower,
    mean = mean,
    sigma = sigma,
    n_points = 4000L,
    n_shifts = 16L
  )

  beta <- sigma[2, 1] / sigma[1, 1]
  cond_mean <- mean[2] + beta * (upper[1] - mean[1])
  cond_sd <- sqrt(sigma[2, 2] - sigma[2, 1]^2 / sigma[1, 1])
  g1 <- dnorm(upper[1], mean = mean[1], sd = sqrt(sigma[1, 1])) *
    pnorm(upper[2], mean = cond_mean, sd = cond_sd)
  expect_equal(out$gradient[1], g1, tolerance = 0.02)
})

test_that("pmvn gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  sigma <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  upper <- c(0.4, -0.1)
  f <- function(x) {
    pmvn(x, lower = c(-Inf, -Inf), sigma = sigma, n_points = 2000L, n_shifts = 12L)$value
  }
  out <- pmvn(upper, lower = c(-Inf, -Inf), sigma = sigma, n_points = 2000L, n_shifts = 12L)
  g_num <- numDeriv::grad(f, upper)
  expect_equal(out$gradient, as.numeric(g_num), tolerance = 0.03)
})

test_that("pmvn hessian is symmetric and close to numDeriv", {
  skip_if_not_installed("numDeriv")
  sigma <- matrix(c(1, 0.1, 0.1, 1), 2, 2)
  upper <- c(0.2, 0.1)
  f <- function(x) {
    pmvn(x, lower = c(-Inf, -Inf), sigma = sigma, n_points = 1500L, n_shifts = 10L)$value
  }
  out <- pmvn(upper, lower = c(-Inf, -Inf), sigma = sigma, n_points = 1500L, n_shifts = 10L)
  h_num <- numDeriv::hessian(f, upper)
  expect_equal(out$hessian, unname(h_num), tolerance = 0.15)
  expect_equal(out$hessian, t(out$hessian), tolerance = 1e-10)
})

test_that("pmvn_fun matches pmvn", {
  sigma <- matrix(c(1, 0.3, 0.3, 1.2), 2, 2)
  upper <- c(0.6, 0.2)
  a <- pmvn(upper, lower = c(-Inf, -Inf), sigma = sigma, seed = 42L)
  f <- pmvn_fun(lower = c(-Inf, -Inf), sigma = sigma, seed = 42L)
  b <- f$eval(upper)
  expect_equal(b$value, a$value, tolerance = 1e-10)
  expect_equal(b$gradient, a$gradient, tolerance = 1e-10)
  expect_equal(b$hessian, a$hessian, tolerance = 1e-10)
})

test_that("univariate pmvn matches pnorm", {
  sd <- sqrt(2.5)
  out <- pmvn(0.5, lower = -Inf, sigma = matrix(2.5, 1, 1))
  expect_equal(out$value, pnorm(0.5, sd = sd))
  expect_equal(out$gradient, dnorm(0.5 / sd) / sd)
})
