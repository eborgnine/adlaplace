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

test_that("pmvn_fun reuses tape when mean and sigma change", {
  skip_if_not_installed("numDeriv")
  sigma <- matrix(c(1, 0.2, 0.2, 1), 2, 2)
  sigma2 <- matrix(c(1.1, 0.15, 0.15, 0.9), 2, 2)
  upper <- c(0.4, -0.1)
  f <- pmvn_fun(
    lower = c(-Inf, -Inf),
    mean = c(0, 0),
    sigma = sigma,
    n_points = 2000L,
    n_shifts = 12L,
    seed = 7L
  )
  out1 <- f$eval(upper)
  out2 <- f$eval(upper, mean = c(0.05, -0.1), sigma = sigma2)

  g_num <- numDeriv::grad(
    function(x) pmvn(x, lower = c(-Inf, -Inf), mean = c(0.05, -0.1),
                     sigma = sigma2, n_points = 2000L, n_shifts = 12L,
                     seed = 7L)$value,
    upper
  )
  expect_equal(out2$gradient, as.numeric(g_num), tolerance = 0.04)
  expect_false(isTRUE(all.equal(out1$value, out2$value, tolerance = 1e-6)))
})

test_that("pack_genz_ch round-trip matches pmvn", {
  sigma <- matrix(c(1, 0.3, 0.3, 1.2), 2, 2)
  upper <- c(0.6, 0.2)
  f <- pmvn_fun(lower = c(-Inf, -Inf), mean = c(0, 0), sigma = sigma, seed = 42L)
  direct <- pmvn(upper, lower = c(-Inf, -Inf), mean = c(0, 0), sigma = sigma, seed = 42L)
  via_tape <- f$eval(upper, mean = c(0, 0), sigma = sigma)
  expect_equal(via_tape$value, direct$value, tolerance = 1e-10)
  pack <- pack_genz_ch(sigma, f$perm)
  expect_equal(length(pack$scale), 2L)
  expect_equal(nrow(pack$ch), 2L)
})

test_that("inner = FALSE gradient w.r.t. mean matches finite differences", {
  skip_if_not_installed("numDeriv")
  sigma <- matrix(c(1, 0.15, 0.15, 1.1), 2, 2)
  mean0 <- c(0.05, -0.1)
  upper <- c(0.35, 0.2)
  f <- pmvn_fun(
    lower = c(-Inf, -Inf),
    mean = mean0,
    sigma = sigma,
    n_points = 1500L,
    n_shifts = 10L,
    seed = 11L
  )
  out <- f$eval(upper, inner = FALSE)
  g_num <- numDeriv::grad(
    function(m) f$eval(upper, mean = m, inner = FALSE)$value,
    mean0
  )
  expect_equal(out$gradient_mean, as.numeric(g_num), tolerance = 0.05)
})

test_that("univariate pmvn matches pnorm", {
  sd <- sqrt(2.5)
  out <- pmvn(0.5, lower = -Inf, sigma = matrix(2.5, 1, 1))
  expect_equal(out$value, pnorm(0.5, sd = sd))
  expect_equal(out$gradient, dnorm(0.5 / sd) / sd)
})

test_that("analytic domain grads for upper and scale match finite differences", {
  # Domain treats (scale, chol(R)) as independent; free chol perturbations that
  # break unit-diagonal R are not the SUN outer-tape path (which re-cov2cors).
  set.seed(3)
  A <- matrix(rnorm(9), 3, 3)
  sigma0 <- crossprod(A) + diag(3)
  scale0 <- sqrt(diag(sigma0))
  L0 <- t(chol(cov2cor(sigma0)))
  upper0 <- c(0.4, -0.1, 0.2)
  mean0 <- c(0.1, 0.0, -0.05)
  g <- admvn:::pmvn_domain_grad_cpp(upper0, mean0, scale0, L0)
  n <- 3L
  f_val <- function(upper, scale) {
    R <- L0 %*% t(L0)
    sigma <- diag(scale) %*% R %*% diag(scale)
    pmvn(upper, lower = rep(-Inf, 3), mean = mean0, sigma = sigma)$value
  }
  h <- 1e-6
  for (j in 1:3) {
    uu <- ud <- upper0
    uu[j] <- upper0[j] + h
    ud[j] <- upper0[j] - h
    fd <- (f_val(uu, scale0) - f_val(ud, scale0)) / (2 * h)
    expect_equal(g[j], fd, tolerance = 5e-5)
  }
  for (j in 1:3) {
    su <- sd <- scale0
    su[j] <- scale0[j] + h
    sd[j] <- scale0[j] - h
    fd <- (f_val(upper0, su) - f_val(upper0, sd)) / (2 * h)
    expect_equal(g[2L * n + j], fd, tolerance = 5e-5)
  }
})

test_that("specialized bivariate/trivariate CDF matches mnormt::pmnorm", {
  skip_if_not_installed("mnormt")
  set.seed(11)
  # bivariate
  sigma2 <- matrix(c(1.2, 0.4, 0.4, 0.9), 2, 2)
  upper2 <- c(0.3, -0.2)
  mean2 <- c(0.1, 0.0)
  a2 <- pmvn(upper2, lower = rep(-Inf, 2), mean = mean2, sigma = sigma2)
  b2 <- mnormt::pmnorm(upper2, mean2, sigma2)
  expect_equal(a2$value, as.numeric(b2), tolerance = 5e-6)

  # trivariate
  A <- matrix(rnorm(9), 3, 3)
  sigma3 <- crossprod(A) + diag(3)
  upper3 <- rnorm(3)
  mean3 <- rnorm(3)
  a3 <- pmvn(upper3, lower = rep(-Inf, 3), mean = mean3, sigma = sigma3)
  b3 <- mnormt::pmnorm(upper3, mean3, sigma3)
  expect_equal(a3$value, as.numeric(b3), tolerance = 5e-5)
})

test_that("compiled double value matches taped Forward path", {
  # Specialized d=2 CDF is used for both; values should still agree.
  sigma <- matrix(c(1, 0.25, 0.25, 1.1), 2, 2)
  upper <- c(0.3, -0.1)
  lower <- c(-Inf, -Inf)
  mean <- c(0, 0)
  a <- pmvn(
    upper, lower = lower, mean = mean, sigma = sigma,
    n_points = 2000L, n_shifts = 10L, seed = 9L
  )
  f <- pmvn_fun(
    lower = lower, mean = mean, sigma = sigma,
    n_points = 2000L, n_shifts = 10L, seed = 9L
  )
  b <- f$eval(upper)
  expect_equal(a$value, b$value, tolerance = 1e-12)
})

test_that("analytic trivariate upper gradient matches conditional formula", {
  skip_if_not_installed("mvtnorm")
  sigma <- matrix(c(
    1.0, 0.3, 0.2,
    0.3, 1.2, 0.1,
    0.2, 0.1, 1.1
  ), 3, 3)
  upper <- c(0.4, 0.1, -0.2)
  out <- pmvn(
    upper,
    lower = rep(-Inf, 3),
    sigma = sigma,
    n_points = 1500L,
    n_shifts = 8L
  )
  s11 <- sigma[1, 1]
  beta <- as.numeric(sigma[2:3, 1] / s11)
  cond_mean <- beta * upper[1]
  cond_cov <- sigma[2:3, 2:3] - tcrossprod(sigma[2:3, 1]) / s11
  g1 <- dnorm(upper[1], sd = sqrt(s11)) *
    as.numeric(mvtnorm::pmvnorm(upper = upper[2:3] - cond_mean, sigma = cond_cov))
  expect_equal(out$gradient[1], g1, tolerance = 1e-6)
})

test_that("dsun gradient still matches numDeriv with analytic pmvn reverse", {
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
  set.seed(4)
  x <- sn::rsun(3, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1200L, n_shifts = 8L, n_threads = 1L)
  f <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
  g_num <- numDeriv::grad(f, par, method = "Richardson")
  g_ad <- tape$eval(par, log = TRUE, deriv = 1L)$gradient
  expect_equal(g_ad, as.numeric(g_num), tolerance = 0.4)
})
