test_that("evaluate_sun33_slice matches weighted sn::dsun", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    omega12 = 0.2, omega13 = 0.1, omega23 = 0.2,
    L11 = 0.4, L12 = 0.1, L13 = 0.1,
    L21 = 0.1, L22 = 0.4, L23 = 0.1,
    L31 = 0, L32 = 0.1, L33 = 0.35,
    gamma12 = 0.3, gamma13 = 0.2, gamma23 = 0.3
  )
  dp <- make_sun_params(par)
  set.seed(21)
  x <- sn::rsun(6, dp = dp)
  w_raw <- runif(nrow(x), 0.3, 1.7)
  w <- w_raw / sum(w_raw)
  logp <- rnorm(nrow(x))
  quad <- list(x = x, w = w_raw, logp = logp)

  sn_logq <- sn::dsun(x, dp = dp, log = TRUE)
  sn_neg_Elq <- -sum(w * sn_logq)
  sn_kl <- sum(w * (logp - sn_logq))

  ev <- evaluate_sun33_slice(par, quad, n_points = 64L, n_shifts = 2L)
  expect_equal(unname(ev[["neg_Elq"]]), sn_neg_Elq, tolerance = 1e-4)
  expect_equal(unname(ev[["kl"]]), sn_kl, tolerance = 1e-4)
})

test_that("fit_sun33_quad recovers near-MVN on a simple weighted sample", {
  skip_if_not_installed("sn")
  set.seed(11)
  mu <- c(0.1, -0.2, 0.05)
  Sigma <- diag(c(1.2, 0.8, 1.0))
  start <- make_sun33_start_from_normal(mu, Sigma, skew_strength = 0)
  # GH product under the same normal as a discrete target measure
  gh <- gh_points(mu, Sigma, n = 4)
  # Target logp = MVN log-density so KL-minimising SUN (zero skew) should stay near start
  logp <- as.numeric(mvtnorm::dmvnorm(gh$x, mean = mu, sigma = Sigma, log = TRUE))
  # Importance weights for the GH measure under the reference are already gh$w;
  # use them directly as the discrete mass.
  quad <- list(x = gh$x, w = gh$w, logp = logp)

  fit <- fit_sun33_quad(
    quad = quad,
    start = start,
    mc.cores = 1L,
    optim_opts = list(maxit = 40L, factr = 1e7)
  )
  expect_equal(fit$engine, "admvn")
  expect_length(fit$gradient, 21L)
  expect_true(all(is.finite(fit$gradient)))
  expect_true(is.finite(fit$value))
  expect_lte(fit$value, fit$value_start + 1e-6)
  expect_true(max(abs(fit$gradient)) < 50) # soft: scaled objective
})

test_that("print_sun_kl_fit runs", {
  skip_if_not_installed("sn")
  set.seed(12)
  start <- make_sun33_start_from_normal(rep(0, 3), diag(3), 0)
  x <- matrix(rnorm(15), 5, 3)
  quad <- list(x = x, w = rep(1, 5), logp = rep(0, 5))
  fit <- fit_sun33_quad(
    quad, start, mc.cores = 1L,
    optim_opts = list(maxit = 5L)
  )
  expect_invisible(print_sun_kl_fit(fit, label = "test"))
})

test_that("fit_sun33_hs_quad recovers near-MVN on a simple weighted sample", {
  skip_if_not_installed("sn")
  skip_if_not_installed("mvtnorm")
  set.seed(13)
  mu <- c(0.1, -0.2, 0.05)
  Sigma <- diag(c(1.2, 0.8, 1.0))
  start <- make_sun33_hs_start_from_normal(mu, Sigma, skew_strength = 0)
  gh <- gh_points(mu, Sigma, n = 4)
  logp <- as.numeric(mvtnorm::dmvnorm(gh$x, mean = mu, sigma = Sigma, log = TRUE))
  quad <- list(x = gh$x, w = gh$w, logp = logp)

  fit <- fit_sun33_hs_quad(
    quad = quad,
    start = start,
    mc.cores = 1L,
    optim_opts = list(maxit = 40L, factr = 1e7)
  )
  expect_equal(fit$engine, "admvn")
  expect_length(fit$gradient, 21L)
  expect_true(all(is.finite(fit$gradient)))
  expect_true(is.finite(fit$value))
  expect_lte(fit$value, fit$value_start + 1e-6)
})
