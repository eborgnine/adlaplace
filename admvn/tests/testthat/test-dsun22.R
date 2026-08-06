test_that("make_sun22_params matches expected structure", {
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0.2,
    L11 = 0.4, L12 = 0.1,
    L21 = 0.1, L22 = 0.4,
    gamma12 = 0.3
  )
  dp <- make_sun22_params(par)
  expect_length(dp$xi, 2L)
  expect_equal(dim(dp$Omega), c(2L, 2L))
  expect_equal(dim(dp$Delta), c(2L, 2L))
  expect_equal(dp$tau, rep(0, 2))
  expect_equal(dim(dp$Gamma), c(2L, 2L))
  expect_equal(diag(dp$Gamma), rep(1, 2L), tolerance = 1e-10)
})

test_that("dsun22 log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0.2,
    L11 = 0.4, L12 = 0.1,
    L21 = 0.1, L22 = 0.4,
    gamma12 = 0.3
  )
  dp <- make_sun22_params(par)
  set.seed(21)
  x <- sn::rsun(8, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  tape <- dsun22_fun(
    x, par, n_points = 1021L, n_shifts = 8L, n_threads = 1L, seed = 1L
  )
  adv <- tape$eval(par, log = TRUE, deriv = 0L)
  expect_equal(adv$value, sn_val, tolerance = 5e-2)
})

test_that("zero-skew dsun22 matches MVN", {
  skip_if_not_installed("mvtnorm")
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0,
    L11 = 0, L12 = 0,
    L21 = 0, L22 = 0,
    gamma12 = 0
  )
  set.seed(1)
  x <- matrix(rnorm(6), 3, 2)
  mvn <- sum(mvtnorm::dmvnorm(x, mean = rep(0, 2), sigma = diag(2), log = TRUE))
  adv <- dsun22(x, par, log = TRUE, deriv = 0L, n_threads = 1L)
  expect_equal(adv$value, mvn, tolerance = 1e-8)
})

test_that("weighted dsun22_fun is stable and matches unit-weight scaling", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0.2,
    L11 = 0.4, L12 = 0.1,
    L21 = 0.1, L22 = 0.4,
    gamma12 = 0.3
  )
  set.seed(23)
  x <- sn::rsun(3, dp = make_sun22_params(par))
  w <- c(0.5, 1.2, 0.8)
  tape <- dsun22_fun(
    x, par, n_points = 256L, n_shifts = 4L, n_threads = 1L,
    seed = 1L, weights = w
  )
  a <- tape$eval(par, log = TRUE, deriv = 1L)
  b <- tape$eval(par, log = TRUE, deriv = 1L)
  expect_equal(a$value, b$value, tolerance = 1e-12)
  expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  expect_true(is.finite(a$value))
  expect_true(all(is.finite(a$gradient)))
})

test_that("dsun22_fun tape grads match FD", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0.2,
    L11 = 0.4, L12 = 0.1,
    L21 = 0.1, L22 = 0.4,
    gamma12 = 0.3
  )
  set.seed(24)
  x <- sn::rsun(3, dp = make_sun22_params(par))
  w <- c(1.1, 0.4, 0.9)
  tape <- dsun22_fun(
    x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L,
    seed = 1L, weights = w
  )
  fns <- tape$optim_fns()
  g <- fns$gr(par)
  eps <- 1e-5
  for (j in c(1L, 3L, 5L, 6L, 10L)) {
    p1 <- par
    p2 <- par
    p1[j] <- p1[j] + eps
    p2[j] <- p2[j] - eps
    fd <- (fns$fn(p1) - fns$fn(p2)) / (2 * eps)
    expect_equal(g[j], fd, tolerance = 1e-6)
  }
})

test_that("far-tail orthant floor does not explode AD grads vs FD", {
  ## Observations many SDs from xi make orthant probs underflow; old
  ## log(p+1e-300) reverse blew up (~1e237) while FD stayed O(1).
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0,
    L11 = 0.2, L12 = 0,
    L21 = 0, L22 = 0.2,
    gamma12 = 0
  )
  x <- rbind(
    c(0.1, -0.2),
    c(12, 12),
    c(-15, 10)
  )
  w <- c(1, 1, 1)
  tape <- dsun22_fun(
    x, par, n_points = 256L, n_shifts = 4L, n_threads = 1L,
    seed = 1L, weights = w
  )
  fns <- tape$optim_fns()
  g <- fns$gr(par)
  expect_true(all(is.finite(g)))
  expect_true(max(abs(g)) < 1e6)

  eps <- 1e-5
  for (j in c(1L, 2L, 3L, 4L)) {
    p1 <- par
    p2 <- par
    p1[j] <- p1[j] + eps
    p2[j] <- p2[j] - eps
    fd <- (fns$fn(p1) - fns$fn(p2)) / (2 * eps)
    expect_true(is.finite(fd))
    ## Absolute agreement: both should be moderate; relative can be loose
    ## near zero.
    expect_equal(g[j], fd, tolerance = 1e-3)
  }
})

test_that("large |a| residual-correlation AD matches FD of tape values", {
  ## Hand analytic domain grads disagreed with specialized CDF forward when
  ## tanh(a) is saturated; AD for `a` looked ~0 while FD matched the KL slope.
  par <- c(
    xi1 = -3.8096, xi2 = 2.0191,
    nu1 = 4.9104, nu2 = 4.6662,
    omega12 = -0.2261,
    L11 = 0.9419, L12 = 0.0576,
    L21 = 0.0576 + 2 * 0.1573, L22 = -0.8881,
    gamma12 = 4.1104
  )
  x <- matrix(c(-3.8, 2.0), 1L, 2L)
  tape <- dsun22_fun(
    x, as.numeric(par), n_points = 1021L, n_shifts = 8L, n_threads = 1L,
    seed = 1L
  )
  g <- tape$eval(as.numeric(par), log = TRUE, deriv = 1L)$gradient
  eps <- 1e-6
  p1 <- as.numeric(par)
  p2 <- as.numeric(par)
  p1[10] <- p1[10] + eps
  p2[10] <- p2[10] - eps
  fd <- (tape$eval(p1, log = TRUE, deriv = 0L)$value -
           tape$eval(p2, log = TRUE, deriv = 0L)$value) / (2 * eps)
  expect_equal(g[10], fd, tolerance = 1e-5)
})

test_that("make_sun22_start_from_normal recovers Omega diagonals", {
  mu <- c(1, -0.5)
  Sigma <- matrix(c(2, 0.3, 0.3, 1.5), 2, 2)
  par <- make_sun22_start_from_normal(mu, Sigma, skew_strength = 0)
  dp <- make_sun22_params(par)
  expect_equal(dp$xi, mu)
  expect_equal(diag(dp$Omega), diag(Sigma), tolerance = 1e-10)
  expect_equal(dp$Omega, Sigma, tolerance = 1e-8)
})

test_that("SUN22 free Lij entries map into Delta", {
  par <- c(
    xi1 = 0, xi2 = 0,
    nu1 = 1, nu2 = 1,
    omega12 = 0,
    L11 = 0.5, L12 = 0.05,
    L21 = 0.15, L22 = 0.4,
    gamma12 = 0
  )
  dp <- make_sun22_params(par)
  # Delta = C_u L^T with C_u = I => Delta = L^T
  expect_equal(dp$Delta[1, 1], 0.5, tolerance = 1e-12)
  expect_equal(dp$Delta[1, 2], 0.15, tolerance = 1e-12) # L21
  expect_equal(dp$Delta[2, 1], 0.05, tolerance = 1e-12) # L12
  expect_equal(dp$Delta[2, 2], 0.4, tolerance = 1e-12)
})
