test_that("make_sun44_params matches expected structure", {
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1, nu4 = 1,
    omega12 = 0.1, omega13 = 0.05, omega23 = 0.1,
    omega14 = 0.02, omega24 = 0.03, omega34 = 0.04,
    L11 = 0.5, L12 = 0.1, L13 = 0, L14 = 0,
    L21 = 0.1, L22 = 0.5, L23 = 0.1, L24 = 0,
    L31 = 0, L32 = 0.1, L33 = 0.5, L34 = 0.1,
    L41 = 0, L42 = 0, L43 = 0.1, L44 = 0.5,
    gamma12 = 0.1, gamma13 = 0.05, gamma23 = 0.1,
    gamma14 = 0, gamma24 = 0.05, gamma34 = 0.1
  )
  dp <- make_sun44_params(par)
  expect_length(dp$xi, 4L)
  expect_equal(dim(dp$Omega), c(4L, 4L))
  expect_equal(dim(dp$Delta), c(4L, 4L))
  expect_equal(dp$tau, rep(0, 4))
  expect_equal(dim(dp$Gamma), c(4L, 4L))
  expect_equal(diag(dp$Gamma), rep(1, 4L), tolerance = 1e-10)
})

test_that("dsun44 log-density matches sn::dsun", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1, nu4 = 1,
    omega12 = 0.1, omega13 = 0.05, omega23 = 0.1,
    omega14 = 0.02, omega24 = 0.03, omega34 = 0.04,
    L11 = 0.5, L12 = 0.1, L13 = 0, L14 = 0,
    L21 = 0.1, L22 = 0.5, L23 = 0.1, L24 = 0,
    L31 = 0, L32 = 0.1, L33 = 0.5, L34 = 0.1,
    L41 = 0, L42 = 0, L43 = 0.1, L44 = 0.5,
    gamma12 = 0.1, gamma13 = 0.05, gamma23 = 0.1,
    gamma14 = 0, gamma24 = 0.05, gamma34 = 0.1
  )
  dp <- make_sun44_params(par)
  set.seed(21)
  x <- sn::rsun(6, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  # Persistent tape (rebuilding one-shot tapes is not bit-stable under QMC AD)
  tape <- dsun44_fun(
    x, par, n_points = 1021L, n_shifts = 8L, n_threads = 1L, seed = 1L
  )
  adv <- tape$eval(par, log = TRUE, deriv = 0L)
  expect_equal(adv$value, sn_val, tolerance = 5e-2)
})

test_that("zero-skew dsun44 matches MVN", {
  skip_if_not_installed("mvtnorm")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1, nu4 = 1,
    omega12 = 0, omega13 = 0, omega23 = 0,
    omega14 = 0, omega24 = 0, omega34 = 0,
    L11 = 0, L12 = 0, L13 = 0, L14 = 0,
    L21 = 0, L22 = 0, L23 = 0, L24 = 0,
    L31 = 0, L32 = 0, L33 = 0, L34 = 0,
    L41 = 0, L42 = 0, L43 = 0, L44 = 0,
    gamma12 = 0, gamma13 = 0, gamma23 = 0,
    gamma14 = 0, gamma24 = 0, gamma34 = 0
  )
  set.seed(1)
  x <- matrix(rnorm(8), 2, 4)
  mvn <- sum(mvtnorm::dmvnorm(x, mean = rep(0, 4), sigma = diag(4), log = TRUE))
  adv <- dsun44(x, par, log = TRUE, deriv = 0L, n_threads = 1L)
  expect_equal(adv$value, mvn, tolerance = 1e-8)
})

test_that("weighted dsun44_fun is stable and matches unit-weight scaling", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1, nu4 = 1,
    omega12 = 0.1, omega13 = 0.05, omega23 = 0.1,
    omega14 = 0.02, omega24 = 0.03, omega34 = 0.04,
    L11 = 0.5, L12 = 0.1, L13 = 0, L14 = 0,
    L21 = 0.1, L22 = 0.5, L23 = 0.1, L24 = 0,
    L31 = 0, L32 = 0.1, L33 = 0.5, L34 = 0.1,
    L41 = 0, L42 = 0, L43 = 0.1, L44 = 0.5,
    gamma12 = 0.1, gamma13 = 0.05, gamma23 = 0.1,
    gamma14 = 0, gamma24 = 0.05, gamma34 = 0.1
  )
  set.seed(23)
  x <- sn::rsun(3, dp = make_sun44_params(par))
  w <- c(0.5, 1.2, 0.8)
  tape <- dsun44_fun(
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

test_that("dsun44_fun tape matches repeated eval; grads match FD", {
  skip_if_not_installed("sn")
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0, xi4 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1, nu4 = 1,
    omega12 = 0.1, omega13 = 0.05, omega23 = 0.1,
    omega14 = 0.02, omega24 = 0.03, omega34 = 0.04,
    L11 = 0.5, L12 = 0.1, L13 = 0, L14 = 0,
    L21 = 0.1, L22 = 0.5, L23 = 0.1, L24 = 0,
    L31 = 0, L32 = 0.1, L33 = 0.5, L34 = 0.1,
    L41 = 0, L42 = 0, L43 = 0.1, L44 = 0.5,
    gamma12 = 0.1, gamma13 = 0.05, gamma23 = 0.1,
    gamma14 = 0, gamma24 = 0.05, gamma34 = 0.1
  )
  set.seed(24)
  x <- sn::rsun(3, dp = make_sun44_params(par))
  w <- c(1.1, 0.4, 0.9)
  tape <- dsun44_fun(
    x, par, n_points = 128L, n_shifts = 2L, n_threads = 1L,
    seed = 1L, weights = w
  )
  fns <- tape$optim_fns()
  g <- fns$gr(par)
  eps <- 1e-5
  for (j in c(1L, 5L, 15L, 30L)) {
    p1 <- par
    p2 <- par
    p1[j] <- p1[j] + eps
    p2[j] <- p2[j] - eps
    fd <- (fns$fn(p1) - fns$fn(p2)) / (2 * eps)
    expect_equal(g[j], fd, tolerance = 1e-6)
  }
})
