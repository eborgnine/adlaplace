test_that("make_sun_hs_params has expected structure", {
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1.2, nu2 = 0.8, nu3 = 1.5,
    z21 = 0.2, z31 = 0.1, z32 = 0.15,
    z41 = 0.3, z42 = 0.05, z43 = -0.1,
    z51 = 0, z52 = 0.25, z53 = 0.1, z54 = 0.05,
    z61 = -0.05, z62 = 0, z63 = 0.2, z64 = 0.1, z65 = -0.05
  )
  dp <- make_sun_hs_params(par)
  expect_length(dp$xi, 3L)
  expect_equal(dim(dp$Omega), c(3L, 3L))
  expect_equal(dim(dp$Delta), c(3L, 3L))
  expect_equal(dp$tau, c(0, 0, 0))
  expect_equal(dim(dp$Gamma), c(3L, 3L))
  expect_equal(diag(dp$Gamma), rep(1, 3L), tolerance = 1e-10)
  Cu <- cov2cor(dp$Omega)
  expect_equal(diag(Cu), rep(1, 3L), tolerance = 1e-10)
  expect_equal(
    unname(diag(dp$Omega)),
    unname(c(par["nu1"], par["nu2"], par["nu3"])^2),
    tolerance = 1e-10
  )
})

test_that("independent-SN pattern isolates pair slots", {
  d <- c(0.4, -0.3, 0.2)
  z <- setNames(numeric(15L), admvn:::.sun_hs_z_names(3L, 3L))
  z[c("z41", "z52", "z63")] <- atanh(d)
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    z
  )
  dp <- make_sun_hs_params(par)
  expect_equal(dp$Delta, diag(d), tolerance = 1e-10)
  expect_equal(dp$Gamma, diag(3), tolerance = 1e-10)
  expect_equal(cov2cor(dp$Omega), diag(3), tolerance = 1e-10)
})

test_that("zero-z recovers MVN density", {
  skip_if_not_installed("mvtnorm")
  par <- c(
    xi1 = 0.1, xi2 = -0.2, xi3 = 0.3,
    nu1 = 1.1, nu2 = 0.9, nu3 = 1.2,
    setNames(numeric(15L), admvn:::.sun_hs_z_names(3L, 3L))
  )
  ## mild UU correlation so Omega is not diagonal
  par[c("z21", "z31", "z32")] <- c(0.2, 0.1, 0.15)
  dp <- make_sun_hs_params(par)
  expect_equal(dp$Delta, matrix(0, 3, 3), tolerance = 1e-12)
  expect_equal(dp$Gamma, diag(3), tolerance = 1e-12)

  set.seed(11)
  x <- matrix(rnorm(15), 5, 3)
  mvn <- sum(mvtnorm::dmvnorm(x, mean = dp$xi, sigma = dp$Omega, log = TRUE))
  adv <- dsun_hs(x, par, log = TRUE, n_points = 64L, n_shifts = 2L)
  expect_equal(adv$value, mvn, tolerance = 1e-6)
})

test_that("lambda_min(J) stays positive at hs bound corners", {
  b <- sun33_hs_bounds()
  ## probe a few extreme z corners (xi/nu fixed)
  base <- c(0, 0, 0, 1, 1, 1)
  corners <- list(
    b$lower[7:21],
    b$upper[7:21],
    {
      z <- b$lower[7:21]
      pair <- admvn:::.sun_hs_pair_z_idx(3L)
      z[pair] <- b$upper[6L + pair]
      z
    }
  )
  for (z in corners) {
    J <- admvn:::.corr_from_free(z, 6L, eps = 1e-6)
    expect_true(min(eigen(J, only.values = TRUE)$values) > 0)
    expect_equal(diag(J), rep(1, 6L), tolerance = 1e-8)
  }
})

test_that("dsun_hs matches sn::dsun on hs-built dp", {
  skip_if_not_installed("sn")
  z <- setNames(numeric(15L), admvn:::.sun_hs_z_names(3L, 3L))
  z[c("z21", "z31", "z32")] <- c(0.2, 0.1, 0.15)
  z[c("z41", "z52", "z63")] <- atanh(c(0.35, 0.25, 0.2))
  z[c("z42", "z53", "z54")] <- c(0.05, -0.05, 0.1)
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    z
  )
  dp <- make_sun_hs_params(par)
  set.seed(3)
  x <- sn::rsun(8, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  adv <- dsun_hs(x, par, log = TRUE, n_points = 64L, n_shifts = 2L)
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
})

test_that("dsun_hs gradient matches numDeriv", {
  skip_if_not_installed("sn")
  skip_if_not_installed("numDeriv")
  z <- setNames(numeric(15L), admvn:::.sun_hs_z_names(3L, 3L))
  z[c("z21", "z31", "z32")] <- c(0.15, 0.1, 0.12)
  z[c("z41", "z52", "z63")] <- atanh(c(0.3, 0.2, 0.25))
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    z
  )
  set.seed(4)
  x <- sn::rsun(4, dp = make_sun_hs_params(par))
  tape <- dsun_hs_fun(x, par, n_points = 1500L, n_shifts = 10L, n_threads = 1L)
  f <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
  g_num <- numDeriv::grad(f, par, method = "Richardson")
  g_ad <- tape$eval(par, log = TRUE, deriv = 1L)$gradient
  expect_equal(g_ad, as.numeric(g_num), tolerance = 0.35)
})

test_that("legacy dsun and dsun_hs agree on a shared dp", {
  skip_if_not_installed("sn")
  ## Independent SN: legacy L = diag(d), Cu = I, residual C = I (gamma = 0)
  ## hs: only pair z nonzero.
  d <- c(0.4, 0.3, 0.25)
  legacy <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    omega12 = 0, omega13 = 0, omega23 = 0,
    L11 = d[1], L12 = 0, L13 = 0,
    L21 = 0, L22 = d[2], L23 = 0,
    L31 = 0, L32 = 0, L33 = d[3],
    gamma12 = 0, gamma13 = 0, gamma23 = 0
  )
  z <- setNames(numeric(15L), admvn:::.sun_hs_z_names(3L, 3L))
  z[c("z41", "z52", "z63")] <- atanh(d)
  hs <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    z
  )
  dp_l <- make_sun_params(legacy)
  dp_h <- make_sun_hs_params(hs)
  expect_equal(dp_l$Omega, dp_h$Omega, tolerance = 1e-10)
  expect_equal(dp_l$Delta, dp_h$Delta, tolerance = 1e-10)
  expect_equal(dp_l$Gamma, dp_h$Gamma, tolerance = 1e-10)

  set.seed(5)
  x <- sn::rsun(6, dp = dp_l)
  v_l <- dsun(x, legacy, log = TRUE, n_points = 64L, n_shifts = 2L)$value
  v_h <- dsun_hs(x, hs, log = TRUE, n_points = 64L, n_shifts = 2L)$value
  expect_equal(v_l, v_h, tolerance = 1e-8)
})
