test_that("make_sun_params uses free Lij entries", {
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    omega12 = 0, omega13 = 0, omega23 = 0,
    L11 = 0.4, L12 = 0.1, L13 = 0.05,
    L21 = -0.1, L22 = 0.4, L23 = 0.2,
    L31 = 0, L32 = 0.15, L33 = 0.35,
    gamma12 = 0, gamma13 = 0, gamma23 = 0
  )
  dp <- make_sun_params(par)
  L <- matrix(
    c(
      par["L11"], par["L12"], par["L13"],
      par["L21"], par["L22"], par["L23"],
      par["L31"], par["L32"], par["L33"]
    ),
    3L, byrow = TRUE
  )
  # Cu = I => Delta = C_u L^T = L^T
  expect_equal(dp$Delta, t(L), tolerance = 1e-12)
})

test_that("make_sun_params matches expected structure", {
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
  expect_length(dp$xi, 3L)
  expect_equal(dim(dp$Omega), c(3L, 3L))
  expect_equal(dim(dp$Delta), c(3L, 3L))
  expect_equal(dp$tau, c(0, 0, 0))
  expect_equal(dim(dp$Gamma), c(3L, 3L))
  expect_equal(diag(dp$Gamma), rep(1, 3L), tolerance = 1e-10)
})

test_that("unit-row Cholesky yields a correlation matrix without cov2cor", {
  Br <- admvn:::.unit_row_chol(c(0.3, 0.2, 0.3), 3L)
  expect_equal(rowSums(Br^2), rep(1, 3L), tolerance = 1e-12)
  R <- Br %*% t(Br)
  expect_equal(diag(R), rep(1, 3L), tolerance = 1e-12)
  expect_true(all(abs(R) <= 1 + 1e-12))
  # Zero free params -> identity
  expect_equal(admvn:::.unit_row_chol(c(0, 0, 0), 3L), diag(3))
})

test_that("make_sun_params Sigma_vv has unit diagonal via residual R", {
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
  expect_equal(diag(dp$Gamma), rep(1, 3L), tolerance = 1e-10)
  expect_equal(unname(diag(dp$Omega)), unname(c(par["nu1"], par["nu2"], par["nu3"])^2),
               tolerance = 1e-10)
  Cu <- admvn:::.corr_from_free(
    c(par["omega12"], par["omega13"], par["omega23"]), 3L
  )
  L <- matrix(
    c(
      par["L11"], par["L12"], par["L13"],
      par["L21"], par["L22"], par["L23"],
      par["L31"], par["L32"], par["L33"]
    ),
    3L, byrow = TRUE
  )
  M <- L %*% Cu %*% t(L)
  R <- admvn:::.residual_corr_R(M, c(par["gamma12"], par["gamma13"], par["gamma23"]))
  expect_equal(dp$Gamma, M + R, tolerance = 1e-10)
  expect_true(all(diag(M) < 1))
  expect_equal(dp$Delta, Cu %*% t(L), tolerance = 1e-10)
})

test_that("make_sun_params with L=0 has Gamma equal to residual correlation C", {
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    omega12 = 0, omega13 = 0, omega23 = 0,
    L11 = 0, L12 = 0, L13 = 0,
    L21 = 0, L22 = 0, L23 = 0,
    L31 = 0, L32 = 0, L33 = 0,
    gamma12 = 0.4, gamma13 = -0.2, gamma23 = 0.1
  )
  dp <- make_sun_params(par)
  Br <- admvn:::.unit_row_chol(c(0.4, -0.2, 0.1), 3L)
  C <- Br %*% t(Br)
  expect_equal(dp$Gamma, C, tolerance = 1e-12)
  expect_equal(diag(dp$Gamma), rep(1, 3L), tolerance = 1e-12)
})

test_that("dsun log-density matches sn::dsun", {
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
  set.seed(1)
  x <- sn::rsun(8, dp = dp)
  sn_val <- sum(sn::dsun(x, dp = dp, log = TRUE))
  # d<=3 orthants use specialized CDF; QMC budget is irrelevant for values.
  adv <- admvn::dsun(x, par, log = TRUE, deriv = 0L, n_points = 64L, n_shifts = 2L)
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
  expect_null(adv$gradient)
})

test_that("dsun gradient matches numDeriv on summed loglik", {
  skip_if_not_installed("sn")
  skip_if_not_installed("numDeriv")
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
  set.seed(2)
  x <- sn::rsun(4, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1500L, n_shifts = 10L, n_threads = 1L)
  f <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
  g_num <- numDeriv::grad(f, par, method = "Richardson")
  g_ad <- tape$eval(par, log = TRUE, deriv = 1L)$gradient
  expect_equal(g_ad, as.numeric(g_num), tolerance = 0.35)
})

test_that("SUN grads w.r.t. L and skewness params match finite differences", {
  skip_if_not_installed("sn")
  # Regression: missing atomic jac_sparsity previously zeroed these grads.
  par <- c(
    xi1 = 0, xi2 = 0, xi3 = 0,
    nu1 = 1, nu2 = 1, nu3 = 1,
    omega12 = 0.2, omega13 = 0.1, omega23 = 0.2,
    L11 = 0.4, L12 = 0.1, L13 = 0.1,
    L21 = 0.1, L22 = 0.4, L23 = 0.1,
    L31 = 0, L32 = 0.1, L33 = 0.35,
    gamma12 = 0.3, gamma13 = 0.2, gamma23 = 0.3
  )
  set.seed(7)
  x <- sn::rsun(12, dp = make_sun_params(par))
  tape <- dsun_fun(x, par, n_points = 64L, n_shifts = 2L, n_threads = 1L)
  g_ad <- tape$eval(par, log = TRUE, deriv = 1L)$gradient
  ll <- function(p) tape$eval(p, log = TRUE, deriv = 0L)$value
  idx <- match(c("L11", "L22", "L33", "L12", "gamma12", "gamma13", "gamma23"), names(par))
  expect_true(all(abs(g_ad[idx]) > 1e-8))
  for (j in idx) {
    h <- pmax(1e-6, 1e-5 * abs(par[j]))
    pu <- pd <- par
    pu[j] <- par[j] + h
    pd[j] <- par[j] - h
    g_fd <- (ll(pu) - ll(pd)) / (2 * h)
    expect_equal(g_ad[j], g_fd, tolerance = 5e-3)
  }
})

test_that("dsun_fun reuses tape with fixed data", {
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
  x <- sn::rsun(3, dp = dp)
  f <- dsun_fun(x, par, n_points = 1500L, n_shifts = 10L, seed = 7L, n_threads = 1L)
  direct <- admvn::dsun(
    x, par, log = TRUE, deriv = 1L,
    n_points = 1500L, n_shifts = 10L, seed = 7L, n_threads = 1L
  )
  via_tape <- f$eval(par, log = TRUE, deriv = 1L)
  expect_equal(via_tape$value, direct$value, tolerance = 1e-10)
  expect_equal(via_tape$gradient, direct$gradient, tolerance = 1e-8)
})

test_that("serial and parallel shard eval agree", {
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
  set.seed(3)
  x <- sn::rsun(6, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1200L, n_shifts = 8L, seed = 11L, n_threads = 1L)
  serial <- tape$eval(par, log = TRUE, deriv = 1L, n_threads = 1L)
  parallel <- tape$eval(par, log = TRUE, deriv = 1L, n_threads = 4L)
  expect_equal(parallel$value, serial$value, tolerance = 1e-10)
  expect_equal(parallel$gradient, serial$gradient, tolerance = 1e-8)
})

test_that("optim_fns caches fn and gr for the same par", {
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
  x <- sn::rsun(3, dp = dp)
  tape <- dsun_fun(x, par, n_points = 1200L, n_shifts = 8L, n_threads = 1L)
  obj <- tape$optim_fns()
  direct <- tape$eval(par, log = TRUE, deriv = 1L)
  expect_equal(obj$fn(par), direct$value, tolerance = 1e-10)
  expect_equal(obj$gr(par), direct$gradient, tolerance = 1e-8)
  expect_equal(obj$fn(par), direct$value, tolerance = 1e-10)
  expect_equal(obj$gr(par), direct$gradient, tolerance = 1e-8)
})

test_that("sun_mle improves loglik from a perturbed start", {
  skip_if_not_installed("sn")
  skip_if_not_installed("trustOptim")
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
  set.seed(5)
  x <- sn::rsun(20, dp = dp)
  start <- par
  start[1:3] <- c(0.2, -0.1, 0.1)
  start[4:6] <- c(1.2, 0.8, 1.1)
  tape <- dsun_fun(x, start, n_points = 400L, n_shifts = 6L, n_threads = 1L)
  ll0 <- tape$eval(start, log = TRUE, deriv = 0L)$value
  fit <- sun_mle(tape, start, control = list(maxit = 20L))
  expect_true(is.finite(fit$value))
  expect_gt(fit$value, ll0 - 1e-6)
  expect_true(all(fit$par[4:6] > 0))
})

test_that("sun_mle reaches the same optimum as optim on a small sample", {
  skip_if_not_installed("sn")
  skip_if_not_installed("trustOptim")
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
  set.seed(6)
  x <- sn::rsun(12, dp = dp)
  start <- par
  start[1:3] <- c(0.15, -0.05, 0.08)
  start[4:6] <- c(1.1, 0.9, 1.05)
  tape <- dsun_fun(x, start, n_points = 400L, n_shifts = 6L, n_threads = 1L)
  obj <- tape$optim_fns()
  lower <- c(rep(-Inf, 3), rep(1e-3, 3), rep(-Inf, 15))
  fit_opt <- optim(
    start, fn = obj$fn, gr = obj$gr,
    method = "L-BFGS-B", lower = lower, upper = rep(Inf, 21L),
    control = list(fnscale = -1, maxit = 40L)
  )
  fit_tr <- sun_mle(tape, start, control = list(maxit = 40L))
  # Same objective; trust-region vs L-BFGS-B may land in different basins.
  expect_equal(fit_tr$value, fit_opt$value, tolerance = 0.5)
  ll_tr <- tape$eval(fit_tr$par, log = TRUE, deriv = 0L)$value
  ll_opt <- tape$eval(fit_opt$par, log = TRUE, deriv = 0L)$value
  expect_equal(ll_tr, fit_tr$value, tolerance = 1e-8)
  expect_equal(ll_opt, fit_opt$value, tolerance = 1e-6)
  expect_gt(ll_tr, tape$eval(start, log = TRUE, deriv = 0L)$value - 1e-6)
  expect_gt(ll_opt, tape$eval(start, log = TRUE, deriv = 0L)$value - 1e-6)
})

test_that("unit weights match unweighted dsun", {
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
  set.seed(11)
  x <- sn::rsun(6, dp = make_sun_params(par))
  unwt <- dsun(x, par, log = TRUE, deriv = 1L, n_points = 64L, n_shifts = 2L)
  unit <- dsun(
    x, par, log = TRUE, deriv = 1L, n_points = 64L, n_shifts = 2L,
    weights = rep(1, nrow(x))
  )
  expect_equal(unit$value, unwt$value, tolerance = 1e-10)
  expect_equal(unit$gradient, unwt$gradient, tolerance = 1e-10)
})

test_that("weighted dsun matches sum of weighted single-row evals", {
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
  set.seed(12)
  x <- sn::rsun(5, dp = make_sun_params(par))
  w <- c(0.5, 1.2, 0.8, 1.5, 0.3)
  joint <- dsun(
    x, par, log = TRUE, deriv = 1L, n_points = 64L, n_shifts = 2L,
    weights = w
  )
  val <- 0
  grad <- numeric(21L)
  for (i in seq_len(nrow(x))) {
    one <- dsun(
      x[i, , drop = FALSE], par, log = TRUE, deriv = 1L,
      n_points = 64L, n_shifts = 2L
    )
    val <- val + w[i] * one$value
    grad <- grad + w[i] * one$gradient
  }
  expect_equal(joint$value, val, tolerance = 1e-8)
  expect_equal(joint$gradient, grad, tolerance = 1e-8)
})

test_that("weighted dsun value matches weighted sn::dsun", {
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
  set.seed(13)
  x <- sn::rsun(8, dp = dp)
  w <- runif(nrow(x), 0.2, 1.8)
  sn_val <- sum(w * sn::dsun(x, dp = dp, log = TRUE))
  adv <- dsun(
    x, par, log = TRUE, deriv = 0L, n_points = 64L, n_shifts = 2L,
    weights = w
  )
  expect_equal(adv$value, sn_val, tolerance = 1e-4)
})

test_that("weighted dsun_fun tape matches dsun weights", {
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
  set.seed(14)
  x <- sn::rsun(4, dp = make_sun_params(par))
  w <- c(1.1, 0.4, 2.0, 0.7)
  tape <- dsun_fun(
    x, par, n_points = 64L, n_shifts = 2L, n_threads = 1L, weights = w
  )
  from_tape <- tape$eval(par, log = TRUE, deriv = 1L)
  from_dsun <- dsun(
    x, par, log = TRUE, deriv = 1L, n_points = 64L, n_shifts = 2L,
    weights = w
  )
  expect_equal(from_tape$value, from_dsun$value, tolerance = 1e-8)
  expect_equal(from_tape$gradient, from_dsun$gradient, tolerance = 1e-8)
})
