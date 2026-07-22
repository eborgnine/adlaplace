test_that("random_fem ssq + det match dense Matern FEM log density", {
  knots_list <- list(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  )
  g <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(g, knots_list, degree = 2L)
  prec <- fem_precision_payload(fem, alpha = 2L)
  nr <- nrow(fem$C)

  # Public theta: (range, sd) with range = sqrt(8*nu)/kappa; alpha=2 => nu=1
  kappa <- 1.5
  nu <- 1
  range <- sqrt(8 * nu) / kappa
  sd <- 0.8
  tau <- 1 / (kappa * sd * sqrt(4 * pi))
  log_range <- log(range)
  log_sd <- log(sd)
  config <- list(
    gamma = rep(0, nr),
    theta = c(log_range, log_sd),
    transform_theta = TRUE
  )

  rand_ssq <- adlaplace::density_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    density = "random_fem_ssq_2",
    package = "adlaplaceFem",
    precision = prec
  )
  rand_det <- adlaplace::density_data(
    beta_map = 0L,
    gamma_map = nr,
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "parameters",
    density = "random_fem_det_2",
    package = "adlaplaceFem",
    precision = prec
  )

  ptr_ssq <- adlaplace::ad_pack_ptr(rand_ssq, config)
  ptr_det <- adlaplace::ad_pack_ptr(rand_det, config)
  expect_equal(adlaplace:::n_groups(ptr_ssq), 1L)
  expect_equal(adlaplace:::n_groups(ptr_det), 1L)

  set.seed(7)
  w <- rnorm(nr)
  x <- c(w, log_range, log_sd)
  Q <- fem_precision(
    kappa = kappa, tau = tau,
    fem$C, fem$G, fem$G2, alpha = 2L
  )
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  qf <- as.numeric(crossprod(w, Q %*% w))
  manual_ssq <- -0.5 * qf
  manual_det <- 0.5 * log_det - 0.5 * nr * log(2 * pi)
  manual <- manual_ssq + manual_det

  expect_equal(
    adlaplace::joint_log_dens(ptr_ssq, x, negative = FALSE),
    manual_ssq,
    tolerance = 1e-8
  )
  expect_equal(
    adlaplace::joint_log_dens(ptr_det, x, negative = FALSE),
    manual_det,
    tolerance = 1e-8
  )
  # c() moves shards and clears inputs
  ptr <- c(ptr_ssq, ptr_det)
  expect_equal(
    adlaplace::joint_log_dens(ptr, x, negative = FALSE),
    manual,
    tolerance = 1e-8
  )

  g_ad <- as.numeric(adlaplace::grad(ptr, x, inner = FALSE, negative = FALSE))
  eps <- 1e-5
  g_fd <- vapply(seq_along(x), function(i) {
    xp <- x
    xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    (adlaplace::joint_log_dens(ptr, xp, negative = FALSE) -
      adlaplace::joint_log_dens(ptr, xm, negative = FALSE)) / (2 * eps)
  }, numeric(1))
  expect_equal(g_ad, g_fd, tolerance = 1e-4)
})

test_that("random_fem_det theta Hessian matches finite differences", {
  knots_list <- list(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  )
  g <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(g, knots_list, degree = 2L)
  prec <- fem_precision_payload(fem, alpha = 2L)
  nr <- nrow(fem$C)

  log_range <- log(1.3)
  log_sd <- log(0.8)
  config <- list(
    gamma = rep(0, nr),
    theta = c(log_range, log_sd),
    transform_theta = TRUE
  )
  rand_det <- adlaplace::density_data(
    beta_map = 0L,
    gamma_map = nr,
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "parameters",
    density = "random_fem_det_2",
    package = "adlaplaceFem",
    precision = prec
  )
  ptr <- adlaplace::ad_pack_ptr(rand_det, config)
  x <- c(rep(0, nr), log_range, log_sd)
  idx <- nr + 1:2

  gr <- function(xx) {
    as.numeric(adlaplace::grad(ptr, xx, inner = FALSE, negative = FALSE))
  }
  H_ad <- as.matrix(adlaplace::hessian(ptr, x, negative = FALSE))[idx, idx]
  eps <- 1e-5
  H_fd <- matrix(0, 2, 2)
  for (a in 1:2) {
    xp <- x
    xm <- x
    xp[idx[a]] <- xp[idx[a]] + eps
    xm[idx[a]] <- xm[idx[a]] - eps
    H_fd[a, ] <- (gr(xp)[idx] - gr(xm)[idx]) / (2 * eps)
  }
  expect_equal(H_ad, H_fd, tolerance = 1e-5)
})

test_that("random_fem alpha=3 det matches dense logdet, grad and Hessian FD", {
  knots_list <- list(
    x = seq(0, 1, length.out = 5),
    y = seq(0, 1, length.out = 5)
  )
  g <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(g, knots_list, degree = 3L)
  prec <- fem_precision_payload(fem, alpha = 3L)
  nr <- nrow(fem$C)

  kappa <- 2.0
  nu <- 2
  range <- sqrt(8 * nu) / kappa
  sd <- 0.6
  tau <- 1 / (kappa^2 * sd * sqrt(8 * pi))
  log_range <- log(range)
  log_sd <- log(sd)
  config <- list(
    gamma = rep(0, nr),
    theta = c(log_range, log_sd),
    transform_theta = TRUE
  )
  rand_det <- adlaplace::density_data(
    beta_map = 0L,
    gamma_map = nr,
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "parameters",
    density = "random_fem_det_3",
    package = "adlaplaceFem",
    precision = prec
  )
  ptr <- adlaplace::ad_pack_ptr(rand_det, config)
  x <- c(rep(0, nr), log_range, log_sd)
  idx <- nr + 1:2

  Q <- fem_precision(
    kappa = kappa, tau = tau,
    fem$C, fem$G, fem$G2, fem$G3, alpha = 3L
  )
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  manual_det <- 0.5 * log_det - 0.5 * nr * log(2 * pi)
  expect_equal(
    adlaplace::joint_log_dens(ptr, x, negative = FALSE),
    manual_det,
    tolerance = 1e-8
  )

  f <- function(xx) adlaplace::joint_log_dens(ptr, xx, negative = FALSE)
  g_ad <- as.numeric(adlaplace::grad(ptr, x, inner = FALSE, negative = FALSE))
  eps <- 1e-5
  g_fd <- vapply(idx, function(i) {
    xp <- x
    xm <- x
    xp[i] <- xp[i] + eps
    xm[i] <- xm[i] - eps
    (f(xp) - f(xm)) / (2 * eps)
  }, numeric(1))
  expect_equal(g_ad[idx], g_fd, tolerance = 1e-4)

  gr <- function(xx) {
    as.numeric(adlaplace::grad(ptr, xx, inner = FALSE, negative = FALSE))
  }
  H_ad <- as.matrix(adlaplace::hessian(ptr, x, negative = FALSE))[idx, idx]
  H_fd <- matrix(0, 2, 2)
  for (a in 1:2) {
    xp <- x
    xm <- x
    xp[idx[a]] <- xp[idx[a]] + eps
    xm[idx[a]] <- xm[idx[a]] - eps
    H_fd[a, ] <- (gr(xp)[idx] - gr(xm)[idx]) / (2 * eps)
  }
  expect_equal(H_ad, H_fd, tolerance = 1e-5)
})
