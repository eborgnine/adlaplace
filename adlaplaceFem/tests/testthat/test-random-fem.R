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

  rand_ssq <- adlaplace::ad_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    ad_fun = "random_fem_ssq_2",
    package = "adlaplaceFem",
    precision = prec
  )
  rand_det <- adlaplace::ad_data(
    beta_map = 0L,
    gamma_map = nr,
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "parameters",
    ad_fun = "random_fem_det_2",
    package = "adlaplaceFem",
    precision = prec
  )

  ptr_ssq <- adlaplace::ad_fun_ptr(rand_ssq, config)
  ptr_det <- adlaplace::ad_fun_ptr(rand_det, config)
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
