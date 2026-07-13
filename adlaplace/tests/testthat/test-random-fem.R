test_that("random_fem_2 matches dense Matérn FEM log density", {
  skip_if_not_installed("adlaplaceGrf")

  g <- expand.grid(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  fem <- adlaplaceGrf::grf_bspline(list(x = g$x, y = g$y), degree = 2L, n_interior = 1L)
  prec <- adlaplaceGrf::fem_precision_payload(fem, alpha = 2L)
  nr <- nrow(fem$C)

  model <- adlaplace:::ad_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    ad_fun = "random_fem_2",
    precision = prec
  )
  log_kappa <- log(1.5)
  log_tau <- log(0.8)
  config <- list(
    gamma = rep(0, nr),
    theta = c(log_kappa, log_tau),
    transform_theta = TRUE
  )
  ptr <- adlaplace::ad_fun_ptr(model, config)

  set.seed(7)
  w <- rnorm(nr)
  x <- c(w, log_kappa, log_tau)
  Q <- adlaplaceGrf::fem_precision(
    kappa = exp(log_kappa), tau = exp(log_tau),
    fem$C, fem$G, fem$G2, alpha = 2L
  )
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  manual <- 0.5 * log_det - 0.5 * as.numeric(crossprod(w, Q %*% w)) -
    0.5 * nr * log(2 * pi)

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
