test_that("1D Matern 1/2 covariance matches the OU kernel in the interior", {
  knots <- seq(0, 20, by = 0.25)
  fem <- fem_bspline_1d(knots)
  kappa <- 1
  sd <- 1.2
  tau <- 1 / (sd * sqrt(2 * kappa))
  Q <- fem_precision(kappa, tau, fem$C, fem$G, alpha = 1L)
  x <- seq(6, 14, by = 0.5)
  A <- adlaplaceFem:::bspline_eval(fem$knots, x, 1L, 0L)
  Sigma <- as.matrix(A %*% Matrix::solve(Q, Matrix::t(A)))
  h <- abs(outer(x, x, "-"))
  target <- sd^2 * exp(-kappa * h)
  expect_equal(diag(Sigma), rep(sd^2, length(x)), tolerance = 0.05)
  expect_equal(Sigma, target, tolerance = 0.08)
  # Homogeneous Neumann boundaries inflate variance at the domain edge.
  x_b <- c(0, 20)
  A_b <- adlaplaceFem:::bspline_eval(fem$knots, x_b, 1L, 0L)
  var_b <- diag(as.matrix(A_b %*% Matrix::solve(Q, Matrix::t(A_b))))
  expect_gt(min(var_b), sd^2)
})

test_that("random_fem alpha=1 gradient matches finite differences", {
  knots <- seq(0, 4, by = 0.5)
  fem <- fem_bspline_1d(knots)
  prec <- fem_precision_payload(fem, alpha = 1L)
  nr <- fem$n
  range <- 1.5
  sd <- 0.7
  config <- list(
    gamma = rep(0, nr),
    theta = c(log(range), log(sd)),
    transform_theta = TRUE
  )
  rand_ssq <- adlaplace::density_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    density = "random_fem_ssq_1",
    package = "adlaplaceFem",
    precision = prec
  )
  rand_det <- adlaplace::density_data(
    beta_map = 0L,
    gamma_map = nr,
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "parameters",
    density = "random_fem_det_1",
    package = "adlaplaceFem",
    precision = prec
  )
  ptr <- c(
    adlaplace::ad_pack_ptr(rand_ssq, config),
    adlaplace::ad_pack_ptr(rand_det, config)
  )
  set.seed(3)
  x <- c(rnorm(nr, sd = 0.2), log(range), log(sd))
  kappa <- 2 / range
  tau <- 1 / (sd * sqrt(2 * kappa))
  Q <- fem_precision(kappa, tau, fem$C, fem$G, alpha = 1L)
  w <- x[seq_len(nr)]
  log_det <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  manual <- -0.5 * as.numeric(crossprod(w, Q %*% w)) +
    0.5 * log_det - 0.5 * nr * log(2 * pi)
  expect_equal(
    adlaplace::joint_log_dens(ptr, x, negative = FALSE),
    manual,
    tolerance = 1e-6
  )
  g_ad <- as.numeric(adlaplace::grad(ptr, x, inner = FALSE, negative = FALSE))
  eps <- 1e-5
  g_fd <- vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- xp[j] + eps
    xm[j] <- xm[j] - eps
    (adlaplace::joint_log_dens(ptr, xp, negative = FALSE) -
      adlaplace::joint_log_dens(ptr, xm, negative = FALSE)) / (2 * eps)
  }, numeric(1))
  expect_equal(g_ad, g_fd, tolerance = 1e-4)
})

test_that("rsmatern recovers range and sd of a simulated FEM field", {
  set.seed(11)
  knots <- seq(0, 10, by = 0.5)
  range <- 2
  sd <- 0.5
  kappa <- 2 / range
  tau <- 1 / (sd * sqrt(2 * kappa))
  fem <- fem_bspline_1d(knots)
  Q <- fem_precision(kappa, tau, fem$C, fem$G, alpha = 1L)
  R <- chol(as.matrix(Q))
  draw_w <- function() backsolve(R, rnorm(nrow(Q)))
  t <- seq(1.5, 8.5, by = 0.25)
  A <- adlaplaceFem:::bspline_eval(fem$knots, t, 1L, 0L)
  dat <- rbind(
    data.frame(
      time = t, region = "A", x = 1,
      y = as.numeric(A %*% draw_w()) + rnorm(length(t), sd = 0.02)
    ),
    data.frame(
      time = t, region = "B", x = 1,
      y = as.numeric(A %*% draw_w()) + rnorm(length(t), sd = 0.02)
    )
  )
  fit <- adlaplace::adlaplace(
    adlaplace::normal(y) ~ rsmatern(
      "time", mult = "x", by = "region", knots = knots,
      levels = c("A", "B"),
      init = c(1.2, 0.3),
      lower = c(0.5, 0.05),
      upper = c(6, 2)
    ),
    data = dat,
    control = list(maxit = 40, trace = 0)
  )
  mle <- stats::setNames(fit$par_info$mle, fit$par_info$label)
  expect_equal(unname(mle["time_x_rsmatern_range"]), range, tolerance = 0.35)
  expect_equal(unname(mle["time_x_rsmatern_sd"]), sd, tolerance = 0.35)
})

test_that("replace_vars retargets the rsmatern exposure column", {
  term <- rsmatern(
    "time", mult = "sqrt_pm", by = "region",
    knots = c(0, 1, 2), levels = "A"
  )
  out <- adlaplace::replace_vars(term, c(sqrt_pm = "sqrt_pm_s3"))
  expect_equal(out@mult, "sqrt_pm_s3")
  expect_equal(out@name, "time")
  expect_true(grepl("sqrt_pm_s3", out@label))
})
