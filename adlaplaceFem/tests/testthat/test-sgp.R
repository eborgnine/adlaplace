test_that("sgp constructor builds the FEM term and its boundary", {
  one <- sgp("x", period = 2, k = 8, region = c(0, 4))
  expect_true(inherits(one$x_sgp, "sgp"))
  expect_true(inherits(one$x_sgp_boundary, "sgp_boundary"))
  expect_equal(one$x_sgp@a, pi)
  expect_equal(one$x_sgp@m, 1L)
  expect_equal(one$x_sgp@ref_value, 0)
  expect_equal(one$x_sgp@fem$rank, 3L * (8L - 2L))
  expect_null(theta_info(one$x_sgp_boundary))

  knots <- seq(0, 4, length.out = 6)
  from_knots <- sgp("x", a = pi, knots = knots)
  expect_equal(from_knots$x_sgp@fem$log_det, one$x_sgp@fem$log_det)
  expect_equal(as.matrix(from_knots$x_sgp@fem$Q), as.matrix(one$x_sgp@fem$Q))

  two <- sgp("x", a = 1, knots = knots, m = 2)
  expect_equal(two$x_sgp@fem$rank, 2L * 3L * length(knots))
  Q <- as.matrix(two$x_sgp@fem$Q)
  n1 <- ncol(Q) / 2
  expect_equal(Q[seq_len(n1), seq_len(n1) + n1], matrix(0, n1, n1))

  renamed <- replace_vars(one, c(x = "z"))
  expect_identical(names(renamed), c("z_sgp", "z_sgp_boundary"))
  expect_identical(renamed$z_sgp@label, "z_sgp")

  expect_error(sgp("x", period = 1, a = 1, k = 8, region = c(0, 1)), "exactly one")
  expect_error(sgp("x", period = 1, k = 3, region = c(0, 1)), "k must")
  expect_error(sgp("x", period = 1), "knots")
})

test_that("sgp design and precision have matching dimensions", {
  k <- 8L
  term <- sgp("x", period = 1, k = k, region = c(0, 3))$x_sgp
  dat <- data.frame(x = seq(0, 3, length.out = 20))
  B <- design(term, dat)
  prec <- precision(term, dat)
  info <- random_info(term, dat)
  expect_equal(ncol(B), 3L * (k - 2L))
  expect_equal(ncol(B), prec$rank)
  expect_equal(nrow(prec$Q), prec$rank)
  expect_true(is.finite(prec$log_det))
  expect_equal(colnames(B), info$gamma_label)
  expect_identical(theta_info(term)$label, "x_sgp")
  expect_error(design(term, data.frame(x = c(0, 4))), "outside")
})

test_that("sgp boundary is sin and cos, flat unless boundary_sd is finite", {
  a <- 1.5
  m <- 2L
  terms <- sgp("x", a = a, knots = seq(1, 5, length.out = 6), m = m)
  bterm <- terms$x_sgp_boundary
  x <- seq(1, 5, length.out = 12)
  D <- as.matrix(design(bterm, data.frame(x = x)))
  u <- x - 1
  expect_equal(ncol(D), 2L * m)
  expect_equal(D[, 1], cos(a * u))
  expect_equal(D[, 2], sin(a * u))
  expect_equal(D[, 3], cos(2 * a * u))
  expect_equal(D[, 4], sin(2 * a * u))
  expect_equal(colnames(D), random_info(bterm, data.frame(x = x))$gamma_label)
  expect_null(precision(bterm, NULL))

  known <- sgp(
    "x", a = a, knots = seq(1, 5, length.out = 6),
    m = m, boundary_sd = 2
  )$x_sgp_boundary
  P <- as.matrix(precision(known, NULL))
  expect_equal(unname(diag(P)), rep(0.25, 2L * m))
  expect_equal(unname(P[upper.tri(P)]), rep(0, sum(upper.tri(P))))
})

test_that("sgp FEM Q and design match BayesGP reference (fda)", {
  skip_if_not_installed("fda")
  region <- c(0, 2)
  k <- 8L
  x <- seq(0, 2, length.out = 17)
  for (a in c(1, pi)) {
    Q_ref <- as.matrix(Compute_Q_sB(
      a = a, k = k, region = region, accuracy = 1e-4, boundary = TRUE
    ))
    basis <- adlaplaceFem:::sgp_basis(seq(region[[1L]], region[[2L]], length.out = k - 2L))
    Q_new <- as.matrix(adlaplaceFem:::sgp_harmonic_precision(
      basis, omega = a, n_quad = 24L
    ))
    rel <- max(abs(Q_new - Q_ref)) / max(abs(Q_ref))
    expect_lt(rel, 1e-3, label = paste("Q a =", a))

    B_ref <- Compute_B_sB(x, a = a, k = k, region = region, boundary = TRUE)
    fem <- adlaplaceFem:::sgp_fem(
      knots = seq(region[[1L]], region[[2L]], length.out = k - 2L),
      a = a, m = 1L, n_quad = 8L
    )
    B_new <- as.matrix(adlaplaceFem:::sgp_fem_design(x, fem))
    expect_equal(unname(B_new), unname(B_ref), tolerance = 1e-8, label = paste("B a =", a))
  }
})

test_that("psd(sgp) matches the seasonal predictive SD", {
  h <- c(0.2, 0.75, 1.5)
  a <- 2 * pi / 5
  term <- sgp("x", a = a, k = 8, region = c(0, 10))$x_sgp
  expect_equal(psd(term, h = h, x = 3), compute_d_step_sgpsd(h, a))
  expect_equal(psd(term, h = h, sd = 2), 2 * compute_d_step_sgpsd(h, a))

  mterm <- sgp("x", a = a, k = 8, region = c(0, 10), m = 2)$x_sgp
  c1 <- compute_d_step_sgpsd(h, a)
  c2 <- compute_d_step_sgpsd(h, 2 * a)
  expect_equal(psd(mterm, h = h), sqrt(c1^2 + c2^2))

  pr <- psd_prior(sgp("x", a = a, k = 8, region = c(0, 10)), u = 0.2, alpha = 0.01, h = 1)
  expect_s4_class(pr, "prior")
  expect_identical(pr@theta, "x_sgp")
  unit <- psd(term, h = 1)
  expect_equal(pr@rate, -log(0.01) * unit / 0.2)
})

test_that("end-to-end sgp fit on a sine has finite logLik", {
  skip_on_cran()
  set.seed(1)
  x <- seq(0, 10, length.out = 40)
  period <- 2
  dat <- data.frame(
    x = x,
    y = sin(2 * pi * x / period) + rnorm(length(x), sd = 0.05)
  )
  fit <- adlaplace(
    normal(y) ~ 0 + sgp(x, period = period, k = 6, region = c(0, 10)),
    data = dat,
    config = list(num_threads = 1L, num_shards = 1L),
    control = list(maxit = 40L),
    verbose = FALSE,
    hessian = FALSE
  )
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})
