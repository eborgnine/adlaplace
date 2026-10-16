test_that("mgp constructor builds one- and two-sided FEM terms", {
  one <- mgp("x", a = 2, knots = seq(1, 5, length.out = 10), ref_value = "left")
  expect_true(inherits(one$x_mgp, "mgp"))
  expect_true(inherits(one$x_mgp_boundary, "mgp_boundary"))
  expect_false(one$x_mgp@two_sided)
  expect_equal(one$x_mgp@ref_value, 1)
  expect_equal(one$x_mgp@shift, 1)

  knots_two <- sort(unique(c(
    seq(8, 40, length.out = 10),
    seq(40, 100, length.out = 10)
  )))
  two <- mgp("x", a = 2, knots = knots_two, ref_value = 40)
  expect_true(two$x_mgp@two_sided)
  expect_equal(two$x_mgp@fem$pos_len, 60)
  expect_equal(two$x_mgp@fem$neg_len, 32)
  expect_equal(ncol(two$x_mgp@fem$Q), 2L * 10L)
})

test_that("mgp design and precision have matching dimensions", {
  term <- mgp("x", a = 2, knots = seq(1, 4, length.out = 8), ref_value = "left")$x_mgp
  dat <- data.frame(x = seq(1, 4, length.out = 20))
  B <- design(term, dat)
  prec <- precision(term, dat)
  expect_equal(ncol(B), prec$rank)
  expect_equal(nrow(prec$Q), prec$rank)
  expect_true(is.finite(prec$log_det))
  expect_equal(colnames(B), paste0("x_mgp_k", formatC(seq_len(ncol(B)), width = 1, flag = "0")))
})

test_that("mgp FEM Q and design match BayesGP reference (fda)", {
  skip_if_not_installed("fda")
  region <- c(0, 2)
  k <- 10
  # Positive side: any cs > 0. Reverse side needs region_max < cs so α_rev stays finite.
  cs_pos <- 0.5
  cs_rev <- 3
  for (a in c(2, 1, -1, 3)) {
    Q_ref <- as.matrix(Compute_Prec(
      a = a, c = cs_pos, k = k, region = region, accuracy = 1e-4, boundary = TRUE
    ))
    basis <- adlaplaceFem:::mgp_side_basis(max(region), k)
    Q_new <- as.matrix(adlaplaceFem:::mgp_side_precision(
      basis$knots, a = a, cs = cs_pos, side = "positive", n_quad = 16L, drop = basis$drop
    ))
    rel <- max(abs(Q_new - Q_ref)) / max(abs(Q_ref))
    expect_lt(rel, 1e-3, label = paste("Q a =", a))

    Q_rev_ref <- as.matrix(Compute_Prec_rev(
      a = a, c = cs_rev, k = k, region = region, accuracy = 1e-4, boundary = TRUE
    ))
    Q_rev <- as.matrix(adlaplaceFem:::mgp_side_precision(
      basis$knots, a = a, cs = cs_rev, side = "reverse", n_quad = 16L, drop = basis$drop
    ))
    rel_rev <- max(abs(Q_rev - Q_rev_ref)) / max(abs(Q_rev_ref))
    expect_lt(rel_rev, 1e-3, label = paste("Q_rev a =", a))
  }

  x <- seq(0, 2, length.out = 17)
  B_ref <- Compute_Design(x, k = k, region = region, boundary = TRUE)
  B_new <- as.matrix(adlaplaceFem:::mgp_side_design(x, basis$knots, drop = basis$drop))
  expect_equal(unname(B_new), unname(B_ref), tolerance = 1e-10)
})

test_that("mgp boundary matches boundary_basis_matrix", {
  knots <- seq(8, 100, length.out = 8)
  for (ref in list("left", "middle", "right", 40)) {
    terms <- mgp("x", a = 2, knots = knots, ref_value = ref)
    bterm <- terms$x_mgp_boundary
    x <- seq(min(knots), max(knots), length.out = 25)
    got <- as.numeric(design(bterm, data.frame(x = x)))
    ref_mat <- boundary_basis_matrix(
      x = x,
      a = 2,
      c = bterm@shift,
      initial_location = bterm@ref_value,
      normalized = TRUE
    )
    expect_equal(got, as.numeric(ref_mat), tolerance = 1e-10, label = paste("ref", ref))
  }
})

test_that("psd(mgp) matches SS_cov", {
  terms <- mgp("x", a = 2, knots = seq(1, 5, length.out = 10), ref_value = "left")
  term <- terms$x_mgp
  for (a in c(2, 1, -1, 3, -3)) {
    term_a <- term
    term_a@a <- a
    x0 <- term@ref_value
    h <- c(0.2, 0.5, 1)
    x <- x0 + c(0, 0.3, 1)
    got <- psd(term_a, h = h, x = x)
    u <- x - x0
    expect_equal(
      got,
      sqrt(mapply(function(ui, hi) SS_cov(ui + hi, ui, c = term@shift, a = a), u, h)),
      tolerance = 1e-5,
      label = paste("a =", a)
    )
  }
})

test_that("mgp FEM prior SD approaches psd as knots densify", {
  knots <- seq(1, 6, length.out = 58)
  x0 <- 1
  h <- 0.75
  a <- 2
  terms <- mgp("x", a = a, knots = knots, ref_value = x0)
  term <- terms$x_mgp
  target <- psd(term, h = h, x = x0)
  B <- as.matrix(design(term, data.frame(x = x0 + h)))
  Q <- as.matrix(precision(term, NULL)$Q)
  fem_sd <- sqrt(as.numeric(B %*% solve(Q, t(B))))
  expect_equal(fem_sd, target, tolerance = 5e-2)
})

test_that("psd_prior(mgp) matches prior_conversion_mgp rate", {
  terms <- mgp("pm10", a = 2, knots = seq(8, 100, length.out = 10), ref_value = "left")
  u <- 0.01
  alpha <- 0.5
  h <- 5
  x <- 8
  pr <- psd_prior(terms, u = u, alpha = alpha, h = h, x = x)
  expect_identical(pr@theta, "pm10_mgp")
  converted <- prior_conversion_mgp(
    prior = list(u = u, alpha = alpha),
    h = h,
    c = terms$pm10_mgp@shift,
    a = 2,
    x = x - terms$pm10_mgp@ref_value
  )
  expect_equal(pr@rate, -log(alpha) / converted$u, tolerance = 1e-8)
})

test_that("random_mult tape from mgp matches dense density", {
  skip_if_not_installed("numDeriv")
  terms <- mgp("x", a = 2, knots = seq(1, 3, length.out = 6), ref_value = "left")
  term <- terms$x_mgp
  prec <- precision(term, NULL)
  nr <- prec$rank
  model <- adlaplace:::density_data(
    gamma_map = Matrix::sparseMatrix(
      i = seq_len(nr), j = seq_len(nr), x = 1, dims = c(nr, nr)
    ),
    theta_map = c(1L, 1L),
    ad_kind = "random",
    density = "random_mult",
    precision = prec
  )
  log_sd <- log(0.4)
  config <- list(gamma = rep(0, nr), theta = log_sd, transform_theta = TRUE)
  ptr <- adlaplace::ad_pack_ptr(model, config)
  set.seed(1)
  u <- rnorm(nr)
  x <- c(u, log_sd)
  tau <- exp(-2 * log_sd)
  Q <- as.matrix(prec$Q)
  manual <- 0.5 * (nr * log(tau) + prec$log_det) -
    0.5 * tau * as.numeric(crossprod(u, Q %*% u)) -
    0.5 * nr * log(2 * pi)
  expect_equal(
    adlaplace::joint_log_dens(ptr, x, negative = FALSE),
    manual,
    tolerance = 1e-8
  )
  g <- as.numeric(adlaplace::grad(ptr, x, inner = FALSE, negative = FALSE))
  g_nd <- numDeriv::grad(
    function(z) adlaplace::joint_log_dens(ptr, z, negative = FALSE),
    x
  )
  expect_equal(g, g_nd, tolerance = 1e-5)
})

test_that("end-to-end mgp fit on a sqrt curve has finite logLik", {
  skip_on_cran()
  set.seed(11)
  x <- seq(1, 9, length.out = 80)
  dat <- data.frame(x = x, y = sqrt(x) + rnorm(length(x), sd = 0.05))
  fit <- adlaplace(
    normal(y) ~ 0 + mgp(x, a = 2, knots = seq(1, 9, length.out = 14), ref_value = "left"),
    data = dat,
    config = list(num_threads = 1L, num_shards = 1L),
    control = list(maxit = 40L),
    verbose = FALSE,
    hessian = FALSE
  )
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("london-2003 case-crossover with mgp returns finite logLik", {
  skip_on_cran()
  data("london", package = "adlaplace")
  london <- london[!is.na(london$pm10) & london$year == 2003, ]
  fit <- adlaplace(
    dirichlet_multinom(all, by = c("year", "month", "dow"), init = 0.02) ~
      0 + rh +
      mgp(pm10, a = 2, knots = seq(8, 100, by = 4), ref_value = "left") +
      iwp(tmean, p = 3, knots = seq(-2, 30, by = 2), ref_value = 20, init = 0.001),
    data = london,
    config = list(num_threads = 1L, num_shards = 4L),
    control = list(maxit = 15L),
    verbose = FALSE,
    hessian = FALSE
  )
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
})

test_that("replace_vars and sim_random work with mgp", {
  skip_on_cran()
  set.seed(2)
  dat <- data.frame(
    x = seq(1, 5, length.out = 40),
    y = sqrt(seq(1, 5, length.out = 40)) + rnorm(40, sd = 0.1)
  )
  terms <- mgp("x", a = 2, knots = seq(1, 5, length.out = 8))
  terms2 <- replace_vars(terms, c(x = "z"))
  expect_identical(terms2$z_mgp@name, "z")
  expect_identical(names(terms2), c("z_mgp", "z_mgp_boundary"))

  fit <- adlaplace(
    normal(y) ~ 0 + mgp(x, a = 2, knots = seq(1, 5, length.out = 8)),
    data = dat,
    config = list(num_threads = 1L, num_shards = 1L),
    control = list(maxit = 25L),
    verbose = FALSE,
    hessian = TRUE
  )
  sims <- sim_random("x", fit = fit, num_sim = 10L, num_grid = 15L)
  expect_equal(ncol(sims$sim), 10L)
  expect_true(all(is.finite(sims$sim)))
})
