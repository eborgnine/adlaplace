test_that("binomial_obs with size matches softplus kernel", {
  set.seed(11)
  n <- 30L
  X <- Matrix::Matrix(cbind(1, rnorm(n)))
  beta <- c(-0.2, 0.4)
  eta <- as.vector(X %*% beta)
  N <- sample(5:20, n, replace = TRUE)
  prob <- 1 / (1 + exp(-eta))
  y <- rbinom(n, N, prob)

  config <- list(
    beta = beta, gamma = numeric(0), theta = numeric(0),
    transform_theta = TRUE
  )
  model <- adlaplace:::ad_data(
    y = y, X = X,
    theta_map = Matrix::Matrix(nrow = 0L, ncol = 0L),
    ad_kind = "observations",
    ad_fun = "binomial_obs",
    weights = as.numeric(N)
  )
  ptr <- adlaplace::ad_fun_ptr(model, config)
  x <- beta
  # dbinom includes choose(N,y); AD density omits that constant
  manual_kernel <- sum(y * eta - N * log1p(exp(eta)))
  expect_equal(
    adlaplace::joint_log_dens(ptr, x, negative = FALSE),
    manual_kernel,
    tolerance = 1e-10
  )
})

test_that("binomial(y, size = N) wires weights through model_data", {
  dat <- data.frame(
    y = c(1L, 2L, 0L),
    N = c(5L, 8L, 3L),
    x1 = c(0.1, -0.2, 0.3)
  )
  md <- adlaplace::model_data(
    binomial(y, size = N) ~ x1,
    data = dat
  )
  expect_equal(md$observations$y@weights, as.numeric(dat$N))
  term <- md$terms[[grep("binomial", names(md$terms))[1L]]]
  expect_s4_class(term, "binomial")
  expect_identical(term@size, "N")
})

test_that("collect_terms evaluates in formula env and coerces size", {
  f <- binomial(y, size = N) ~ x1
  terms <- adlaplace::collect_terms(f)
  bin <- terms[[grep("binomial", names(terms))[1L]]]
  expect_s4_class(bin, "binomial")
  expect_identical(bin@term, "y")
  expect_identical(bin@size, "N")
})

test_that("model_data passes list precision through for random_fem", {
  skip_if_not_installed("adlaplaceGrf")
  knots_list <- list(
    x = seq(0, 1, length.out = 4),
    y = seq(0, 1, length.out = 4)
  )
  g <- do.call(expand.grid, knots_list)
  fem <- adlaplaceGrf::grf_bspline(g, knots_list, degree = 2L)
  prec <- adlaplaceGrf::fem_precision_payload(fem, alpha = 2L)
  nr <- nrow(fem$C)

  # Simulate what matern precision() returns via a tiny hand-built random term path:
  # build ad_data directly and ensure list precision is accepted by the AD kernel.
  model <- adlaplace:::ad_data(
    gamma_map = Matrix::Diagonal(nr),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "random",
    ad_fun = "random_fem_2",
    package = "adlaplaceGrf",
    precision = prec
  )
  expect_true(is.list(model@precision))
  config <- list(
    gamma = rep(0, nr),
    theta = c(0, 0),
    transform_theta = TRUE
  )
  ptr <- adlaplace::ad_fun_ptr(model, config)
  expect_true(is.finite(adlaplace::joint_log_dens(ptr, c(rep(0, nr), 0, 0), negative = FALSE)))
})
