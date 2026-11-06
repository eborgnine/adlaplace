test_that("gamma_obs + gamma_extra match dgamma for two groups", {
  set.seed(3)
  n <- 40L
  n_re <- 3L
  X <- Matrix::Matrix(cbind(1, rnorm(n)))
  A <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = sample(n_re, n, replace = TRUE),
    x = 1,
    dims = c(n, n_re)
  )
  beta <- c(2.5, -0.2)
  gamma <- rnorm(n_re, sd = 0.15)
  cv <- c(0.08, 0.12)
  group <- sample(c(0L, 1L), n, replace = TRUE)
  eta <- as.vector(X %*% beta + A %*% gamma)
  shape <- 1 / cv[group + 1L]^2
  y <- rgamma(n, shape = shape, rate = shape / exp(eta))

  config <- list(
    beta = beta,
    gamma = gamma,
    theta = log(cv),
    transform_theta = c(TRUE, TRUE),
    obs_groups = adlaplace::obs_groups(A, num_shards = 2L)
  )
  obs <- adlaplace:::density_data(
    y = y, A = A, X = X,
    theta_map = Matrix::Diagonal(2L),
    weights = as.numeric(group),
    ad_kind = "observations",
    density = "gamma_obs"
  )
  extra <- adlaplace::density_data(
    obs, ad_kind = "parameters", density = "gamma_extra"
  )
  ptr_obs <- adlaplace::ad_pack_ptr(obs, config)
  ptr_extra <- adlaplace::ad_pack_ptr(extra, config)

  x <- c(beta, gamma, log(cv))
  manual <- sum(stats::dgamma(
    y, shape = shape, rate = shape / exp(eta), log = TRUE
  ))
  ad_total <- adlaplace::joint_log_dens(ptr_obs, x, negative = FALSE) +
    adlaplace::joint_log_dens(ptr_extra, x, negative = FALSE)
  expect_equal(ad_total, manual, tolerance = 1e-8)
})

test_that("biv_iid log density matches a block-diagonal bivariate normal", {
  set.seed(4)
  n_station <- 5L
  s1 <- 0.4
  s2 <- 0.15
  rho <- -0.35
  g <- rnorm(n_station * 2L)
  config <- list(
    beta = numeric(0),
    gamma = g,
    theta = c(log(s1), log(s2), rho),
    transform_theta = c(TRUE, TRUE, FALSE)
  )
  model <- adlaplace:::density_data(
    y = rep(1, length(g)),
    A = Matrix::Diagonal(length(g)),
    theta_map = Matrix::Diagonal(3L),
    precision = list(n_pollutant = 2L),
    ad_kind = "random",
    density = "biv_iid"
  )
  ptr <- adlaplace::ad_pack_ptr(model, config)
  manual <- 0
  for (s in seq_len(n_station)) {
    g1 <- g[2L * s - 1L]
    g2 <- g[2L * s]
    quad <- (g1^2 / s1^2 + g2^2 / s2^2 - 2 * rho * g1 * g2 / (s1 * s2)) /
      (1 - rho^2)
    manual <- manual - 0.5 * quad - log(s1) - log(s2) -
      0.5 * log(1 - rho^2) - log(2 * pi)
  }
  expect_equal(
    adlaplace::joint_log_dens(
      ptr, c(g, log(s1), log(s2), rho), negative = FALSE
    ),
    manual,
    tolerance = 1e-8
  )
})

test_that("bivariate gamma fit recovers the sign of the correlation", {
  set.seed(5)
  n_s <- 30L
  n_t <- 20L
  s1 <- 0.45
  s2 <- 0.2
  rho <- 0.7
  sigma <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2L, 2L)
  b <- matrix(rnorm(n_s * 2L), n_s, 2L) %*% chol(sigma)
  dat <- expand.grid(
    station = factor(seq_len(n_s)),
    day = seq_len(n_t),
    pollutant = factor(c("PM2.5", "PM10"), levels = c("PM2.5", "PM10")),
    KEEP.OUT.ATTRS = FALSE
  )
  beta <- c(2.6, 3.3)
  cv <- 0.07
  pol <- as.integer(dat$pollutant)
  eta <- beta[pol] + b[cbind(as.integer(dat$station), pol)]
  shape <- 1 / cv^2
  dat$value <- rgamma(nrow(dat), shape = shape, rate = shape / exp(eta))
  dat$pm25_ind <- as.numeric(dat$pollutant == "PM2.5")
  dat$pm10_ind <- as.numeric(dat$pollutant == "PM10")

  fit <- adlaplace::adlaplace(
    gamma(value, by = "pollutant") ~ 0 + pm25_ind + pm10_ind +
      biv_iid(station, by = "pollutant", init = 0.3, rho_init = 0.2),
    data = dat,
    config = list(num_threads = 1L, num_shards = 4L),
    control = list(maxit = 80L),
    hessian = FALSE
  )
  expect_identical(fit$details$outer_opt$convergence, 0L)
  co <- coef(fit)
  expect_gt(unname(co[["station_biv_rho"]]), 0.2)
  expect_equal(unname(co[["pm25_ind_linear_pm25_ind"]]), 2.6, tolerance = 0.15)
  expect_equal(unname(co[["pm10_ind_linear_pm10_ind"]]), 3.3, tolerance = 0.15)
  expect_equal(unname(co[["value_gamma_cv_PM2.5"]]), cv, tolerance = 0.04)
  expect_equal(unname(co[["value_gamma_cv_PM10"]]), cv, tolerance = 0.04)
  expect_lt(abs(unname(co[["station_biv_sd_PM2.5"]]) - s1), 0.15)
  expect_lt(abs(unname(co[["station_biv_sd_PM10"]]) - s2), 0.1)
})
