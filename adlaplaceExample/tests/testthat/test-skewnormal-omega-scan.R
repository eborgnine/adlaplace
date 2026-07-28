test_that("skewnormal obs and extra log densities vary with log(omega)", {
  skip_if_not_installed("sn")

  set.seed(0L)
  Nobs <- 200L
  Nrandom1 <- 5L
  Nrandom2 <- 8L

  X <- Matrix::Matrix(cbind(1, stats::rbinom(Nobs, 1, 0.5)))
  Adf <- data.frame(
    r1 = sample(Nrandom1, Nobs, replace = TRUE),
    r2 = sample(Nrandom2, Nobs, replace = TRUE)
  )
  Amat <- cbind(
    Matrix::sparseMatrix(i = seq_len(Nobs), j = Adf$r1),
    Matrix::sparseMatrix(i = seq_len(Nobs), j = Adf$r2)
  )

  beta <- rep(1, ncol(X))
  thetaOrig <- c(sd1 = 0.1, sd2 = 0.1, omega = 0.7, alpha = 0.8)
  gamma <- rnorm(ncol(Amat), sd = rep(thetaOrig[1:2], c(Nrandom1, Nrandom2)))
  eta_true <- as.vector(X %*% beta + Amat %*% gamma)
  y <- sn::rsn(Nobs, xi = eta_true, omega = thetaOrig["omega"], alpha = thetaOrig["alpha"])

  config <- list(
    beta = beta,
    theta = c(log(thetaOrig[c("sd1", "sd2", "omega")]), thetaOrig["alpha"]),
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    obs_groups = adlaplace::obs_groups(Amat, num_shards = 10L),
    num_threads = 1L,
    verbose = FALSE
  )

  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  n_theta <- length(config$theta)

  data_obs <- adlaplace::density_data(
    y = y,
    A = Amat,
    X = X,
    beta_map = Matrix::Diagonal(n_beta),
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(c(
      which(names(thetaOrig) == "omega"),
      which(names(thetaOrig) == "alpha")
    ), n_theta),
    ad_kind = "observations",
    density = "skewnormal_obs",
    package = "adlaplaceExample"
  )
  data_extra <- adlaplace::density_data(
    y = data_obs@y,
    beta_map = data_obs@beta_map,
    gamma_map = data_obs@gamma_map,
    theta_map = data_obs@theta_map,
    ad_kind = "parameters",
    density = "skewnormal_extra",
    package = "adlaplaceExample"
  )

  ad_fun_obs <- adlaplace::ad_pack_ptr(data = data_obs, config = config)
  ad_fun_extra <- adlaplace::ad_pack_ptr(data = data_extra, config = config)

  xx <- c(beta, config$gamma, config$theta)
  omega_idx <- n_beta + n_gamma + which(names(thetaOrig) == "omega")
  log_omega_true <- config$theta[which(names(thetaOrig) == "omega")]
  seq_log_omega <- log_omega_true + c(-0.2, 0, 0.2)

  eta_ad <- as.vector(X %*% beta + Amat %*% config$gamma)

  sn_obs_ref <- function(log_omega) {
    omega <- exp(log_omega)
    z <- (y - eta_ad) / (omega * sqrt(2))
    tt <- thetaOrig["alpha"] * z
    sum(-z^2 + log(pracma::erfc(-tt)))
  }

  ad_obs <- vapply(seq_log_omega, function(lo) {
    x <- xx
    x[omega_idx] <- lo
    adlaplace::joint_log_dens(ad_fun_obs, x, negative = FALSE)
  }, numeric(1))
  ref_obs <- vapply(seq_log_omega, sn_obs_ref, numeric(1))
  ad_extra <- vapply(seq_log_omega, function(lo) {
    x <- xx
    x[omega_idx] <- lo
    adlaplace::joint_log_dens(ad_fun_extra, x, negative = FALSE)
  }, numeric(1))

  expect_true(max(ad_obs) - min(ad_obs) > 1e-6)
  expect_true(max(ad_extra) - min(ad_extra) > 1e-6)
  expect_equal(ad_obs, ref_obs, tolerance = 1e-5)

  g_obs <- adlaplace::grad(ad_fun_obs, xx, negative = FALSE)
  expect_true(is.finite(g_obs[omega_idx]))
  expect_true(abs(g_obs[omega_idx]) > 1e-8)

  eps <- 1e-5
  x0 <- xx
  x0[omega_idx] <- log_omega_true
  fd <- (adlaplace::joint_log_dens(ad_fun_extra, {
    x <- x0
    x[omega_idx] <- log_omega_true + eps
    x
  }, negative = FALSE) -
    adlaplace::joint_log_dens(ad_fun_extra, {
      x <- x0
      x[omega_idx] <- log_omega_true - eps
      x
    }, negative = FALSE)) / (2 * eps)
  g <- adlaplace::grad(ad_fun_extra, x0, negative = FALSE)
  expect_equal(g[omega_idx], fd, tolerance = 1e-4)
})
