# Mixed-DSO skew-normal: obs/extra in adlaplaceExample.so, random in adlaplace.so.

build_skewnormal_ad_fun <- function(num_threads = 4L, num_shards = 20L) {
  skip_if_not_installed("sn")

  set.seed(42)
  Nobs <- 120L
  Nrandom1 <- 4L
  Nrandom2 <- 5L

  X <- Matrix::Matrix(cbind(1, stats::rbinom(Nobs, 1, 0.5)))
  Amat <- cbind(
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom1, Nobs, replace = TRUE)),
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom2, Nobs, replace = TRUE))
  )

  beta <- rep(0, ncol(X))
  thetaOrig <- c(sd1 = 0.1, sd2 = 0.1, omega = 0.7, alpha = 0.8)
  gamma <- rep(0, ncol(Amat))
  eta_true <- as.vector(X %*% beta + Amat %*% gamma)
  y <- sn::rsn(Nobs, xi = eta_true, omega = thetaOrig["omega"], alpha = thetaOrig["alpha"])

  config <- list(
    beta = beta,
    theta = c(log(thetaOrig[c("sd1", "sd2", "omega")]), thetaOrig["alpha"]),
    transform_theta = TRUE,
    gamma = gamma,
    shards = adlaplace::ad_shards(Amat, num_shards = num_shards),
    num_threads = as.integer(num_threads),
    verbose = FALSE
  )

  n_beta <- length(beta)
  n_gamma <- length(gamma)
  n_theta <- length(config$theta)

  obs_theta_idx <- (n_theta - 1L):n_theta
  rand_theta_idx <- seq_len(n_theta - 2L)

  obs_shard <- adlaplace::ad_data(
    y = y,
    A = Amat,
    X = X,
    beta_map = Matrix::Diagonal(n_beta),
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(obs_theta_idx, n_theta),
    ad_kind = "observations",
    ad_fun = "skewnormal_obs",
    package = "adlaplaceExample"
  )
  random_shard <- adlaplace::ad_data(
    beta_map = n_beta,
    gamma_map = Matrix::Diagonal(n_gamma),
    theta_map = list(rand_theta_idx, n_theta),
    ad_kind = "random",
    ad_fun = "random_diagonal",
    precision = rep(1, n_gamma)
  )
  extra_shard <- adlaplace::ad_data(
    y = y,
    beta_map = n_beta,
    gamma_map = n_gamma,
    theta_map = list(obs_theta_idx, n_theta),
    ad_kind = "parameters",
    ad_fun = "skewnormal_extra",
    package = "adlaplaceExample"
  )

  ptrs <- c(
    adlaplace::ad_fun_ptr(obs_shard, config),
    adlaplace::ad_fun_ptr(random_shard, config),
    adlaplace::ad_fun_ptr(extra_shard, config)
  )

  list(
    config = config,
    ad_fun = adlaplace::ad_fun(ptrs, num_threads = num_threads)
  )
}

skip_parallel_on_macos <- function() {
  skip_if(
    Sys.info()[["sysname"]] == "Darwin" &&
      Sys.getenv("ADLAPLACE_TEST_PARALLEL_INNER_OPT", "") != "1",
    "parallel deriv+trace test (set ADLAPLACE_TEST_PARALLEL_INNER_OPT=1 on macOS)"
  )
}

test_that("serial log_lik_laplace deriv=TRUE with mixed-DSO skewnormal", {
  skip_if_not_installed("sn")
  fx <- build_skewnormal_ad_fun(num_threads = 1L, num_shards = 20L)
  ll <- adlaplace::log_lik_laplace(
    x = c(fx$config$beta, fx$config$theta),
    config = list(verbose = FALSE),
    gamma = fx$config$gamma,
    ad_fun = fx$ad_fun,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  expect_true(is.finite(ll$log_lik))
  expect_true(all(is.finite(ll$deriv$d_neg_log_lik)))
  expect_true(all(is.finite(ll$extra$trace3)))
})

test_that("parallel log_lik_laplace deriv=TRUE with mixed-DSO skewnormal", {
  skip_if_not_installed("sn")
  skip_parallel_on_macos()
  fx <- build_skewnormal_ad_fun(num_threads = 4L, num_shards = 20L)
  ll <- adlaplace::log_lik_laplace(
    x = c(fx$config$beta, fx$config$theta),
    config = list(verbose = FALSE),
    gamma = fx$config$gamma,
    ad_fun = fx$ad_fun,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  expect_true(is.finite(ll$log_lik))
  expect_true(all(is.finite(ll$deriv$d_neg_log_lik)))
  expect_true(all(is.finite(ll$extra$trace3)))
})

test_that("back-to-back parallel deriv=TRUE reuses trace prep on mixed-DSO skewnormal", {
  skip_if_not_installed("sn")
  skip_parallel_on_macos()
  fx <- build_skewnormal_ad_fun(num_threads = 4L, num_shards = 20L)
  args <- list(
    x = c(fx$config$beta, fx$config$theta),
    config = list(verbose = FALSE),
    gamma = fx$config$gamma,
    ad_fun = fx$ad_fun,
    control = list(maxit = 4L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  ll1 <- do.call(adlaplace::log_lik_laplace, args)
  ll2 <- do.call(adlaplace::log_lik_laplace, args)
  expect_true(is.finite(ll1$log_lik))
  expect_true(is.finite(ll2$log_lik))
  expect_true(all(is.finite(ll2$extra$trace3)))
})
