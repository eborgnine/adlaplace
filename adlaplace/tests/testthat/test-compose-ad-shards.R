test_that("combine matches expected shard count for small GLMM", {
  set.seed(2)
  Nobs <- 40L
  Nrandom1 <- 3L
  X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
  Amat <- Matrix::sparseMatrix(
    i = seq_len(Nobs),
    j = sample(Nrandom1, Nobs, replace = TRUE),
    x = 1
  )
  config <- list(
    beta = rep(0, 2),
    theta = -1,
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    groups = adlaplace::adFun_groups(as(Matrix::t(Amat), "dMatrix"), Ngroups = 10L),
    num_threads = 1L,
    verbose = FALSE
  )
  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  data_obs <- list(
    y = rpois(Nobs, 2),
    ATp = as(Matrix::t(Amat), "dMatrix"),
    XTp = as(Matrix::t(X), "CsparseMatrix"),
    theta_map = n_beta + n_gamma
  )
  data_r <- list(
    Q = rep(1, ncol(Amat)),
    theta_map = n_beta + n_gamma,
    gamma_map = seq.int(n_beta, length.out = ncol(Amat))
  )

  ad_ptr <- adlaplace::combine(
    list(
      adlaplace::get_ad_fun_raw(data_obs, config, kind = "obs", name = "neg_binom_obs"),
      adlaplace::get_ad_fun_raw(data_r, config, kind = "single", name = "random_diagonal"),
      adlaplace::get_ad_fun_raw(data_obs, config, kind = "single", name = "neg_binom_extra")
    ),
    config
  )
  expect_equal(
    adlaplace::n_groups(ad_ptr),
    ncol(config$groups) + 2L
  )
  x <- c(config$beta, config$gamma, config$theta)
  expect_true(is.finite(adlaplace::joint_log_dens(ad_ptr, x)))
})

test_that("get_ad_fun_raw obs-only builds observation groups only", {
  set.seed(3)
  Nobs <- 20L
  X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
  Amat <- Matrix::sparseMatrix(
    i = seq_len(Nobs),
    j = sample(3L, Nobs, replace = TRUE),
    x = 1
  )
  config <- list(
    beta = rep(0, 2),
    theta = -1,
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    groups = adlaplace::adFun_groups(as(Matrix::t(Amat), "dMatrix"), Ngroups = 5L),
    verbose = FALSE
  )
  data_obs <- list(
    y = rpois(Nobs, 2),
    ATp = as(Matrix::t(Amat), "dMatrix"),
    XTp = as(Matrix::t(X), "CsparseMatrix"),
    theta_map = length(config$beta) + length(config$gamma)
  )
  ad_obs <- adlaplace::get_ad_fun_raw(data_obs, config, kind = "obs", name = "neg_binom_obs")

  expect_error(
    adlaplace::get_ad_fun_raw(data_obs, config, kind = "obs"),
    "`name` is required",
    fixed = TRUE
  )
  expect_equal(adlaplace::n_groups(ad_obs), ncol(config$groups))
})
