test_that("get_ad_fun and derivatives run on small GLMM data", {
  set.seed(0)
  Nobs <- 80L
  Nrandom1 <- 4L
  Nrandom2 <- 5L
  X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
  AmatList <- list(
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom1, Nobs, replace = TRUE), x = 1),
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom2, Nobs, replace = TRUE), x = 1)
  )
  Amat <- do.call(cbind, AmatList)
  config <- list(
    beta = rep(0, ncol(X)),
    theta = c(-1, -1, -1),
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    groups = adlaplace::adFun_groups(as(Matrix::t(Amat), "dMatrix"), Ngroups = 20L),
    num_threads = 1L,
    verbose = FALSE,
    package = "adlaplace"
  )
  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  data_obs <- list(
    y = rpois(Nobs, 2),
    ATp = as(Matrix::t(Amat), "dMatrix"),
    XTp = as(Matrix::t(X), "CsparseMatrix"),
    theta_map = n_beta + n_gamma + length(config$theta) - 1L
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
  ad_fun <- adlaplace::get_ad_fun(data_obs, config, ad_ptr)
  x <- c(config$beta, config$gamma, config$theta)

  dens <- adlaplace::joint_log_dens(ad_fun, x)
  expect_true(is.finite(dens))

  g <- adlaplace::grad(ad_fun, x)
  expect_equal(length(g), length(x))

  H <- adlaplace::hessian(ad_fun, x)
  expect_true(inherits(H, "Matrix"))

  inner_res <- adlaplace::inner_opt(
    parameters = c(config$beta, config$theta),
    gamma = config$gamma,
    config = config,
    ad_fun = ad_fun,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = FALSE
  )
  expect_true(is.finite(inner_res$log_lik))
  expect_equal(inner_res$neg_log_lik, -inner_res$log_lik)
})

test_that("model_setup builds data for iwp formula", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  n <- 30L
  df <- data.frame(
    y = rnorm(n),
    x = runif(n),
    id = rep(1:10, each = 3L)
  )
  model_stuff <- adlaplace::model_setup(
    data = df,
    formula = y ~ intercept() + linear(x),
    verbose = FALSE
  )
  expect_true(is.list(model_stuff$data))
  expect_true(nrow(model_stuff$data$XTp) >= 1L)
  expect_true(length(model_stuff$info$beta$init) >= 1L)
})
