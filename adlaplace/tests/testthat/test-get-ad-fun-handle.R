test_that("get_ad_fun attaches hessian_map to combined handle", {
  set.seed(4)
  Nobs <- 30L
  X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
  Amat <- Matrix::sparseMatrix(
    i = seq_len(Nobs),
    j = sample(4L, Nobs, replace = TRUE),
    x = 1
  )
  config <- list(
    beta = rep(0, 2),
    theta = -1,
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    groups = adlaplace::adFun_groups(as(Matrix::t(Amat), "dMatrix"), Ngroups = 8L),
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
  ad_backend <- adlaplace::get_ad_fun(data_obs, config, ad_ptr)
  x <- c(config$beta, config$gamma, config$theta)
  expect_equal(
    adlaplace::joint_log_dens(ad_backend, x, negative = FALSE),
    adlaplace::joint_log_dens(ad_ptr, x, negative = FALSE)
  )
  expect_equal(
    adlaplace::joint_log_dens(ad_ptr, x, negative = FALSE),
    -adlaplace::joint_log_dens(ad_ptr, x, negative = TRUE)
  )
  inner_res <- adlaplace::inner_opt(
    c(config$beta, config$theta), config$gamma,
    config = config, ad_fun = ad_backend,
    control = list(maxit = 3L, report.level = 0, report.freq = 0),
    deriv = FALSE
  )
  expect_true(is.finite(inner_res$log_lik))
  expect_equal(inner_res$neg_log_lik, -inner_res$log_lik)
  expect_length(inner_res$gradient$inner, length(config$gamma))
  expect_length(inner_res$gradient$outer, 0L)
  expect_equal(nrow(inner_res$hessian$inner), length(config$gamma))
  expect_equal(nrow(inner_res$hessian$outer), 0L)
  expect_true(is.list(inner_res$hessian$chol_inner))
  expect_true(is.list(ad_backend$chol_inner_list))
  expect_true(inherits(ad_backend$chol_inner_list$L1, "dtCMatrix"))
})
