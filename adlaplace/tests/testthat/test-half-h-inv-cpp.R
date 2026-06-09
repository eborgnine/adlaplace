test_that("hessian_map chol_inner_list includes half_H_inv and H_inv patterns", {
  set.seed(2)
  Nobs <- 20L
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
    shards = adlaplace::ad_shards(Amat, num_shards = 4L),
    num_threads = 1L,
    verbose = FALSE
  )
  model <- test_ad_data(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  random_shard <- test_random_shard(
    data = model,
    config = config,
    gamma_ids = seq.int(0L, length.out = ncol(Amat)),
    theta_id = 0L,
    Q = rep(1, ncol(Amat))
  )
  ad_ptr <- do.call(c, list(
    adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_fun_ptr(random_shard, config),
    adlaplace::ad_fun_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  af <- adlaplace::ad_fun(ad_ptr, num_threads = 1L)
  cil <- af@chol_inner_list
  expect_true(!is.null(cil$half_H_inv))
  expect_true(!is.null(cil$H_inv))
  expect_equal(nrow(cil$half_H_inv), ncol(Amat))
  expect_equal(nrow(cil$H_inv), ncol(Amat))
})

test_that("inner_opt deriv half_H_inv and H_inv match R reference", {
  set.seed(41)
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
    shards = adlaplace::ad_shards(Amat, num_shards = 8L),
    num_threads = 1L,
    verbose = FALSE
  )
  model <- test_ad_data(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  random_shard <- test_random_shard(
    data = model,
    config = config,
    gamma_ids = seq.int(0L, length.out = ncol(Amat)),
    theta_id = 0L,
    Q = rep(1, ncol(Amat))
  )
  ad_ptr <- do.call(c, list(
    adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_fun_ptr(random_shard, config),
    adlaplace::ad_fun_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  af <- adlaplace::ad_fun(ad_ptr, num_threads = 1L)
  inner_deriv <- adlaplace::inner_opt(
    c(config$beta, config$theta), config$gamma,
    ad_fun = af,
    control = list(maxit = 3L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  ci <- inner_deriv$hessian$chol_inner
  expect_true(!is.null(inner_deriv$hessian$half_H_inv))
  expect_true(!is.null(inner_deriv$hessian$H_inv))

  half_ref <- adlaplace:::half_H_inv_from_ldl(adlaplace:::as_ldl_list(ci))
  expect_true(max(abs(inner_deriv$hessian$half_H_inv - half_ref)) < 1e-8)

  H_inv_ref <- Matrix::tcrossprod(half_ref)
  expect_true(max(abs(inner_deriv$hessian$H_inv - H_inv_ref)) < 1e-8)

  H_inner <- inner_deriv$hessian$inner
  expect_true(max(abs(inner_deriv$hessian$H_inv - Matrix::solve(H_inner))) < 1e-6)
})
