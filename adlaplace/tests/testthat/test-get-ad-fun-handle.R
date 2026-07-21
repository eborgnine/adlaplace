test_that("ad_fun attaches hessian_map to combined handle", {
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
  x <- c(config$beta, config$gamma, config$theta)
  expect_equal(
    adlaplace::joint_log_dens(af, x, negative = FALSE),
    adlaplace::joint_log_dens(ad_ptr, x, negative = FALSE)
  )
  expect_equal(
    adlaplace::joint_log_dens(ad_ptr, x, negative = FALSE),
    -adlaplace::joint_log_dens(ad_ptr, x, negative = TRUE)
  )
  inner_res <- adlaplace::inner_opt(
    c(config$beta, config$theta), config$gamma,
    ad_fun = af,
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
  expect_true(is.list(af@chol_inner_list))
  expect_true(inherits(af@chol_inner_list$L1, "dtCMatrix"))

  ci <- inner_res$hessian$chol_inner
  cil <- af@chol_inner_list
  expect_equal(ci$L1@i, cil$L1@i)
  expect_equal(ci$L1@p, cil$L1@p)
  expect_equal(ci$perm, cil$perm)
  expect_equal(ci$perm_inv, cil$perm_inv)

  H_perm <- ci$L1 %*% Matrix::Diagonal(length(ci$D), ci$D) %*% Matrix::t(ci$L1)
  reconstructH <- H_perm[order(ci$perm), order(ci$perm)]
  expect_true(max(abs(reconstructH - inner_res$hessian$inner)) < 1e-8)
})

test_that("inner_opt deriv exports Linv matching chol_inner_list pattern", {
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
  cil <- af@chol_inner_list
  expect_true(!is.null(ci$Linv))
  expect_equal(ci$Linv@i, cil$Linv@i)
  expect_equal(ci$Linv@p, cil$Linv@p)
  Linv_ref <- Matrix::solve(ci$L1)
  expect_true(max(abs(ci$Linv - Linv_ref)) < 1e-8)

  laplace <- adlaplace::log_lik_laplace(
    x = c(config$beta, config$theta),
    config = list(verbose = FALSE),
    gamma = config$gamma,
    ad_fun = af,
    control = list(maxit = 3L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  expect_true(!is.null(laplace$hessian$chol_inner$Linv))
  expect_true(!is.null(laplace$hessian$half_H_inv))
})

test_that("ad_fun works without config when layout is on ptr", {
  set.seed(5)
  Nobs <- 25L
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
    shards = adlaplace::ad_shards(Amat, num_shards = 6L),
    verbose = FALSE
  )
  model <- test_ad_data(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  ad_ptr <- adlaplace::ad_fun_ptr(
    as_shard(model, "observations", "nbinom_obs"),
    config
  )
  af <- adlaplace::ad_fun(ad_ptr, num_threads = 1L)
  sz <- adlaplace:::get_sizes(ad_ptr, 0L)
  expect_true(methods::is(af, "ad_fun"))
  expect_equal(af@sizes[["beta"]], sz$n_beta)
  expect_equal(af@sizes[["gamma"]], sz$n_outer - sz$n_beta - sz$n_theta)
  expect_equal(af@sizes[["theta"]], sz$n_theta)
})

test_that("ad_fun variadic ad_fun_ptr matches explicit c() composition", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(14)
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
    num_threads = 2L,
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

  ad_ptr_explicit <- do.call(c, list(
    adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config),
    adlaplace::ad_fun_ptr(random_shard, config),
    adlaplace::ad_fun_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  ))
  dens_explicit <- adlaplace::clone_ad_fun_ptr(ad_ptr_explicit)
  af_explicit <- adlaplace::ad_fun(ad_ptr_explicit, num_threads = 2L)

  make_three <- function() {
    list(
      adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config),
      adlaplace::ad_fun_ptr(random_shard, config),
      adlaplace::ad_fun_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
    )
  }
  dens_variadic <- do.call(c, make_three())
  af_variadic <- do.call(
    adlaplace::ad_fun,
    c(make_three(), list(num_threads = 2L))
  )

  x <- c(config$beta, config$gamma, config$theta)
  expect_equal(
    adlaplace::joint_log_dens(dens_explicit, x, negative = FALSE),
    adlaplace::joint_log_dens(dens_variadic, x, negative = FALSE)
  )
  expect_equal(
    as.matrix(af_explicit@parallel_map),
    as.matrix(af_variadic@parallel_map)
  )
  expect_equal(af_explicit@sizes, af_variadic@sizes)
})

test_that("ad_fun variadic composition clears source pointers like c()", {
  skip_if_not(adlaplace:::has_openmp(), "OpenMP not available in this build")
  set.seed(15)
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
    shards = adlaplace::ad_shards(Amat, num_shards = 5L),
    num_threads = 2L,
    verbose = FALSE
  )
  model <- test_ad_data(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )

  obs <- adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config)
  extra <- adlaplace::ad_fun_ptr(as_shard(model, "parameters", "nbinom_extra"), config)
  dens_ptr <- c(
    adlaplace::clone_ad_fun_ptr(obs),
    adlaplace::clone_ad_fun_ptr(extra)
  )
  af <- adlaplace::ad_fun(obs, extra, num_threads = 2L)
  x <- c(config$beta, config$gamma, config$theta)

  expect_true(is(af, "ad_fun"))
  expect_true(is.finite(adlaplace::joint_log_dens(dens_ptr, x, negative = FALSE)))
  expect_error(
    adlaplace::joint_log_dens(af, x, negative = FALSE),
    "serial debug APIs|num_threads"
  )
  expect_error(
    adlaplace::joint_log_dens(obs, x, negative = FALSE),
    "NULL|cleared|invalid",
    ignore.case = TRUE
  )
  expect_error(
    adlaplace::joint_log_dens(extra, x, negative = FALSE),
    "NULL|cleared|invalid",
    ignore.case = TRUE
  )
})

test_that("ad_fun variadic input validates additional arguments", {
  set.seed(16)
  Nobs <- 10L
  X <- Matrix::Matrix(cbind(1, rbinom(Nobs, 1, prob = 0.5)))
  Amat <- Matrix::sparseMatrix(
    i = seq_len(Nobs),
    j = sample(2L, Nobs, replace = TRUE),
    x = 1
  )
  config <- list(
    beta = rep(0, 2),
    theta = -1,
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    shards = adlaplace::ad_shards(Amat, num_shards = 3L),
    verbose = FALSE
  )
  model <- test_ad_data(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  obs <- adlaplace::ad_fun_ptr(as_shard(model, "observations", "nbinom_obs"), config)

  expect_error(
    adlaplace::ad_fun(obs, 1L),
    "additional arguments must be ad_fun_ptr",
    fixed = TRUE
  )
})
