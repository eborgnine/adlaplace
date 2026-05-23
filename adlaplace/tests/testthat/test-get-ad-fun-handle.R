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
  n_beta <- length(config$beta)
  model <- test_ad_model(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  random_shard <- test_random_shard(
    model = model,
    config = config,
    gamma_ids = seq.int(0L, length.out = ncol(Amat)),
    theta_id = 0L,
    Q = rep(1, ncol(Amat))
  )
  ad_ptr <- do.call(c, list(
    adlaplace::ad_fun_ptr(config, "observations", "neg_binom_obs", model),
    adlaplace::ad_fun_ptr(
      config, "random", "random_diagonal",
      random_shard$model, precision = random_shard$precision
    ),
    adlaplace::ad_fun_ptr(
      config, "parameters", "neg_binom_extra", model
    )
  ))
  af <- adlaplace::ad_fun(ad_ptr)
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
    config = config, ad_fun = af,
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
  model <- test_ad_model(
    y = rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config
  )
  ad_ptr <- adlaplace::ad_fun_ptr(config, "observations", "neg_binom_obs", model)
  af <- adlaplace::ad_fun(ad_ptr)
  sz <- adlaplace::get_sizes(ad_ptr, 0L)
  expect_true(methods::is(af, "ad_fun"))
  expect_equal(af@sizes[["beta"]], sz$n_beta)
  expect_equal(af@sizes[["gamma"]], sz$n_outer - sz$n_beta - sz$n_theta)
  expect_equal(af@sizes[["theta"]], sz$n_theta)
})
