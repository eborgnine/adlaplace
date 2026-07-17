test_that("repeated OpenMP log_lik sessions do not abort", {
  skip_if_not(adlaplace::has_openmp(), "OpenMP not available in this build")
  set.seed(1)
  Nobs <- 200L
  Nrandom1 <- 6L
  Nrandom2 <- 8L
  X <- Matrix::Matrix(cbind(1, stats::rbinom(Nobs, 1, 0.5)))
  Amat <- cbind(
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom1, Nobs, replace = TRUE)),
    Matrix::sparseMatrix(i = seq_len(Nobs), j = sample(Nrandom2, Nobs, replace = TRUE))
  )
  config <- list(
    beta = rep(0, ncol(X)),
    theta = c(-1, -1, -1),
    transform_theta = TRUE,
    gamma = rep(0, ncol(Amat)),
    shards = adlaplace::ad_shards(Amat, num_shards = 80L),
    verbose = FALSE,
    package = "adlaplace"
  )
  model <- test_ad_data(
    y = stats::rpois(Nobs, 2),
    A = Amat,
    X = X,
    config = config,
    theta_local_row = length(config$theta) - 1L
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
  ad_fun <- adlaplace::ad_fun(ad_ptr, num_threads = 2L)
  x <- c(config$beta, config$theta)
  for (i in seq_len(25L)) {
    ll <- adlaplace::log_lik_laplace(
      x = x,
      config = list(verbose = FALSE),
      gamma = config$gamma,
      ad_fun = ad_fun,
      control = list(maxit = 4L, report.level = 0, report.freq = 0),
      deriv = TRUE
    )
    expect_true(is.finite(ll$log_lik), info = paste("iteration", i))
  }
})

test_that("OpenMP sessions work in a fresh Rscript process", {
  skip_if_not(adlaplace::has_openmp(), "OpenMP not available in this build")
  skip_on_cran()

  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(
    c(
      "library(adlaplace)",
      "stopifnot(isTRUE(adlaplace::has_openmp()))",
      "set.seed(2)",
      "n <- 60L",
      "dat <- data.frame(",
      "  y = rpois(n, 2),",
      "  x = runif(n),",
      "  g = factor(rep(seq_len(6L), length.out = n))",
      ")",
      "md <- adlaplace::model_data(",
      "  data = dat,",
      "  formula = adlaplace::nbinom(y, init = 0.2) ~ x + adlaplace::iid(g, init = 0.2)",
      ")",
      "config <- list(",
      "  transform_theta = TRUE,",
      "  shards = adlaplace::ad_shards(md$data$A, num_shards = 8L),",
      "  verbose = FALSE",
      ")",
      "af <- adlaplace::ad_fun(md, config, num_threads = 2L)",
      "for (i in 1:8) {",
      "  ll <- adlaplace::log_lik_laplace(",
      "    x = md$data$info$parameters$init,",
      "    ad_fun = af,",
      "    config = config,",
      "    control = list(maxit = 3L, report.level = 0, report.freq = 0),",
      "    deriv = TRUE",
      "  )",
      "  stopifnot(is.finite(ll$log_lik))",
      "}",
      "quit(save = 'no', status = 0)"
    ),
    script
  )

  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("--vanilla", script),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  expect_equal(status, 0L, info = paste(out, collapse = "\n"))
})
