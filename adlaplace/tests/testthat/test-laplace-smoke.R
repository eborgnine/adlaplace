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
  map <- Matrix::sparseMatrix(
    i = seq(0L, len = ncol(Amat)),
    j = rep(seq(0L, len = length(AmatList)), vapply(AmatList, ncol, integer(1))),
    x = 1L,
    index1 = FALSE,
    dims = c(ncol(Amat), length(AmatList))
  )
  data <- list(
    y = rpois(Nobs, 2),
    ATp = as(Matrix::t(Amat), "dMatrix"),
    XTp = as(Matrix::t(X), "CsparseMatrix"),
    map = map,
    Qdiag = rep(1, ncol(Amat))
  )
  config <- list(
    beta = rep(0, nrow(data$XTp)),
    theta = c(-1, -1, -1),
    transform_theta = TRUE,
    gamma = rep(0, nrow(data$ATp)),
    groups = adlaplace::adFun_groups(data$ATp, Ngroups = 20L),
    num_threads = 1L,
    verbose = FALSE,
    package = "adlaplace"
  )

  ad_fun <- adlaplace::get_ad_fun(data, config)
  x <- c(config$beta, config$gamma, config$theta)

  dens <- adlaplace::jointLogDens(ad_fun, x)
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
