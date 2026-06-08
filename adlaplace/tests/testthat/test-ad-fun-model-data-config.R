test_that("ad_fun(model_data) keeps user-supplied config beta and theta", {
  skip_if_not_installed("adlaplaceExample")
  skip_if_not_installed("sn")

  set.seed(1L)
  Nobs <- 40L
  X <- Matrix::Matrix(cbind(1, stats::rnorm(Nobs)))
  r1 <- sample(3L, Nobs, replace = TRUE)
  Amat <- Matrix::sparseMatrix(
    i = seq_len(Nobs),
    j = r1,
    dims = c(Nobs, 3L)
  )
  beta <- c(0.1, -0.2)
  thetaOrig <- c(sd1 = 0.2, omega = 0.5, alpha = 0.3)
  y <- sn::rsn(Nobs, xi = as.vector(X %*% beta), omega = thetaOrig["omega"], alpha = thetaOrig["alpha"])

  dat <- data.frame(y = y, x1 = X[, 2], r1 = r1)
  md <- adlaplace::model_data(
    adlaplaceExample::skewnormal(y) ~ x1 + adlaplace::f(r1, model = "iid"),
    data = dat
  )

  custom_beta <- beta + c(0.5, -0.5)
  custom_theta <- c(log(0.15), log(0.6), 0.9)

  config_default <- list(
    transform_theta = TRUE,
    shards = adlaplace::ad_shards(md$data$A, num_shards = 4L),
    verbose = FALSE
  )
  config_custom <- modifyList(config_default, list(beta = custom_beta, theta = custom_theta))

  ad_default <- adlaplace::ad_fun(md, config_default, num_threads = 1L)
  ad_custom <- adlaplace::ad_fun(md, config_custom, num_threads = 1L)

  x_default <- c(md$data$info$beta$init, rep(0, nrow(md$data$info$gamma)), md$data$info$theta$init)
  x_custom <- c(custom_beta, rep(0, nrow(md$data$info$gamma)), custom_theta)

  ld_default <- adlaplace::joint_log_dens(ad_default, x_default, negative = FALSE)
  ld_custom <- adlaplace::joint_log_dens(ad_custom, x_custom, negative = FALSE)

  expect_true(is.finite(ld_default))
  expect_true(is.finite(ld_custom))
  expect_false(isTRUE(all.equal(ld_default, ld_custom, tolerance = 1e-8)))
})
