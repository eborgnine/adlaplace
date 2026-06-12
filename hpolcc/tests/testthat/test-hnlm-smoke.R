test_that("hnlm for_dev builds ad_fun via ad_fun", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_shards = 2L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )
  expect_s3_class(forres, "hnlm")
  expect_s3_class(forres, "hnlm_dev")
  expect_true(methods::is(forres$ad_fun, "ad_fun"))
  expect_true(is.list(forres$model_data))
  expect_true("observations" %in% names(forres$model_data))

  sz <- forres$ad_fun@sizes
  x <- c(
    forres$config$opt$init[seq_len(sz["beta"])],
    rep(0, sz["gamma"]),
    forres$config$opt$init[seq_len(sz["theta"]) + sz["beta"]]
  )
  dens <- adlaplace::joint_log_dens(forres$ad_fun, x)
  expect_true(is.finite(dens))
})

test_that("hnlm fit returns flat hnlm object", {
  skip_if_not_installed("adlaplace")
  skip_if_not_installed("adlaplaceHgp")
  td <- make_hpolcc_test_data()
  fit <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_shards = 4L,
      num_threads = 1L
    ),
    control = list(maxit = 3L, trace = 0),
    control_inner = list(report.level = 0, maxit = 5L)
  )
  expect_s3_class(fit, "hnlm")
  expect_false("hnlm_dev" %in% class(fit))
  expect_null(fit$objects)
  expect_true(is.list(fit$coefficients))
  expect_true(is.finite(fit$log_lik))
  expect_true(is.list(fit$extra))
  expect_true(is.list(fit$info))
  expect_true(is.list(fit$call))
  expect_true(is.list(fit$hessian))
  expect_true(is.list(fit$optim))
  expect_true(is.logical(fit$converged))
})
