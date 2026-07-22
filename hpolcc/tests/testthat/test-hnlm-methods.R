test_that("hnlm S3 methods work on fitted object", {
  skip_if_not_installed("adlaplace")
  skip_if_not_installed("adlaplaceHgp")
  td <- make_hpolcc_test_data()
  fit <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_groups = 4L,
      num_threads = 1L
    ),
    control = list(maxit = 3L, trace = 0),
    control_inner = list(report.level = 0, maxit = 5L)
  )

  expect_output(print(fit), "hnlm fit")
  expect_true(is.finite(as.numeric(stats::logLik(fit))))
  expect_equal(
    length(stats::coef(fit)),
    nrow(fit$coefficients$parameters)
  )
  expect_true(is.finite(fit$log_lik))
  expect_equal(as.numeric(stats::logLik(fit)), fit$log_lik)

  sim <- stats::simulate(fit, nsim = 5L)
  expect_true(is.list(sim))
  expect_true("sim" %in% names(sim))
})

test_that("summary and vcov use outer Hessian", {
  skip_if_not_installed("adlaplace")
  skip_if_not_installed("adlaplaceHgp")
  td <- make_hpolcc_test_data()
  fit <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_groups = 4L,
      num_threads = 1L
    ),
    control = list(maxit = 5L, trace = 0),
    control_inner = list(report.level = 0, maxit = 8L)
  )

  V <- stats::vcov(fit)
  expect_equal(dim(V), c(length(fit$optim$par), length(fit$optim$par)))
  expect_true(all(is.finite(V)))
  expect_equal(V, t(V), tolerance = 1e-8)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.hnlm")
  expect_equal(nrow(sm$coefficients), nrow(fit$coefficients$parameters))
  expect_true("mle" %in% names(sm$coefficients))
  expect_true("se" %in% names(sm$coefficients))
  # Natural-scale SEs for log-transformed parameters are left as NA.
  non_log <- !(sm$coefficients$log %in% TRUE)
  expect_true(all(is.finite(sm$coefficients$se[non_log])))
  expect_equal(sm$coefficients$mle, fit$coefficients$parameters$mle)
  expect_output(print(sm), "Parameters:")
})

test_that("coef and logLik fail on for_dev bundle", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  dev <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(verbose = FALSE, num_groups = 2L, num_threads = 1L),
    for_dev = TRUE
  )
  expect_error(stats::coef(dev), "fitted")
  expect_error(stats::logLik(dev), "fitted")
  expect_error(summary(dev), "fitted")
  expect_error(stats::vcov(dev), "fitted")
})
