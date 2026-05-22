test_that("hnlm for_dev builds ad_fun via get_ad_fun", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    cc_design = td$cc_design,
    config = list(
      verbose = FALSE,
      num_groups = 2L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )
  expect_true(is.data.frame(forres$data))
  expect_true(is.list(forres$ad_fun))
  expect_true("ad_fun" %in% names(forres$ad_fun))

  ad_fun2 <- adlaplace::get_ad_fun(
    forres$model$data,
    forres$config,
    package = "hpolcc"
  )
  expect_true(is.list(ad_fun2))
  x <- c(forres$config$beta, forres$config$theta)
  ll <- adlaplace::log_lik_laplace(
    x = x,
    config = forres$config,
    ad_fun = ad_fun2,
    deriv = FALSE,
    control = list(maxit = 20, report.level = 0)
  )
  expect_true(is.finite(ll$log_lik))
  expect_equal(ll$neg_log_lik, -ll$log_lik)
})
