test_that("getAdFun_r returns full backend list and derivatives run", {
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
  ad_fun <- hpolcc::getAdFun_r(forres$model$data, forres$config)
  expect_true(is.list(ad_fun))
  expect_true("ad_fun" %in% names(ad_fun))
  expect_true("outer" %in% names(ad_fun))
  expect_true("map_outer" %in% names(ad_fun))

  x <- c(forres$config$beta, forres$config$gamma, forres$config$theta)
  dens <- adlaplace::joint_log_dens(ad_fun, x)
  expect_true(is.finite(dens))

  g <- adlaplace::grad(ad_fun, x)
  expect_equal(length(g), length(x))

  H <- adlaplace::hessian(ad_fun, x)
  expect_true(inherits(H, "Matrix") || inherits(H, "matrix"))
})
