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
