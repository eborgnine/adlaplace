# Mixed-DSO dirichlet_multinom: obs/extra in hpolcc.so, random (iid) in adlaplace.so.

skip_parallel_on_macos <- function() {
  skip_if(
    Sys.info()[["sysname"]] == "Darwin" &&
      Sys.getenv("ADLAPLACE_TEST_PARALLEL_INNER_OPT", "") != "1",
    "parallel inner_opt test (set ADLAPLACE_TEST_PARALLEL_INNER_OPT=1 on macOS)"
  )
}

test_that("parallel log_lik_laplace deriv=FALSE via hnlm for_dev bundle", {
  skip_if_not_installed("adlaplace")
  skip_parallel_on_macos()
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_shards = 8L,
      num_threads = 4L
    ),
    for_dev = TRUE
  )
  sz <- forres$ad_fun@sizes
  n_beta <- as.integer(sz["beta"])
  n_theta <- as.integer(sz["theta"])
  ll <- adlaplace::log_lik_laplace(
    x = forres$config$opt$init,
    config = list(verbose = FALSE),
    gamma = forres$cache$gamma,
    ad_fun = forres$ad_fun,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = FALSE
  )
  expect_true(is.finite(ll$log_lik))
  expect_equal(length(forres$config$opt$init), n_beta + n_theta)
})

test_that("parallel log_lik_laplace deriv=TRUE via hnlm for_dev bundle", {
  skip_if_not_installed("adlaplace")
  skip_parallel_on_macos()
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_shards = 8L,
      num_threads = 4L
    ),
    for_dev = TRUE
  )
  ll <- adlaplace::log_lik_laplace(
    x = forres$config$opt$init,
    config = list(verbose = FALSE),
    gamma = forres$cache$gamma,
    ad_fun = forres$ad_fun,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  expect_true(is.finite(ll$log_lik))
  expect_true(all(is.finite(ll$grad)))
  expect_true(all(is.finite(ll$extra$trace3)))
})
