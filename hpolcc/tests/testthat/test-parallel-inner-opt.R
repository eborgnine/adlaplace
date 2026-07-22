# Dirichlet-multinom obs/extra live in adlaplace.so; random (iid) also in adlaplace.so.
# Parallel OpenMP + CppAD teardown is flaky under full R CMD check on macOS
# (thread_alloc returns from the wrong thread). Coverage for the parallel
# Laplace path lives in adlaplace's test-cppad-teardown.R /
# test-log-lik-deriv-parallel.R; this file keeps a serial for_dev smoke check.

test_that("for_dev bundle log_lik_laplace deriv=TRUE (serial)", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_groups = 8L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )
  ll <- adlaplace::log_lik_laplace(
    x = forres$config$opt$init,
    config = list(verbose = FALSE),
    gamma = forres$cache$gamma,
    ad_pack = forres$ad_pack,
    control = list(maxit = 5L, report.level = 0, report.freq = 0),
    deriv = TRUE
  )
  expect_true(is.finite(ll$log_lik))
  expect_true(all(is.finite(ll$deriv$d_neg_log_lik)))
  expect_true(all(is.finite(ll$extra$trace3)))
})
