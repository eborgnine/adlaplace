test_that("hnlm fit stores Woodbury var_iid when package available", {
  skip_if_not_installed("adlaplace")
  skip_if_not_installed("adlaplaceHgp")
  skip_if_not_installed("WoodburyMatrix")

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

  expect_s3_class(fit, "hnlm")
  expect_false(is.null(fit$hessian$inner))
  expect_false(is.null(fit$hessian$var_iid))
  expect_equal(nrow(fit$hessian$var_iid), nrow(fit$hessian$inner))
  expect_equal(ncol(fit$hessian$var_iid), ncol(fit$hessian$inner))
  expect_true(all(is.finite(fit$hessian$var_iid)))
})
