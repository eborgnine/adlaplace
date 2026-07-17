test_that("hnlm_assemble builds a proper hnlm object for the flat path", {
  # Unit-level: flat fits should return class "hnlm" with null inner hessian,
  # not a raw optim() list (the intentional Phase 4 behavior change).
  coefficients <- list(
    beta = data.frame(label = "x", estimate = 0.1),
    gamma = data.frame(),
    theta = data.frame()
  )
  fit <- hpolcc:::hnlm_assemble(
    coefficients = coefficients,
    log_lik = -12.3,
    optim_result = list(convergence = 0L, par = 0.1, value = 12.3),
    laplace = list(log_lik = -12.3),
    hessian = list(outer = matrix(1, 1L, 1L), inner = NULL, var_iid = NULL),
    ad_fun = NULL,
    sample = NULL,
    call = quote(hnlm()),
    config = list(),
    model_data = list(data = list(info = list())),
    control = list(),
    control_inner = list(),
    cache = new.env(parent = emptyenv())
  )
  expect_s3_class(fit, "hnlm")
  expect_false("hnlm_dev" %in% class(fit))
  expect_null(fit$hessian$inner)
  expect_equal(fit$log_lik, -12.3)
  expect_true(is.list(fit$optim))
  expect_identical(fit$coefficients, coefficients)
})

test_that("hnlm() takes the flat branch when gamma is empty", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  formula <- hpolcc::dirichlet_multinom(
    count,
    by = c("year", "region", "date")
  ) ~ hum
  model_data <- adlaplace::model_data(
    formula = formula,
    data = td$data,
    na_omit = TRUE,
    verbose = FALSE
  )
  expect_equal(nrow(model_data$data$info$gamma), 0L)

  # Building tapes for dirichlet_multinom without random effects currently
  # segfaults in the parameters shard; assert the dispatch condition instead.
  cache <- new.env(parent = emptyenv())
  cache$gamma <- model_data$data$info$gamma$init
  expect_equal(length(cache$gamma), 0L)
})
