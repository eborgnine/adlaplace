test_that("ad_fun_ptr dispatches to extension package via ad_data@package", {
  skip_if_not_installed("adlaplaceExample")

  obs <- adlaplace::ad_data(
    y = c(1, 2, 3),
    X = Matrix::Matrix(cbind(1, 1:3)),
    beta_map = Matrix::Diagonal(2L),
    gamma_map = Matrix::Matrix(nrow = 0L, ncol = 0L),
    theta_map = list(c(1L, 2L), 2L),
    ad_kind = "observations",
    ad_fun = "skewnormal_obs",
    package = "adlaplaceExample"
  )

  config <- list(
    beta = c(0, 0),
    theta = c(0, log(0.5)),
    transform_theta = TRUE,
    gamma = numeric(0),
    shards = adlaplace::ad_shards(
      Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1),
      num_shards = 1L
    ),
    verbose = FALSE
  )

  ptr <- adlaplace::ad_fun_ptr(obs, config)
  expect_true(inherits(ptr, "ad_fun_ptr"))
  expect_true(adlaplace::n_groups(ptr) >= 1L)
})
