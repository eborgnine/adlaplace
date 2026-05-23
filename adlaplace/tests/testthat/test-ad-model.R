test_that("length-1 integer map args become zero-column matrices", {
  model <- adlaplace::ad_model(
    beta_map = 2L,
    gamma_map = 3L,
    theta_map = 1L
  )
  expect_equal(sizes(model), list(n_beta = 2L, n_gamma = 3L, n_theta = 1L))
  expect_equal(ncol(model@beta_map), 0L)
  expect_equal(ncol(model@gamma_map), 0L)
  expect_equal(ncol(model@theta_map), 0L)
})

test_that("length-2 map shorthand builds one-column maps", {
  model <- adlaplace::ad_model(
    beta_map = c(1L, 2L),
    theta_map = c(2L, 3L)
  )
  expect_equal(sizes(model), list(n_beta = 2L, n_gamma = 0L, n_theta = 3L))
  expect_equal(nrow(model@beta_map), 2L)
  expect_equal(ncol(model@beta_map), 1L)
  expect_equal(model@beta_map@i + 1L, 1L)
  expect_equal(nrow(model@theta_map), 3L)
  expect_equal(ncol(model@theta_map), 1L)
  expect_equal(model@theta_map@i + 1L, 2L)
})

test_that("length-2 map shorthand validates row index", {
  expect_error(
    adlaplace::ad_model(theta_map = c(0L, 1L)),
    "row index"
  )
  expect_error(
    adlaplace::ad_model(theta_map = c(4L, 3L)),
    "row index"
  )
})

test_that("random_diagonal ad_fun_ptr defaults precision$Q", {
  nr <- 3L
  model <- adlaplace::ad_model(
    gamma_map = Matrix::sparseMatrix(
      i = seq_len(nr),
      j = seq_len(nr),
      x = 1,
      dims = c(5L, nr)
    ),
    theta_map = c(1L, 1L)
  )
  config <- list(
    beta = numeric(0),
    gamma = rep(0, 5),
    theta = 0.1,
    transform_theta = FALSE
  )
  ptr <- adlaplace::ad_fun_ptr(
    config,
    kind = "random",
    name = "random_diagonal",
    model = model
  )
  expect_true(is(ptr, "ad_fun_ptr"))
})

test_that("gamma_map nrow defaults to ncol(A)", {
  A <- Matrix::sparseMatrix(
    i = c(0L, 1L),
    j = c(0L, 1L),
    x = 1,
    dims = c(2L, 3L),
    index1 = FALSE,
    giveCsparse = TRUE
  )
  X <- matrix(1, 2, 1)
  model <- adlaplace::ad_model(y = c(1, 2), A = A, X = X, theta_map = Matrix::Diagonal(1L))
  expect_equal(sizes(model)$n_gamma, 3L)
})

test_that("infer n_theta from info$parameters and info$gamma", {
  info <- list(
    parameters = data.frame(a = 1:2),
    gamma = data.frame(b = 1:4),
    theta = data.frame(c = 1)
  )
  model <- adlaplace::ad_model(
    y = 1:4,
    A = matrix(0, 4, 4),
    X = matrix(0, 4, 2),
    info = info
  )
  expect_equal(sizes(model)$n_beta, 2L)
  expect_equal(sizes(model)$n_gamma, 4L)
  expect_equal(sizes(model)$n_theta, 1L)
})

test_that("ad_model stops when n_theta cannot be inferred", {
  A <- Matrix::sparseMatrix(
    i = 0L,
    j = 0L,
    x = 1,
    dims = c(1L, 1L),
    index1 = FALSE,
    giveCsparse = TRUE
  )
  expect_error(
    adlaplace::ad_model(y = 1, A = A, X = matrix(1, 1, 1)),
    "cannot infer n_theta"
  )
})

test_that("validate_config_layout infers missing config$theta length", {
  A <- Matrix::sparseMatrix(
    i = c(0L, 1L),
    j = c(0L, 1L),
    x = 1,
    dims = c(2L, 2L),
    index1 = FALSE,
    giveCsparse = TRUE
  )
  info <- list(
    parameters = data.frame(x = 1:3),
    gamma = data.frame(y = 1:2),
    theta = data.frame(z = 1)
  )
  model <- adlaplace::ad_model(y = 1:2, A = A, X = matrix(1, 2, 2), info = info)
  config <- list(
    beta = rep(0, 2),
    gamma = rep(0, 2)
  )
  expect_silent(adlaplace:::validate_config_layout(model, config))
})
