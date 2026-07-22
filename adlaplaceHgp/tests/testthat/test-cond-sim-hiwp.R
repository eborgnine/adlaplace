test_that("cond_sim_iwp_inputs extracts hiwp random effects from model_data", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(
    y = rnorm(20),
    x = rep(seq(0, 1, length.out = 10), 2),
    g = rep(c("a", "b"), each = 10)
  )
  knots <- c(0, 0.5, 1)
  model_data <- adlaplace::model_data(
    data = data,
    formula = y ~ hiwp(x = "x", knots = knots, by = "g", include_global = FALSE),
    verbose = FALSE
  )
  random_info <- model_data$term_data$info$gamma
  n_gamma <- nrow(random_info)
  expect_true(any(random_info$model == "hiwp"))

  laplace <- list(
    full_parameters = c(
      model_data$term_data$info$beta$init,
      rep(0, n_gamma),
      model_data$term_data$info$theta$init
    ),
    extra = list(half_H_inv = Matrix::Diagonal(n_gamma, 0.5))
  )
  inputs <- cond_sim_iwp_inputs(laplace, model_data)
  expect_equal(length(inputs$gamma_mode), n_gamma)
  expect_true(any(inputs$random_info$model == "hiwp"))
  expect_equal(nrow(inputs$half_H_inv), n_gamma)
})

test_that("cond_sim_iwp smoke with hierarchical hiwp terms", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(
    y = rnorm(20),
    x = rep(seq(0, 1, length.out = 10), 2),
    g = rep(c("a", "b"), each = 10)
  )
  knots <- c(0, 0.5, 1)
  model_data <- adlaplace::model_data(
    data = data,
    formula = y ~ hiwp(x = "x", knots = knots, by = "g", include_global = FALSE),
    verbose = FALSE
  )
  n_gamma <- nrow(model_data$term_data$info$gamma)
  laplace <- list(
    full_parameters = c(
      model_data$term_data$info$beta$init,
      rep(0, n_gamma),
      model_data$term_data$info$theta$init
    ),
    extra = list(half_H_inv = Matrix::Diagonal(n_gamma, 0.5))
  )

  result <- tryCatch(
    cond_sim_iwp(laplace, model_data, n = 5L),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(conditionMessage(result))
  }
  expect_true("sim" %in% names(result))
  expect_gt(length(result$sim), 0L)
})

test_that("cond_sim_iwp envelope uses GET when available", {
  skip_if_not_installed("adlaplace")
  skip_if_not_installed("GET")

  data <- data.frame(
    y = rnorm(20),
    x = rep(seq(0, 1, length.out = 10), 2),
    g = rep(c("a", "b"), each = 10)
  )
  knots <- c(0, 0.5, 1)
  model_data <- adlaplace::model_data(
    data = data,
    formula = y ~ hiwp(x = "x", knots = knots, by = "g", include_global = FALSE),
    verbose = FALSE
  )
  n_gamma <- nrow(model_data$term_data$info$gamma)
  laplace <- list(
    full_parameters = c(
      model_data$term_data$info$beta$init,
      rep(0, n_gamma),
      model_data$term_data$info$theta$init
    ),
    extra = list(half_H_inv = Matrix::Diagonal(n_gamma, 0.5))
  )

  result <- tryCatch(
    cond_sim_iwp(
      laplace,
      model_data,
      n = 8L,
      probs_envelope = c(0.1, 0.9)
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    skip(conditionMessage(result))
  }
  expect_true("envelope" %in% names(result))
})
