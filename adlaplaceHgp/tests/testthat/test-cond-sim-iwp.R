test_that("cond_sim_iwp_at builds summaries from flat arguments", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(x = seq(0, 1, length.out = 10))
  knots <- c(0, 0.5, 1)
  iwp_terms <- adlaplace::iwp(x = "x", knots = knots, include_poly = FALSE)
  iwp_term <- iwp_terms[[1]]
  random_info <- adlaplace::random_info(iwp_term, data)
  n_gamma <- nrow(random_info)

  result <- cond_sim_iwp_at(
    terms = unlist(iwp_terms, recursive = FALSE),
    random_info = random_info,
    beta = numeric(0),
    gamma_mode = stats::setNames(rep(0, n_gamma), random_info$gamma_label),
    half_H_inv = Matrix::Diagonal(n_gamma, 0.5),
    n = 5L
  )
  expect_true("sim" %in% names(result))
  expect_true("quantiles" %in% names(result))
  expect_gt(length(result$sim), 0L)
})

test_that("cond_sim_iwp_inputs and cond_sim_iwp use laplace + model_data", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(x = seq(0, 1, length.out = 10))
  knots <- c(0, 0.5, 1)
  iwp_terms <- adlaplace::iwp(x = "x", knots = knots, include_poly = FALSE)
  iwp_term <- iwp_terms[[1]]
  random_info <- adlaplace::random_info(iwp_term, data)
  n_gamma <- nrow(random_info)

  model_data <- list(
    terms = unlist(iwp_terms, recursive = FALSE),
    term_data = list(info = list(
      beta = data.frame(beta_label = character(0)),
      gamma = random_info,
      theta = data.frame(label = "x_iwp")
    ))
  )
  laplace <- list(
    full_parameters = c(rep(0, n_gamma), -1),
    extra = list(
      half_H_inv = Matrix::Diagonal(n_gamma, 0.5)
    )
  )

  inputs <- cond_sim_iwp_inputs(laplace, model_data)
  expect_equal(
    names(inputs),
    c("terms", "random_info", "beta", "gamma_mode", "half_H_inv")
  )
  expect_equal(length(inputs$beta), 0L)
  expect_equal(length(inputs$gamma_mode), n_gamma)

  result <- cond_sim_iwp(laplace, model_data, n = 5L)
  expect_true("sim" %in% names(result))
  expect_gt(length(result$sim), 0L)
})
