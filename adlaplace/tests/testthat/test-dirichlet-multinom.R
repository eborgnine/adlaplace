test_that("model_data attaches elgm_matrix for dirichlet_multinom", {
  set.seed(0)
  n_strata <- 6L
  n_per <- 4L
  n <- n_strata * n_per
  data <- data.frame(
    count = rpois(n, lambda = 2),
    hum = rnorm(n),
    region = rep(1:2, length.out = n),
    date = rep(seq_len(n_strata), each = n_per),
    year = 2002L
  )
  formula <- dirichlet_multinom(
    count,
    by = c("year", "region", "date")
  ) ~ hum + iid(date)

  md <- model_data(
    formula = formula,
    data = data,
    na_omit = TRUE,
    verbose = FALSE
  )
  expect_gt(ncol(md$term_data$elgm_matrix), 0L)
  expect_gt(ncol(md$observations$count@elgm_matrix), 0L)
  expect_identical(md$parameters$count_extra@density, "dirichlet_multinomial_extra")
  expect_identical(md$observations$count@package, "adlaplace")
})
