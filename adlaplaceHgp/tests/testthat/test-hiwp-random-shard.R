test_that("hiwp builds a random_diagonal shard and matches design labels", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(
    y = rnorm(20),
    x = rep(seq(0, 1, length.out = 10), 2),
    g = rep(c("a", "b"), each = 10)
  )
  knots <- c(0, 0.5, 1)
  terms <- hiwp(
    x = "x", knots = knots, by = "g",
    include_global = FALSE, include_poly = FALSE
  )
  expect_equal(terms$x_hiwp@density, "random_diagonal")
  expect_equal(terms$x_hiwp@ad_kind, "random")

  a_cols <- colnames(adlaplace::design(terms$x_hiwp, data))
  q_cols <- colnames(adlaplace::precision(terms$x_hiwp, data))
  r_labels <- adlaplace::random_info(terms$x_hiwp, data)$gamma_label
  expect_identical(a_cols, q_cols)
  expect_identical(a_cols, r_labels)

  md <- adlaplace::model_data(
    c(list(y = adlaplace::gaussian("y")), terms),
    data = data,
    verbose = FALSE
  )
  expect_true("x_hiwp" %in% names(md$random))
  expect_equal(md$random$x_hiwp@density, "random_diagonal")
  expect_equal(ncol(md$random$x_hiwp@gamma_map), length(a_cols))
})
