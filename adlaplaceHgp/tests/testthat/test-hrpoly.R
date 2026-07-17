test_that("hrpoly design and precision match group structure", {
  skip_if_not_installed("adlaplace")

  data <- data.frame(
    age = rep(seq(0, 40, by = 10), each = 3),
    site = rep(c("A", "B", "C"), times = 5)
  )
  term <- hrpoly(x = "age", p = 2, ref_value = 0, by = "site")
  expect_s4_class(term, "hrpoly")

  A <- design(term, data)
  P <- precision(term, data)
  expect_s4_class(A, "dgCMatrix")
  expect_equal(nrow(A), nrow(data))
  expect_equal(ncol(A), length(unique(data$site)))
  expect_equal(nrow(P), ncol(A))
  expect_equal(ncol(P), ncol(A))
  expect_true(all(Matrix::diag(P) > 0))

  gamma_info <- random_info(term, data)
  expect_equal(nrow(gamma_info), ncol(A))
  expect_equal(gamma_info$gamma_label, colnames(A))
  expect_null(beta_info(term, data))
})
