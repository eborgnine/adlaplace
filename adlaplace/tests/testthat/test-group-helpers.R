# model_term has no by slot; use a local subclass for add_by_levels.
methods::setClass(
  "group_helper_term",
  slots = c(by = "by_group"),
  contains = "model_term"
)

test_that("add_by_levels fills empty levels from data and rejects extras", {
  term <- methods::new(
    "group_helper_term",
    name = "x",
    by = by_group(term = "g")
  )
  filled <- add_by_levels(term, data.frame(g = c("B", "A", "B")))
  expect_equal(filled@by@levels, c("B", "A"))

  fact <- add_by_levels(
    term,
    data.frame(g = factor(c("A", "A"), levels = c("B", "A")))
  )
  expect_equal(fact@by@levels, c("B", "A"))

  locked <- methods::new(
    "group_helper_term",
    name = "x",
    by = by_group(term = "g", levels = c("A", "B"))
  )
  ok <- add_by_levels(locked, data.frame(g = c("B", "A")))
  expect_equal(ok@by@levels, c("A", "B"))
  expect_error(
    add_by_levels(locked, data.frame(g = c("A", "C"))),
    "levels not in the term"
  )
  expect_error(
    add_by_levels(locked, data.frame(g = c("A", NA))),
    "missing values"
  )
  expect_identical(add_by_levels(methods::new("model_term", name = "x"), data.frame(g = "A"))@name, "x")
})

test_that("group_index matches levels and errors on NA or unknown levels", {
  data <- data.frame(g = c("B", "A", "B"))
  expect_equal(group_index(data, "g", c("A", "B")), c(2L, 1L, 2L))
  expect_error(group_index(data, "missing", "A"), "not found")
  expect_error(group_index(data.frame(g = c("A", NA)), "g", "A"), "missing values")
  expect_error(group_index(data, "g", "A"), "levels not in the term")
})

test_that("scatter_by_group writes one block per group and labels match", {
  B <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L),
    j = c(1L, 1L, 2L),
    x = c(1, 2, 3),
    dims = c(2L, 2L)
  )
  group <- c(2L, 1L)
  A <- scatter_by_group(B, group, 2L)
  expect_equal(dim(A), c(2L, 4L))
  # Group 1 (row 2) occupies columns 1:2; group 2 (row 1) occupies columns 3:4.
  expect_equal(as.numeric(A[2L, 1L]), 2)
  expect_equal(as.numeric(A[2L, 2L]), 3)
  expect_equal(as.numeric(A[1L, 3L]), 1)
  expect_equal(as.numeric(A[1L, c(1L, 2L, 4L)]), c(0, 0, 0))

  dense <- scatter_by_group(as.matrix(B), group, 2L)
  expect_equal(as.matrix(dense), as.matrix(A))

  labels <- group_gamma_labels(c("k1", "k2"), c("A", "B"))
  expect_equal(labels, c("k1_gA", "k2_gA", "k1_gB", "k2_gB"))
  expect_equal(length(labels), ncol(A))
})
