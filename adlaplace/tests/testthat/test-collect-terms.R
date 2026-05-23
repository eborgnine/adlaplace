test_that("collect_terms parses linear(x) without x in caller env", {
  f <- y ~ intercept() + linear(x)
  terms <- adlaplace::collect_terms(f)
  expect_true(inherits(terms$`linear(x)`, "linear"))
  expect_equal(terms$`linear(x)`@term, "x")
})

test_that("collect_terms parses bare covariate as linear", {
  f <- y ~ x1
  terms <- adlaplace::collect_terms(f)
  expect_true(inherits(terms$x1, "linear"))
  expect_equal(terms$x1@term, "x1")
})
