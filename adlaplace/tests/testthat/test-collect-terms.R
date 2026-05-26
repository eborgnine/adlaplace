test_that("collect_terms parses linear(x) without x in caller env", {
  f <- y ~ intercept() + linear(x)
  terms <- adlaplace::collect_terms(f)
  expect_true(inherits(terms$`linear(x)`, "linear"))
  expect_equal(terms$`linear(x)`@term, "x")
  expect_equal(terms$`linear(x)`@label, "x_linear")
})

test_that("collect_terms parses bare covariate as linear", {
  f <- y ~ x1
  terms <- adlaplace::collect_terms(f)
  expect_true(inherits(terms$x1, "linear"))
  expect_equal(terms$x1@term, "x1")
})

test_that("collect_terms keeps bare outcome LHS out of model terms", {
  terms <- adlaplace::collect_terms(y ~ x1)
  expect_false("response" %in% names(terms))
  expect_false(any(vapply(terms, function(x) {
    methods::is(x, "nbinom")
  }, logical(1L))))
})

test_that("random_info reuses precomputed term label", {
  term <- adlaplace::iid("grp")[[1L]]
  expect_equal(term@label, "grp_iid")
  info <- adlaplace::random_info(term, data.frame(grp = factor(c("a", "b", "a"))))
  expect_true(all(info$label == term@label))
})
