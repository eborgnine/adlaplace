test_that("replace_vars rewrites rsiid @mult", {
  term <- adlaplaceHgp::rsiid("state", mult = "sqrt_pm", ref_mult = 4)[[1]]
  out <- adlaplace::replace_vars(term, c(sqrt_pm = "sqrt_pm_s5"))
  expect_equal(out@name, "state")
  expect_equal(out@mult, "sqrt_pm_s5")
  expect_true(grepl("sqrt_pm_s5", out@label))
})

test_that("replace_vars rewrites hiwp name", {
  kn <- seq(0, 10, by = 2)
  term <- adlaplaceHgp::hiwp(
    "sqrt_pm",
    by = "state",
    ref_value = 4,
    p = 2,
    knots = kn,
    include_poly = FALSE,
    include_global = FALSE
  )[[1]]
  out <- adlaplace::replace_vars(term, c(sqrt_pm = "sqrt_pm_s2"))
  expect_equal(out@name, "sqrt_pm_s2")
  expect_equal(unname(out@by@term), "state")
})
