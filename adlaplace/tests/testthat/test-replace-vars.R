test_that("replace_vars rewrites model_term name/label/formula", {
  term <- adlaplace::fpoly("sqrt_pm", p = 1, ref_value = 4)
  out <- adlaplace::replace_vars(term, c(sqrt_pm = "sqrt_pm_s3"))
  expect_equal(out@name, "sqrt_pm_s3")
  expect_true(grepl("sqrt_pm_s3", out@label))
  expect_true(grepl("sqrt_pm_s3", deparse(out@formula)))
})

test_that("replace_vars on list of terms preserves structure", {
  terms <- list(
    adlaplace::fpoly("sqrt_pm", p = 1, ref_value = 4),
    adlaplace::dirichlet_multinom("mort_all", by = c("geo_id", "yearMonthDow"))
  )
  out <- adlaplace::replace_vars(
    terms,
    c(sqrt_pm = "sqrt_pm_s1", geo_id = "region")
  )
  expect_equal(out[[1]]@name, "sqrt_pm_s1")
  expect_equal(out[[2]]@by, c("region", "yearMonthDow"))
})

test_that("replace_vars on formula rewrites symbols", {
  f <- y ~ sqrt_pm + temperature
  out <- adlaplace::replace_vars(f, c(sqrt_pm = "sqrt_pm_s2"))
  expect_equal(all.vars(out), c("y", "sqrt_pm_s2", "temperature"))
})
