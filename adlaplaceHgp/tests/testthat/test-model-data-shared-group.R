test_that("model_data keeps both random terms that share a grouping variable", {
  skip_if_not_installed("adlaplace")
  set.seed(1)
  dat <- data.frame(
    y = rpois(40, 2),
    grp = factor(rep(letters[1:4], each = 10)),
    x1 = runif(40),
    x2 = runif(40)
  )
  md <- adlaplace::model_data(
    list(
      adlaplace::gaussian("y"),
      rsiid("grp", mult = "x1", ref_mult = 0)[[1]],
      rsiid("grp", mult = "x2", ref_mult = 0)[[1]]
    ),
    data = dat,
    verbose = FALSE
  )
  expect_equal(
    sort(names(md$random)),
    sort(c("grp_x1_rsiid", "grp_x2_rsiid"))
  )
  n_mapped <- sum(vapply(
    md$random,
    function(r) Matrix::nnzero(r@gamma_map),
    integer(1)
  ))
  expect_equal(n_mapped, nrow(md$data$info$gamma))
  expect_equal(ncol(md$data$A), nrow(md$data$info$gamma))
  expect_false(anyDuplicated(colnames(md$data$A)) > 0)
})
