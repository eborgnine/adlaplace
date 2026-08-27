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
  expect_equal(n_mapped, nrow(md$term_data$info$gamma))
  expect_equal(ncol(md$term_data$A), nrow(md$term_data$info$gamma))
  expect_false(anyDuplicated(colnames(md$term_data$A)) > 0)
})

test_that("rsiwp + iwp model_data A is sparse (IWP-1 shape)", {
  skip_if_not_installed("adlaplace")
  set.seed(5)
  n <- 60L
  dat <- data.frame(
    y = rpois(n, 2),
    time_years = seq(2000, 2016, length.out = n) / 365.25 * 365.25,
    sqrt_pm = sqrt(runif(n, 1, 25)),
    temperature_div10 = runif(n, -2, 3)
  )
  time_knots <- seq(
    min(dat$time_years), max(dat$time_years), length.out = 5
  )
  temp_knots <- seq(-2.5, 3.5, by = 1)
  terms <- unlist(list(
    adlaplace::gaussian("y"),
    rsiwp(
      "time_years",
      mult = "sqrt_pm",
      p = 1,
      ref_value = mean(time_knots),
      ref_mult = 4,
      knots = time_knots,
      include_linear = TRUE
    ),
    adlaplace::iwp(
      "temperature_div10",
      p = 3,
      ref_value = 2,
      knots = temp_knots,
      init = 0.01
    )
  ), recursive = FALSE)
  md <- adlaplace::model_data(terms, data = dat, verbose = FALSE)
  expect_true(inherits(md$term_data$A, "sparseMatrix"))
  expect_s4_class(md$term_data$A, "dgCMatrix")
  expect_equal(ncol(md$term_data$A), nrow(md$term_data$info$gamma))
  # Dense rpoly from iwp(p = 3) is cbind-ed but A stays sparse.
  expect_true(any(grepl("rpoly", colnames(md$term_data$A))))
  expect_true(any(grepl("rsiwp", colnames(md$term_data$A))))
})
