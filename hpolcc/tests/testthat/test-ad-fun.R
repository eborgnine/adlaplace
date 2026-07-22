test_that("model_data attaches elgm_matrix for dirichlet_multinom", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  model_terms <- adlaplace::collect_terms(
    stats::update.formula(td$formula, . ~ . - 1)
  )
  md <- adlaplace::model_data(
    formula = model_terms,
    data = td$data,
    na_omit = TRUE,
    verbose = FALSE
  )
  expect_gt(ncol(md$term_data$elgm_matrix), 0L)
  expect_gt(ncol(md$observations$count@elgm_matrix), 0L)
  expect_identical(md$parameters$count_extra@density, "dirichlet_multinomial_extra")
})

test_that("ad_pack builds from hnlm for_dev model_data bundle", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_groups = 2L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )
  expect_true(methods::is(forres$ad_pack, "ad_pack"))
  expect_true("outer" %in% methods::slotNames(forres$ad_pack))

  sz <- forres$ad_pack@sizes
  n_beta <- as.integer(sz["beta"])
  n_gamma <- as.integer(sz["gamma"])
  n_theta <- as.integer(sz["theta"])
  x <- c(
    forres$config$opt$init[seq_len(n_beta)],
    rep(0, n_gamma),
    forres$config$opt$init[seq(n_beta + 1L, length.out = n_theta)]
  )

  dens <- adlaplace::joint_log_dens(forres$ad_pack, x)
  expect_true(is.finite(dens))

  g <- adlaplace::grad(forres$ad_pack, x)
  expect_equal(length(g), length(x))

  H <- adlaplace::hessian(forres$ad_pack, x)
  expect_true(inherits(H, "Matrix") || inherits(H, "matrix"))
})

test_that("ad_pack builds from hnlm for_dev with no fixed effects", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  formula <- adlaplace::dirichlet_multinom(
    count,
    by = c("year", "region", "date")
  ) ~ adlaplace::iid(date)
  forres <- hpolcc::hnlm(
    formula = formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_groups = 2L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )

  expect_equal(as.integer(forres$ad_pack@sizes["beta"]), 0L)
  expect_equal(nrow(forres$model_data$term_data$info$beta), 0L)

  n_gamma <- as.integer(forres$ad_pack@sizes["gamma"])
  n_theta <- as.integer(forres$ad_pack@sizes["theta"])
  x <- c(
    rep(0, n_gamma),
    forres$config$opt$init[seq_len(n_theta)]
  )

  dens <- adlaplace::joint_log_dens(forres$ad_pack, x)
  expect_true(is.finite(dens))

  g <- adlaplace::grad(forres$ad_pack, x)
  expect_equal(length(g), length(x))
})
