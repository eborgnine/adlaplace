test_that("model_data attaches elgm_matrix for dirichlet_multinom", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  model_terms <- adlaplace::collect_terms(
    stats::update.formula(td$formula, . ~ . - 1),
    package = c("hpolcc", "adlaplaceHgp")
  )
  md <- adlaplace::model_data(
    formula = model_terms,
    data = td$data,
    na_omit = TRUE,
    verbose = FALSE
  )
  expect_gt(ncol(md$data$elgm_matrix), 0L)
  expect_gt(ncol(md$observations$count@elgm_matrix), 0L)
  expect_identical(md$parameters$count_extra@ad_fun, "dirichlet_multinomial_extra")
})

test_that("ad_fun builds from hnlm for_dev model_data bundle", {
  skip_if_not_installed("adlaplace")
  td <- make_hpolcc_test_data()
  forres <- hpolcc::hnlm(
    formula = td$formula,
    data = td$data,
    config = list(
      verbose = FALSE,
      num_shards = 2L,
      num_threads = 1L
    ),
    for_dev = TRUE
  )
  expect_true(methods::is(forres$ad_fun, "ad_fun"))
  expect_true("outer" %in% methods::slotNames(forres$ad_fun))

  sz <- forres$ad_fun@sizes
  n_beta <- as.integer(sz["beta"])
  n_gamma <- as.integer(sz["gamma"])
  n_theta <- as.integer(sz["theta"])
  x <- c(
    forres$config$opt$init[seq_len(n_beta)],
    rep(0, n_gamma),
    forres$config$opt$init[seq(n_beta + 1L, length.out = n_theta)]
  )

  dens <- adlaplace::joint_log_dens(forres$ad_fun, x)
  expect_true(is.finite(dens))

  g <- adlaplace::grad(forres$ad_fun, x)
  expect_equal(length(g), length(x))

  H <- adlaplace::hessian(forres$ad_fun, x)
  expect_true(inherits(H, "Matrix") || inherits(H, "matrix"))
})
