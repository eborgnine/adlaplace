test_that("model_data skewnormal theta has selective transform column", {
  skip_if_not_installed("sn")
  set.seed(0L)
  Nobs <- 50L
  dat <- data.frame(
    y = sn::rsn(Nobs, xi = 0, omega = 0.7, alpha = 0.8),
    x = stats::rnorm(Nobs),
    r1 = sample(3L, Nobs, replace = TRUE)
  )
  md <- adlaplace::model_data(
    adlaplaceExample::skewnormal(y, init = c(0.7, 0.8)) ~
      x + adlaplace::iid(r1, init = 0.1),
    data = dat,
    verbose = FALSE
  )
  theta <- md$data$info$theta
  expect_true("log" %in% names(theta))
  expect_false(any(is.na(theta$log)))

  alpha_row <- grepl("_alpha$", theta$label)
  omega_row <- grepl("_omega$", theta$label)
  sd_row <- theta$model == "iid"
  expect_true(all(theta$log[sd_row]))
  expect_true(all(theta$log[omega_row]))
  expect_false(any(theta$log[alpha_row]))
})

test_that("ad_fun(model_data) leaves alpha untransformed in config theta", {
  skip_if_not_installed("sn")
  set.seed(1L)
  Nobs <- 40L
  dat <- data.frame(
    y = sn::rsn(Nobs, xi = 0, omega = 0.7, alpha = 0.8),
    x = stats::rnorm(Nobs),
    r1 = sample(3L, Nobs, replace = TRUE)
  )
  md <- adlaplace::model_data(
    adlaplaceExample::skewnormal(y, init = c(0.7, 0.8)) ~
      x + adlaplace::iid(r1, init = 0.1),
    data = dat,
    verbose = FALSE
  )
  expected_theta <- adlaplace::apply_theta_log(md$data$info$theta, cols = "init")$init

  config <- list(
    transform_theta = TRUE,
    shards = adlaplace::ad_shards(md$data$A, num_shards = 5L),
    verbose = FALSE
  )
  ad_fun <- adlaplace::ad_fun(md, config, num_threads = 1L)

  alpha_label <- md$data$info$theta$label[grepl("_alpha$", md$data$info$theta$label)]
  alpha_idx <- match(alpha_label, md$data$info$theta$label)
  expect_equal(expected_theta[alpha_idx], md$data$info$theta$init[alpha_idx])

  x_full <- c(
    md$data$info$beta$init,
    rep(0, nrow(md$data$info$gamma)),
    expected_theta
  )
  expect_true(is.finite(adlaplace::joint_log_dens(ad_fun, x_full, negative = FALSE)))
})
