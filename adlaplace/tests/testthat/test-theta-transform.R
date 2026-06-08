test_that("data_setup defaults transform to TRUE for all theta rows", {
  skip_if_not_installed("mgcv")
  dat <- mgcv::gamSim(6, n = 80, scale = 0.2, dist = "poisson")
  md <- adlaplace::model_data(
    adlaplace::f(y, model = "nbinom", lower = 1e-9) ~
      x1 +
      adlaplace::iid(fac, init = 0.25),
    data = dat,
    verbose = FALSE
  )
  theta <- md$data$info$theta
  expect_true("transform" %in% names(theta))
  expect_false(any(is.na(theta$transform)))
  expect_true(all(theta$transform))
})

test_that("apply_theta_log logs only rows with transform TRUE", {
  theta_info <- data.frame(
    label = c("sd", "alpha"),
    init = c(0.5, 0.8),
    lower = c(1e-9, -Inf),
    upper = c(Inf, Inf),
    transform = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  out <- adlaplace::apply_theta_log(
    theta_info,
    cols = c("init", "lower", "upper")
  )
  expect_equal(out$init, c(log(0.5), 0.8))
  expect_equal(out$lower[1], log(1e-9))
  expect_equal(out$lower[2], -Inf)
  expect_equal(out$upper, c(Inf, Inf))
})

test_that("apply_theta_log is a no-op when active is FALSE", {
  theta_info <- data.frame(
    label = "sd",
    init = 0.5,
    transform = TRUE,
    stringsAsFactors = FALSE
  )
  out <- adlaplace::apply_theta_log(theta_info, active = FALSE)
  expect_equal(out$init, 0.5)
})
