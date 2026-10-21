test_that("aic recovers a short autoregression and returns its polynomial", {
  set.seed(1)
  y <- as.numeric(arima.sim(n = 400, list(ar = c(0.6, -0.25)), sd = 1))
  poly <- aic(y, 6)

  expect_equal(poly[1], 1)
  expect_equal(length(poly) - 1, 2)
  expect_equal(poly[-1], c(-0.6, 0.25), tolerance = 0.15)
  expect_equal(length(attr(poly, "eef")), length(y))
  expect_s3_class(attr(poly, "aic"), "data.frame")
  expect_error(aic(y, 4, out = 2), "out must be 0 or 1")
})

test_that("subset least squares expands excluded lags to zero", {
  set.seed(2)
  y <- as.numeric(arima.sim(n = 300, list(ar = c(0.4, rep(0, 10), 0.3)), sd = 1))
  fit <- captaintvp:::ar_subset(y, c(1, 12))

  expect_equal(length(fit$ARp), 12)
  expect_equal(fit$ARp[c(2:11)], rep(0, 10))
  expect_equal(fit$ARp[c(1, 12)], c(-0.4, -0.3), tolerance = 0.15)
  expect_true(all(is.na(fit$ARpse[2:11])))
  expect_true(all(is.finite(fit$ARpse[c(1, 12)])))
})

test_that("missing observations are not used as autoregressive responses", {
  y <- sin(2 * pi * (1:40) / 12)
  y[10:15] <- NA
  fit <- captaintvp:::ar_subset(y, 1:2)
  expect_true(fit$n < sum(is.finite(y)))
  expect_true(fit$sigma2 > 0)
})
