published_arp <- c(
  -0.6076, -0.1137, -0.0526, 0.1670, -0.1447,
  0.2348, -0.1281, 0.1891, -0.1225, 0.2281,
  -0.1813, -0.8041, 0.4586, 0.3524
)

air_demo_series <- function() {
  fcast(as.numeric(AirPassengers), rbind(c(84, 94), c(131, 144)))
}

test_that("aic puts seasonal orders well ahead of short perturbations", {
  y <- air_demo_series()
  initial <- univ(y, numeric(), 1, 1e-4)
  residual <- y - initial$trend
  poly <- aic(residual, 20)
  table <- attr(poly, "aic")

  expect_gte(length(poly) - 1, 12)
  order14 <- table$aic[table$order == 14]
  expect_lt(order14, table$aic[table$order == 1])
  expect_lt(order14, table$aic[table$order == 8])
})

test_that("the published air-demo AR(14) model is recovered up to smoother tolerance", {
  y <- air_demo_series()
  expect_output(
    est <- univopt(y, 1:14, 1, 1e-4, out = 1),
    "Final trend NVR estimate"
  )

  large <- abs(published_arp) > 0.1
  expect_equal(sign(est$ARp[large]), sign(published_arp[large]))
  expect_lt(max(abs(est$ARp - published_arp)), 0.12)
  expect_lt(sqrt(mean((est$ARp - published_arp)^2)), 0.08)
  expect_lt(abs(est$nvr - 2.5934e-3) / 2.5934e-3, 0.40)
  expect_equal(est$ARt, numeric())
  expect_true(all(is.finite(est$ARpse)))

  sm <- univ(y, est$ARp, 1, est$nvr)
  expect_equal(as.numeric(sm$fit), as.numeric(sm$trend) + as.numeric(sm$comp))
  expect_true(all(is.finite(sm$fit)))
  expect_true(all(is.finite(sm$y0)))
  expect_gt(as.numeric(sm$fitse)[144], as.numeric(sm$fitse)[131])
  held_out <- !is.finite(y)
  expect_equal(as.numeric(sm$y0)[held_out], as.numeric(sm$fit)[held_out])
  expect_equal(as.numeric(sm$y0)[!held_out], y[!held_out])
})
