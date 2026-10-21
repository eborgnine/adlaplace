test_that("trend-only smoothing forecasts with increasing uncertainty", {
  y <- fcast(as.numeric(AirPassengers), c(0, 12))
  sm <- univ(y, numeric(), TVP = 1, nvr = 1e-4)

  expect_equal(sm$fit, sm$trend)
  expect_equal(sm$fitse, sm$trendse)
  expect_equal(as.numeric(sm$comp), rep(0, length(y)))
  observed <- is.finite(y)
  expect_equal(as.numeric(sm$y0)[observed], y[observed])
  expect_equal(as.numeric(sm$y0)[!observed], as.numeric(sm$fit)[!observed])
  expect_true(all(diff(tail(as.numeric(sm$trendse), 12)) > 0))
  expect_true(all(is.finite(sm$fitse)))
})

test_that("a joint trend and AR model interpolates and passes through observations", {
  air <- as.numeric(AirPassengers)
  y <- fcast(air, rbind(c(84, 94), c(131, 144)))
  sm <- univ(y, c(-0.4, rep(0, 10), -0.5), TVP = 1, nvr = 1e-3)
  observed <- is.finite(y)

  expect_equal(as.numeric(sm$fit), as.numeric(sm$trend) + as.numeric(sm$comp))
  expect_lt(max(abs(sm$fit[observed] - y[observed])), 1e-6)
  expect_true(all(is.finite(sm$fit[!observed])))
  expect_gt(sm$fitse[144], sm$fitse[131])
  expect_gt(min(sm$fitse[!observed]), max(sm$fitse[observed]))
  expect_equal(ncol(sm$state), 2 + 12)
})

test_that("filtering, random walks, and unsupported options are explicit", {
  y <- as.numeric(AirPassengers)[1:36]
  filtered <- univ(y, numeric(), 0, 0.05, sm = 0)
  expect_equal(length(filtered$fit), length(y))
  expect_true(all(is.finite(filtered$fit)))

  expect_error(univ(y, numeric(), 1, 1e-3, ARt = -0.2), "not implemented")
  expect_error(univ(y, numeric(), 1, 1e-3, Int = c(1, rep(0, length(y) - 1))), "not implemented")
  expect_error(univopt(y, 1:2, tar = 1:2, out = 0), "not implemented")
  expect_error(univ(y, numeric(), TVP = 3, nvr = 1), "TVP")
})
