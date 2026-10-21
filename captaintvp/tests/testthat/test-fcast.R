test_that("fcast masks inclusive intervals and extends the series", {
  y <- 1:30
  x <- fcast(y, rbind(c(14, 19), c(22, 22), c(0, 10), c(5, 0)))

  expect_equal(length(x), 45)
  expect_true(all(is.na(x[c(1:5, 19:24, 27, 36:45)])))
  expect_equal(x[6:18], 1:13)
  expect_equal(x[25:26], 20:21)
  expect_equal(x[28:35], 23:30)

  forecast <- fcast(1:5, 4)
  expect_equal(forecast, c(1:5, rep(NA_real_, 4)))
  expect_equal(fcast(1:5, c(0, 4)), forecast)
  expect_equal(fcast(1:5), c(1:5, rep(NA_real_, 10)))

  backcast <- fcast(1:5, c(3, 0))
  expect_equal(backcast, c(rep(NA_real_, 3), 1:5))
  expect_error(fcast(1:10, c(6, 4)), "p <= q")
})

test_that("fcast preserves and shifts ts time bases", {
  y <- ts(1:24, start = c(1949, 1), frequency = 12)
  masked <- fcast(y, c(5, 6))
  expect_equal(start(masked), c(1949, 1))
  expect_equal(frequency(masked), 12)
  expect_true(all(is.na(masked[5:6])))
  expect_equal(as.numeric(masked[-c(5, 6)]), as.numeric(y[-c(5, 6)]))

  extended <- fcast(y, c(2, 0))
  expect_equal(start(extended), c(1948, 11))
  expect_equal(as.numeric(extended[3:26]), as.numeric(y))
})

test_that("the air-demo mask removes samples 84:94 and 131:144", {
  air <- as.numeric(AirPassengers)
  y <- fcast(air, rbind(c(84, 94), c(131, 144)))
  missing <- which(is.na(y))
  expect_equal(missing, c(84:94, 131:144))
  expect_equal(y[-missing], air[-missing])
  expect_equal(length(y), length(air))
})
