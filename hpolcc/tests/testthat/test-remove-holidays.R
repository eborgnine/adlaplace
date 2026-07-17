test_that("remove_holidays rm_none returns data unchanged", {
  dat <- data.frame(
    date = as.Date("2002-06-15"),
    y = 1
  )
  out <- hpolcc::remove_holidays(dat, type = "rm_none")
  expect_equal(out, dat)
})

test_that("remove_holidays stops when no date column", {
  dat <- data.frame(year = 2002L, y = 1)
  expect_error(
    hpolcc::remove_holidays(dat, type = "rm_all"),
    "no date column found",
    fixed = TRUE
  )
})

test_that("remove_holidays rm_all filters holidays and school-year window", {
  skip_if_not_installed("timeDate")
  dat <- data.frame(
    date = as.Date(c(
      "2002-01-05", # before school-year window
      "2002-07-01", # Canada Day
      "2002-06-15"  # ordinary weekday in window
    )),
    y = 1:3
  )
  out <- hpolcc::remove_holidays(dat, type = "rm_all")
  expect_false(any(out$date == as.Date("2002-01-05")))
  expect_false(any(out$date == as.Date("2002-07-01")))
  expect_true(any(out$date == as.Date("2002-06-15")))
})
