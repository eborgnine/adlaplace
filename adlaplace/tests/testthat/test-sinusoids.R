test_that("sinusoid_basis matches cos and sin at the given harmonics", {
  period <- 365.25
  x <- c(0, 10, 40, 180, 365.25)
  harmonics <- c(1L, 2L)
  D <- sinusoid_basis(x, period = period, harmonics = harmonics)
  u <- x
  expect_equal(D[, 1], cos(2 * pi * 1 * u / period))
  expect_equal(D[, 2], sin(2 * pi * 1 * u / period))
  expect_equal(D[, 3], cos(2 * pi * 2 * u / period))
  expect_equal(D[, 4], sin(2 * pi * 2 * u / period))
  expect_null(colnames(D))

  shifted <- sinusoid_basis(x, period = period, harmonics = 1L, ref_value = 10)
  expect_equal(shifted[, 1], cos(2 * pi * (x - 10) / period))
})

test_that("sinusoids design scales by the exposure column", {
  period <- 365.25
  dat <- data.frame(
    day = c(0, 30, 60, 90),
    exposure = c(1, 2, 0, 4)
  )
  plain <- sinusoids("day", period = period, harmonics = c(1, 2))
  expect_identical(plain@label, "day_sinusoids")
  expect_true(is.na(plain@mult))
  D <- design(plain, dat)
  expect_identical(
    colnames(D),
    c("day_sinusoids_c1", "day_sinusoids_s1", "day_sinusoids_c2", "day_sinusoids_s2")
  )
  basis <- sinusoid_basis(dat$day, period = period, harmonics = c(1, 2))
  expect_equal(unname(D), basis)

  scaled <- sinusoids(
    "day",
    period = period,
    harmonics = c(1, 2),
    mult = "exposure",
    ref_mult = 1
  )
  expect_identical(scaled@label, "day_exposure_sinusoids")
  S <- design(scaled, dat)
  expect_identical(
    colnames(S),
    c(
      "day_exposure_sinusoids_c1", "day_exposure_sinusoids_s1",
      "day_exposure_sinusoids_c2", "day_exposure_sinusoids_s2"
    )
  )
  expect_equal(unname(S), basis * (dat$exposure - 1))

  info <- beta_info(scaled, dat)
  expect_identical(info$beta_label, colnames(S))
  expect_equal(nrow(info), 4L)
  expect_identical(info$model, rep("sinusoids", 4L))
  expect_null(precision(scaled, dat))
  expect_null(theta_info(scaled))
  expect_null(random_info(scaled, dat))
})

test_that("collect_terms treats mult as a column name", {
  terms <- collect_terms(
    y ~ sinusoids(day, period = 365.25, harmonics = c(1, 2), mult = exposure)
  )
  term <- terms$sinusoids_day_exposure
  expect_true(inherits(term, "sinusoids"))
  expect_identical(term@name, "day")
  expect_identical(term@mult, "exposure")
  expect_identical(term@harmonics, c(1L, 2L))
  expect_equal(term@period, 365.25)
})

test_that("replace_vars rewrites the sinusoid covariate and the multiplier", {
  term <- sinusoids("day", harmonics = c(1, 4), mult = "exposure")
  out <- replace_vars(term, c(day = "day_s1", exposure = "exposure_s1"))
  expect_identical(out@name, "day_s1")
  expect_identical(out@mult, "exposure_s1")
  expect_identical(out@label, "day_s1_exposure_s1_sinusoids")
  expect_identical(out@harmonics, c(1L, 4L))

  only_mult <- replace_vars(term, c(exposure = "pm"))
  expect_identical(only_mult@name, "day")
  expect_identical(only_mult@mult, "pm")
  expect_identical(only_mult@label, "day_pm_sinusoids")
})

test_that("sinusoids rejects a bad period or harmonic", {
  expect_error(sinusoids("day", period = 0), "period")
  expect_error(sinusoids("day", harmonics = c(1, 1)), "unique")
  expect_error(sinusoids("day", harmonics = 0), "positive")
  expect_error(sinusoid_basis(c(0, NA), period = 1), "finite")
})
