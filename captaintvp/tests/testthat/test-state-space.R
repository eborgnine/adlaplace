test_that("IRW plus AR state matrices have the CAPTAIN structure", {
  arp <- c(-0.5, 0.25)
  system <- captaintvp:::trend_ar_system(arp, TVP = 1, nvr = 0.1)

  expect_equal(system$T[1:2, 1:2], rbind(c(1, 1), c(0, 1)))
  expect_equal(system$Z[1, 1:4], c(1, 0, 1, 0))
  expect_equal(system$T[3, 3:4], -arp)
  expect_equal(system$T[4, 3:4], c(1, 0))
  expect_equal(diag(system$Q), c(0.1, 1))
  expect_equal(system$H, 0)
  expect_true(all(eigen(system$P1, symmetric = TRUE, only.values = TRUE)$values > 0))
  expect_equal(system$names, c("level", "slope", "ar1", "ar2"))
})

test_that("a stationary AR(1) block uses its marginal variance", {
  system <- captaintvp:::trend_ar_system(0.5, TVP = 1, nvr = 0.2)
  # ARp = 0.5 is the polynomial coefficient, so the pole is -0.5.
  expect_equal(system$P1[3, 3], 1 / (1 - 0.5^2), tolerance = 1e-8)
  expect_equal(diag(system$P1)[1:2], rep(captaintvp:::CAPTAIN_P0, 2))
})

test_that("random-walk trend-only and diffuse initializations are positive", {
  rw <- captaintvp:::trend_ar_system(NULL, TVP = 0, nvr = 0.3)
  expect_equal(rw$T, matrix(1))
  expect_equal(rw$H, 1)
  expect_equal(rw$P1, matrix(CAPTAIN_P0))
  expect_equal(rw$P1inf, matrix(0))

  diffuse <- captaintvp:::trend_ar_system(NULL, TVP = 1, nvr = 1e-4, diffuse = TRUE)
  expect_equal(diffuse$P1inf, diag(2))
  expect_equal(diffuse$P1, diag(0, 2))
})
