test_that("refine_matrix_1d preserves coarse basis values", {
  kn0x <- adlaplaceFem:::axis_to_open_knots(seq(0, 1, length.out = 5), 2L)
  kn1x <- adlaplaceFem:::axis_to_open_knots(seq(0, 1, length.out = 9), 2L)
  R <- adlaplaceFem:::refine_matrix_1d(kn0x, kn1x, 2L)
  xs <- seq(0, 1, length.out = 20)
  err <- vapply(xs, function(x) {
    Bf <- as.matrix(adlaplaceFem:::bspline_eval(kn1x, x, 2L, 0L))
    Bc <- as.matrix(adlaplaceFem:::bspline_eval(kn0x, x, 2L, 0L))
    max(abs(as.numeric(Bf %*% R - Bc)))
  }, numeric(1))
  expect_true(max(err) < 1e-10)
})

test_that("single-level hierarchy matches tensor product", {
  skip_if_not_installed("terra")
  coarse <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.2), y = seq(0.1, 0.9, by = 0.2))
  kn <- adlaplaceFem:::knots_from_spatraster(coarse, 2L)
  fem_tensor <- fem_bspline(sites, kn, degree = 2L)
  fem_hb <- fem_bspline(sites, list(coarse), degree = 2L)
  expect_equal(max(abs(fem_tensor$C - fem_hb$C)), 0, tolerance = 1e-9)
  expect_equal(max(abs(fem_tensor$G - fem_hb$G)), 0, tolerance = 1e-9)
})

test_that("hierarchical THB partition of unity holds", {
  skip_if_not_installed("terra")
  coarse <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  fine <- terra::rast(terra::ext(0.2, 0.8, 0.2, 0.8), resolution = 0.1)
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.1), y = seq(0.1, 0.9, by = 0.1))
  fem <- fem_bspline(sites, list(coarse, fine), degree = 2L)
  rs <- Matrix::rowSums(fem$A)
  expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("hierarchical Q is SPD", {
  skip_if_not_installed("terra")
  coarse <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  fine <- terra::rast(terra::ext(0.25, 0.75, 0.25, 0.75), resolution = 0.1)
  fem <- fem_bspline(
    expand.grid(x = 0.5, y = 0.5),
    list(coarse, fine),
    degree = 2L
  )
  Q <- fem_precision(2, 1, fem$C, fem$G, fem$G2, alpha = 2L)
  ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 1e-10))
})

test_that("hierarchical basis improves Matern correlation fit", {
  skip_if_not_installed("terra")
  kappa <- 3
  tau <- 1
  nu <- 1
  range <- sqrt(8 * nu) / kappa
  sites <- expand.grid(x = seq(0.05, 0.95, by = 0.05), y = seq(0.05, 0.95, by = 0.05))
  coarse <- terra::rast(terra::ext(-0.2, 1.2, -0.2, 1.2), resolution = 0.2)
  fine <- terra::rast(terra::ext(0.2, 0.8, 0.2, 0.8), resolution = 0.1)

  fem_cov <- function(fem) {
    Q <- fem_precision(kappa, tau, fem$C, fem$G, fem$G2, alpha = 2L)
    A <- fem$A
    X <- Matrix::solve(Q, Matrix::t(A))
    as.matrix(A %*% X)
  }

  fem_c <- fem_bspline(sites, adlaplaceFem:::knots_from_spatraster(coarse, 2L), degree = 2L)
  fem_h <- fem_bspline(sites, list(coarse, fine), degree = 2L)
  C_c <- fem_cov(fem_c)
  C_h <- fem_cov(fem_h)

  xy <- as.matrix(sites)
  d <- as.matrix(dist(xy))
  d[d == 0] <- NA
  r <- range / sqrt(8 * nu)
  K <- 2^(1 - nu) / gamma(nu) * (sqrt(2 * nu) * d / r)^nu * besselK(sqrt(2 * nu) * d / r, nu)

  err_c <- mean(abs(C_c[lower.tri(C_c)] - K[lower.tri(K)]), na.rm = TRUE)
  err_h <- mean(abs(C_h[lower.tri(C_h)] - K[lower.tri(K)]), na.rm = TRUE)
  expect_lte(err_h, err_c + 1e-6)
  fine_full <- terra::rast(terra::ext(-0.2, 1.2, -0.2, 1.2), resolution = 0.1)
  fem_f <- fem_bspline(
    sites,
    adlaplaceFem:::knots_from_spatraster(fine_full, 2L),
    degree = 2L
  )
  expect_lt(ncol(fem_h$A), ncol(fem_f$A))
})

test_that("matern accepts hierarchical knots list", {
  skip_if_not_installed("terra")
  coarse <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  fine <- terra::rast(terra::ext(0.2, 0.8, 0.2, 0.8), resolution = 0.1)
  term <- matern("geometry", knots = list(coarse, fine), shape = 1L)
  expect_s4_class(term, "matern")
  expect_true(inherits(term@fem$knots, "hb_knots"))
  n_active <- nrow(term@fem$C)
  xy <- cbind(x = c(0.3, 0.5), y = c(0.3, 0.5))
  d <- design(term, data.frame(geometry = I(xy)))
  expect_equal(ncol(d), n_active)
})

test_that("hb_summary reports per-level counts", {
  skip_if_not_installed("terra")
  coarse <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  fine <- terra::rast(terra::ext(0.2, 0.8, 0.2, 0.8), resolution = 0.1)
  hb <- hb_knots(list(coarse, fine), degree = 2L)
  sm <- hb_summary(hb)
  expect_true(all(c("level", "n_active", "total_active") %in% names(sm)))
  expect_equal(sum(sm$n_active), sm$total_active[1])
})
