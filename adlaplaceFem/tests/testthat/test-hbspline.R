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
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.2), y = seq(0.1, 0.9, by = 0.2))
  kn_lines <- adlaplaceFem:::knots_from_spatraster(outer, 2L)
  fem_tensor <- fem_bspline(sites, kn_lines, degree = 2L)
  fem_hb <- fem_bspline(sites, hb_knots(outer), degree = 2L)
  expect_equal(max(abs(fem_tensor$C - fem_hb$C)), 0, tolerance = 1e-9)
  expect_equal(max(abs(fem_tensor$G - fem_hb$G)), 0, tolerance = 1e-9)
})

test_that("list(ext) and list(list(ext)) are the same inner level", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  e <- terra::ext(0.25, 0.75, 0.25, 0.75)
  k1 <- hb_knots(outer, list(e), fact = 2)
  k2 <- hb_knots(outer, list(list(e)), fact = 2)
  k3 <- hb_knots(outer, e, fact = 2)
  expect_equal(k1$n_levels, 2L)
  expect_equal(k1$levels, k2$levels)
  expect_equal(k1$levels, k3$levels)
})

test_that("siblings share a level; two extents are two levels", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  e1 <- terra::ext(0, 0.5, 0, 0.5)
  e2 <- terra::ext(0.5, 1, 0.5, 1)
  sib <- hb_knots(outer, list(list(e1, e2)), fact = 2)
  seqn <- hb_knots(
    outer,
    list(terra::ext(0, 0.75, 0, 0.75), terra::ext(0.25, 0.5, 0.25, 0.5)),
    fact = c(2, 2)
  )
  expect_equal(sib$n_levels, 2L)
  expect_equal(length(sib$levels[[2]]$regions), 2L)
  expect_equal(seqn$n_levels, 3L)
  expect_equal(length(seqn$levels[[2]]$regions), 1L)
  expect_equal(length(seqn$levels[[3]]$regions), 1L)
})

test_that("inner extents snap to outer cell edges", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(outer, terra::ext(0.2, 0.8, 0.2, 0.8), fact = 2)
  reg <- kn$levels[[2]]$regions[[1]]
  expect_equal(reg$xmin, 0)
  expect_equal(reg$xmax, 1)
  expect_equal(reg$ymin, 0)
  expect_equal(reg$ymax, 1)
})

test_that("fact vector sets nested resolutions", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(
    outer,
    list(terra::ext(0, 1, 0, 1), terra::ext(0.25, 0.75, 0.25, 0.75)),
    fact = c(2, 4)
  )
  res0 <- terra::res(outer)
  expect_equal(kn$levels[[1]]$resolution, res0)
  expect_equal(kn$levels[[2]]$resolution, res0 / 2)
  expect_equal(kn$levels[[3]]$resolution, res0 / 8)
  expect_equal(kn$fact, c(2L, 4L))
})

test_that("snapped inner boxes are clipped to the parent", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(
    outer,
    list(terra::ext(0.25, 0.75, 0.25, 0.75), terra::ext(0.2, 0.8, 0.2, 0.8)),
    fact = c(2, 2)
  )
  parent <- kn$levels[[2]]$regions[[1]]
  child <- kn$levels[[3]]$regions[[1]]
  expect_true(child$xmin >= parent$xmin - 1e-9)
  expect_true(child$xmax <= parent$xmax + 1e-9)
  expect_true(child$ymin >= parent$ymin - 1e-9)
  expect_true(child$ymax <= parent$ymax + 1e-9)
})

test_that("refinement outside the parent domain errors", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  expect_error(
    hb_knots(
      outer,
      list(terra::ext(0, 0.5, 0, 0.5), terra::ext(0.5, 1, 0.5, 1)),
      fact = c(2, 2)
    ),
    "inside the parent domain"
  )
})

test_that("legacy raster-level lists are rejected", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  expect_error(
    matern("geometry", knots = list(outer), shape = 1L),
    "hb_knots\\(outer, inner, fact\\)"
  )
})

test_that("hierarchical THB partition of unity holds", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(outer, terra::ext(0.25, 0.75, 0.25, 0.75), fact = 2)
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.1), y = seq(0.1, 0.9, by = 0.1))
  fem <- fem_bspline(sites, kn, degree = 2L)
  rs <- Matrix::rowSums(fem$A)
  expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("hierarchical Q is SPD", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(outer, terra::ext(0.25, 0.75, 0.25, 0.75), fact = 2)
  fem <- fem_bspline(expand.grid(x = 0.5, y = 0.5), kn, degree = 2L)
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
  sites <- expand.grid(
    x = seq(0.05, 0.95, by = 0.05),
    y = seq(0.05, 0.95, by = 0.05)
  )
  outer <- terra::rast(terra::ext(-0.2, 1.2, -0.2, 1.2), resolution = 0.2)
  kn <- hb_knots(outer, terra::ext(0.2, 0.8, 0.2, 0.8), fact = 2)

  fem_cov <- function(fem) {
    Q <- fem_precision(kappa, tau, fem$C, fem$G, fem$G2, alpha = 2L)
    A <- fem$A
    X <- Matrix::solve(Q, Matrix::t(A))
    as.matrix(A %*% X)
  }

  fem_c <- fem_bspline(
    sites,
    adlaplaceFem:::knots_from_spatraster(outer, 2L),
    degree = 2L
  )
  fem_h <- fem_bspline(sites, kn, degree = 2L)
  C_c <- fem_cov(fem_c)
  C_h <- fem_cov(fem_h)

  xy <- as.matrix(sites)
  d <- as.matrix(dist(xy))
  d[d == 0] <- NA
  r <- range / sqrt(8 * nu)
  K <- 2^(1 - nu) / gamma(nu) * (sqrt(2 * nu) * d / r)^nu *
    besselK(sqrt(2 * nu) * d / r, nu)

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

test_that("matern accepts hb_knots", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(outer, terra::ext(0.25, 0.75, 0.25, 0.75), fact = 2)
  term <- matern("geometry", knots = kn, shape = 1L)
  expect_s4_class(term, "matern")
  expect_true(inherits(term@fem$knots, "hb_knots"))
  n_active <- nrow(term@fem$C)
  xy <- cbind(x = c(0.3, 0.5), y = c(0.3, 0.5))
  d <- design(term, data.frame(geometry = I(xy)))
  expect_equal(ncol(d), n_active)
})

test_that("hb_summary reports per-level counts on hb_basis", {
  skip_if_not_installed("terra")
  outer <- terra::rast(terra::ext(0, 1, 0, 1), resolution = 0.25)
  kn <- hb_knots(outer, terra::ext(0.25, 0.75, 0.25, 0.75), fact = 2)
  sm <- hb_summary(hb_basis(kn, degree = 2L))
  expect_true(all(c("level", "n_active", "total_active") %in% names(sm)))
  expect_equal(sum(sm$n_active), sm$total_active[1])
})
