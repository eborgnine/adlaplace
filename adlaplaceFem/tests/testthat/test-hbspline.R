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

unit_box <- function() {
  list(xmin = 0, xmax = 1, ymin = 0, ymax = 1, resolution = 0.25)
}

test_that("hb_basis works for fact-4 local refinement", {
  kn <- hb_knots(
    unit_box(),
    list(c(0, 1, 0, 1), c(0.25, 0.75, 0.25, 0.75)),
    fact = 4
  )
  hb <- hb_basis(kn, degree = 3L)
  expect_true(inherits(hb, "hb_basis"))
  expect_gt(nrow(hb$S), 0)
})

test_that("single-level hierarchy matches tensor product", {
  outer <- unit_box()
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.2), y = seq(0.1, 0.9, by = 0.2))
  kn_lines <- adlaplaceFem:::knots_from_spec(outer)
  fem_tensor <- fem_bspline(sites, kn_lines, degree = 2L)
  fem_hb <- fem_bspline(sites, hb_knots(outer), degree = 2L)
  expect_equal(max(abs(fem_tensor$C - fem_hb$C)), 0, tolerance = 1e-9)
  expect_equal(max(abs(fem_tensor$G - fem_hb$G)), 0, tolerance = 1e-9)
})

test_that("list(ext) and list(list(ext)) are the same inner level", {
  outer <- unit_box()
  e <- c(0.25, 0.75, 0.25, 0.75)
  k1 <- hb_knots(outer, list(e), fact = 2)
  k2 <- hb_knots(outer, list(list(e)), fact = 2)
  k3 <- hb_knots(outer, e, fact = 2)
  k4 <- hb_knots(
    outer,
    list(xmin = 0.25, xmax = 0.75, ymin = 0.25, ymax = 0.75),
    fact = 2
  )
  expect_equal(length(k1), 2L)
  expect_equal(unclass(k1), unclass(k2))
  expect_equal(unclass(k1), unclass(k3))
  expect_equal(unclass(k1), unclass(k4))
})

test_that("siblings share a level; two extents are two levels", {
  outer <- unit_box()
  e1 <- c(0, 0.5, 0, 0.5)
  e2 <- c(0.5, 1, 0.5, 1)
  sib <- hb_knots(outer, list(list(e1, e2)), fact = 2)
  seqn <- hb_knots(
    outer,
    list(c(0, 0.75, 0, 0.75), c(0.25, 0.5, 0.25, 0.5)),
    fact = c(2, 2)
  )
  expect_equal(length(sib), 2L)
  expect_equal(length(sib[[2]]$rasters), 2L)
  expect_equal(length(seqn), 3L)
  expect_equal(length(seqn[[2]]$rasters), 1L)
  expect_equal(length(seqn[[3]]$rasters), 1L)
})

test_that("inner extents snap to that level's nested grid", {
  kn <- hb_knots(unit_box(), c(0.2, 0.8, 0.2, 0.8), fact = 2)
  reg <- kn[[2]]$rasters[[1]]
  expect_equal(reg$xmin, 0.125)
  expect_equal(reg$xmax, 0.875)
  expect_equal(reg$ymin, 0.125)
  expect_equal(reg$ymax, 0.875)
})

test_that("later levels keep tight boxes on the fine grid", {
  kn <- hb_knots(
    unit_box(),
    list(c(0, 1, 0, 1), c(0.3, 0.45, 0.3, 0.45)),
    fact = c(2, 4)
  )
  child <- kn[[3]]$rasters[[1]]
  expect_lte(child$xmax - child$xmin, 0.25)
  expect_lte(child$xmin, 0.3)
  expect_gte(child$xmax, 0.45)
})

test_that("fact vector sets nested resolutions", {
  res0 <- c(0.25, 0.25)
  kn <- hb_knots(
    unit_box(),
    list(c(0, 1, 0, 1), c(0.25, 0.75, 0.25, 0.75)),
    fact = c(2, 4)
  )
  expect_equal(kn[[1]]$rasters[[1]]$resolution, res0)
  expect_equal(kn[[2]]$rasters[[1]]$resolution, res0 / 2)
  expect_equal(kn[[3]]$rasters[[1]]$resolution, res0 / 8)
})

test_that("snapped inner boxes are clipped to the parent", {
  kn <- hb_knots(
    unit_box(),
    list(c(0.25, 0.75, 0.25, 0.75), c(0.2, 0.8, 0.2, 0.8)),
    fact = c(2, 2)
  )
  parent <- kn[[2]]$rasters[[1]]
  child <- kn[[3]]$rasters[[1]]
  expect_true(child$xmin >= parent$xmin - 1e-9)
  expect_true(child$xmax <= parent$xmax + 1e-9)
  expect_true(child$ymin >= parent$ymin - 1e-9)
  expect_true(child$ymax <= parent$ymax + 1e-9)
})

test_that("refinement outside the parent domain errors", {
  expect_error(
    hb_knots(
      unit_box(),
      list(c(0, 0.5, 0, 0.5), c(0.5, 1, 0.5, 1)),
      fact = c(2, 2)
    ),
    "inside the parent domain"
  )
})

test_that("SpatRaster outer matches raster-spec list", {
  skip_if_not_installed("terra")
  spec <- unit_box()
  r <- do.call(terra::rast, spec)
  e <- terra::ext(0.25, 0.75, 0.25, 0.75)
  expect_equal(
    unclass(hb_knots(r, e, fact = 2)),
    unclass(hb_knots(spec, c(0.25, 0.75, 0.25, 0.75), fact = 2))
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
  kn <- hb_knots(unit_box(), c(0.25, 0.75, 0.25, 0.75), fact = 2)
  sites <- expand.grid(x = seq(0.1, 0.9, by = 0.1), y = seq(0.1, 0.9, by = 0.1))
  fem <- fem_bspline(sites, kn, degree = 2L)
  rs <- Matrix::rowSums(fem$A)
  expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("hierarchical Q is SPD", {
  kn <- hb_knots(unit_box(), c(0.25, 0.75, 0.25, 0.75), fact = 2)
  fem <- fem_bspline(expand.grid(x = 0.5, y = 0.5), kn, degree = 2L)
  Q <- fem_precision(2, 1, fem$C, fem$G, fem$G2, alpha = 2L)
  ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 1e-10))
})

test_that("hierarchical basis improves Matern correlation fit", {
  kappa <- 3
  tau <- 1
  nu <- 1
  range <- sqrt(8 * nu) / kappa
  sites <- expand.grid(
    x = seq(0.05, 0.95, by = 0.05),
    y = seq(0.05, 0.95, by = 0.05)
  )
  outer <- list(xmin = -0.2, xmax = 1.2, ymin = -0.2, ymax = 1.2, resolution = 0.2)
  kn <- hb_knots(outer, c(0.2, 0.8, 0.2, 0.8), fact = 2)

  fem_cov <- function(fem) {
    Q <- fem_precision(kappa, tau, fem$C, fem$G, fem$G2, alpha = 2L)
    A <- fem$A
    X <- Matrix::solve(Q, Matrix::t(A))
    as.matrix(A %*% X)
  }

  fem_c <- fem_bspline(
    sites,
    adlaplaceFem:::knots_from_spec(outer),
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
  fine_full <- list(
    xmin = -0.2, xmax = 1.2, ymin = -0.2, ymax = 1.2, resolution = 0.1
  )
  fem_f <- fem_bspline(
    sites,
    adlaplaceFem:::knots_from_spec(fine_full),
    degree = 2L
  )
  expect_lt(ncol(fem_h$A), ncol(fem_f$A))
})

test_that("matern accepts hb_knots", {
  kn <- hb_knots(unit_box(), c(0.25, 0.75, 0.25, 0.75), fact = 2)
  term <- matern("geometry", knots = kn, shape = 1L)
  expect_s4_class(term, "matern")
  expect_true(inherits(term@fem$knots, "hb_knots"))
  n_active <- nrow(term@fem$C)
  xy <- cbind(x = c(0.3, 0.5), y = c(0.3, 0.5))
  d <- design(term, data.frame(geometry = I(xy)))
  expect_equal(ncol(d), n_active)
})

test_that("hb_summary reports per-level counts on hb_basis", {
  kn <- hb_knots(unit_box(), c(0.25, 0.75, 0.25, 0.75), fact = 2)
  sm <- hb_summary(hb_basis(kn, degree = 2L))
  expect_true(all(c("level", "n_active", "total_active") %in% names(sm)))
  expect_equal(sum(sm$n_active), sm$total_active[1])
})
