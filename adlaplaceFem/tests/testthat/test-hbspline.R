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

# Reference THB map built with explicit Kronecker products, i.e. the direct
# transcription of S's definition. hb_basis_map() instead uses
# (Ty %x% Tx) vec(X) = vec(Tx X t(Ty)); this pins that identity down.
naive_hb_basis_map <- function(hb) {
  levs <- hb$levels
  n_levels <- length(levs)
  active <- hb$active
  nx_f <- levs[[n_levels]]$n_basis[["x"]]
  ny_f <- levs[[n_levels]]$n_basis[["y"]]
  refine_x <- adlaplaceFem:::hb_pairwise_refine_1d(hb, "x")
  refine_y <- adlaplaceFem:::hb_pairwise_refine_1d(hb, "y")
  cum_x <- adlaplaceFem:::hb_cumulative_from_pairwise(refine_x, nx_f)
  cum_y <- adlaplaceFem:::hb_cumulative_from_pairwise(refine_y, ny_f)
  cols <- list()
  for (lev in seq_len(n_levels)) {
    act <- active[[lev]]
    for (k in seq_along(act$idx)) {
      i <- act$i[k]
      j <- act$j[k]
      if (lev < n_levels) {
        v <- numeric(act$n_x * act$n_y)
        v[i + (j - 1L) * act$n_x] <- 1
        w <- as.numeric(
          Matrix::kronecker(refine_y[[lev]], refine_x[[lev]]) %*% v
        )
        mask <- numeric(length(w))
        mask[active[[lev + 1L]]$idx] <- 1
        col <- as.numeric(
          Matrix::kronecker(cum_y[[lev + 1L]], cum_x[[lev + 1L]]) %*%
            (w - w * mask)
        )
      } else {
        col <- numeric(nx_f * ny_f)
        col[i + (j - 1L) * nx_f] <- 1
      }
      cols[[length(cols) + 1L]] <- col
    }
  }
  do.call(cbind, cols)
}

test_that("hb_basis_map matches the naive Kronecker construction", {
  for (deg in c(2L, 3L)) {
    kn <- hb_knots(
      unit_box(),
      list(c(0, 1, 0, 1), c(0.25, 0.75, 0.25, 0.75)),
      fact = c(2, 4)
    )
    hb <- hb_basis(kn, degree = deg)
    expect_equal(
      max(abs(as.matrix(hb$S) - naive_hb_basis_map(hb))),
      0,
      tolerance = 1e-10
    )
  }
})

test_that("S is sparse, not dense with roundoff", {
  kn <- hb_knots(
    unit_box(),
    list(c(0, 1, 0, 1), c(0.25, 0.75, 0.25, 0.75)),
    fact = c(2, 4)
  )
  S <- hb_basis(kn, degree = 3L)$S
  expect_lt(Matrix::nnzero(S) / prod(dim(S)), 0.05)
})

test_that("untruncated map represents coarse functions exactly", {
  degree <- 2L
  kn <- hb_knots(unit_box(), c(0.25, 0.75, 0.25, 0.75), fact = 2)
  hb <- hb_basis(kn, degree = degree)
  S <- adlaplaceFem:::hb_basis_map(
    hb, degree,
    truncate = FALSE, active = hb$active
  )
  levs <- hb$levels
  kf <- levs[[length(levs)]]$open_knots
  k0 <- levs[[1L]]$open_knots
  pts <- expand.grid(x = seq(0.05, 0.95, by = 0.1), y = seq(0.05, 0.95, by = 0.1))
  A <- adlaplaceFem:::tensor_design(pts$x, pts$y, kf$x, kf$y, degree)
  Bx <- as.matrix(adlaplaceFem:::bspline_eval(k0$x, pts$x, degree, 0L))
  By <- as.matrix(adlaplaceFem:::bspline_eval(k0$y, pts$y, degree, 0L))
  act <- hb$active[[1L]]
  # Each untruncated column is the level-0 tensor function itself, re-expressed
  # in the finest basis, so evaluating it must reproduce Bx[, i] * By[, j].
  for (k in seq_len(min(8L, length(act$idx)))) {
    expect_equal(
      as.numeric(A %*% S[, k, drop = FALSE]),
      Bx[, act$i[k]] * By[, act$j[k]],
      tolerance = 1e-8
    )
  }
})

test_that("hb_active_sets matches a naive support scan", {
  degree <- 3L
  kn <- hb_knots(
    unit_box(),
    list(c(0, 1, 0, 1), c(0.25, 0.75, 0.25, 0.75)),
    fact = c(2, 2)
  )
  hb <- hb_basis(kn, degree = degree)
  levs <- hb$levels
  inside <- function(sx, sy, omega) {
    any(vapply(omega, function(o) {
      sx[1] >= o$xmin - 1e-6 && sx[2] <= o$xmax + 1e-6 &&
        sy[1] >= o$ymin - 1e-6 && sy[2] <= o$ymax + 1e-6
    }, logical(1L)))
  }
  for (lev in seq_along(levs)) {
    kx <- levs[[lev]]$open_knots$x
    ky <- levs[[lev]]$open_knots$y
    nx <- levs[[lev]]$n_basis[["x"]]
    ny <- levs[[lev]]$n_basis[["y"]]
    om_l <- levs[[lev]]$rasters
    om_n <- if (lev < length(levs)) levs[[lev + 1L]]$rasters else list()
    want <- integer(0)
    for (j in seq_len(ny)) {
      sy <- c(ky[j + 1L], ky[j + degree + 1L])
      for (i in seq_len(nx)) {
        sx <- c(kx[i + 1L], kx[i + degree + 1L])
        in_n <- length(om_n) > 0L && inside(sx, sy, om_n)
        ok <- if (lev == 1L) !in_n else inside(sx, sy, om_l) && !in_n
        if (ok) {
          want <- c(want, i + (j - 1L) * nx)
        }
      }
    }
    expect_equal(hb$active[[lev]]$idx, want)
  }
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
