test_that("matern builds design from a matrix geometry column", {
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  n <- 9L
  geom <- cbind(runif(n), runif(n))
  dat <- data.frame(y = rnorm(n), geometry = I(geom))
  # I(matrix) stores as list-like; use explicit matrix column via structure
  dat$geometry <- geom
  term <- matern("geometry", knots = knots_list)
  expect_s4_class(term, "matern")
  expect_equal(term@term, "geometry")
  expect_equal(term@p.order, 2L)
  expect_null(term@fem$A)
  A <- design(term, dat)
  expect_equal(nrow(A), n)
  expect_equal(ncol(A), nrow(term@fem$C))
})

test_that("matern parses WKT and HEX point columns", {
  skip_if_not_installed("terra")
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  pts <- cbind(c(0.2, 0.5, 0.8), c(0.3, 0.4, 0.7))
  v <- terra::vect(pts, type = "points")
  wkt <- as.data.frame(v, geom = "WKT")$geometry
  hex <- as.data.frame(v, geom = "HEX")$geometry

  dat_wkt <- data.frame(y = rnorm(3), geometry = wkt)
  term <- matern("geometry", knots = knots_list, shape = 1L)
  A_wkt <- design(term, dat_wkt)
  expect_equal(nrow(A_wkt), 3L)

  dat_hex <- data.frame(y = rnorm(3), geometry = hex)
  A_hex <- design(term, dat_hex)
  expect_equal(as.matrix(A_hex), as.matrix(A_wkt))
})

test_that("matern shape 2 selects random_fem_ssq_3", {
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  term <- matern("geometry", knots = knots_list, shape = 2L)
  expect_identical(term@ad_fun, "random_fem_ssq_3")
  expect_identical(extra_ad_fun(term), "random_fem_det_3")
  expect_identical(term@package, "adlaplaceGrf")
  expect_equal(term@p.order, 3L)
  expect_equal(precision(term, data.frame())$alpha, 3L)
})

test_that("matern accepts SpatRaster knots via knots_from_spatraster", {
  skip_if_not_installed("terra")
  set.seed(3)
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  term <- matern("geometry", knots = r)
  expect_s4_class(term, "matern")
  kn <- adlaplaceGrf:::knots_from_spatraster(r, degree = 2L)
  expect_equal(sort(unique(term@fem$knots$x)), sort(unique(kn$x)))
  expect_equal(sort(unique(term@fem$knots$y)), sort(unique(kn$y)))
  expect_equal(length(unique(kn$x)), terra::ncol(r) + 2L)
  expect_equal(length(unique(kn$y)), terra::nrow(r) + 2L)
  expect_null(term@fem$A)
  expect_silent(precision(term, data.frame()))
})

test_that("matern formula term works via :: and model_data", {
  set.seed(2)
  n <- 16L
  geom <- cbind(runif(n), runif(n))
  colnames(geom) <- c("x", "y")
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  dat <- data.frame(y = rnorm(n), z = rnorm(n), geometry = I(geom))
  dat$geometry <- geom
  f <- y ~ z + adlaplaceGrf::matern(geometry, knots = knots_list)
  environment(f) <- environment()
  terms <- adlaplace::collect_terms(f)
  mat <- terms[[grep("matern", names(terms))[1L]]]
  expect_s4_class(mat, "matern")
  expect_identical(mat@ad_fun, "random_fem_ssq_2")
  expect_identical(extra_ad_fun(mat), "random_fem_det_2")
  expect_equal(mat@term, "geometry")

  md <- adlaplace::model_data(f, data = dat)
  expect_length(md$random, 1L)
  prec <- md$random[[1L]]@precision
  expect_true(is.list(prec))
  expect_true(all(c("Q_p", "Q_i", "C_x", "chol", "alpha") %in% names(prec)))
  expect_identical(md$random[[1L]]@ad_fun, "random_fem_ssq_2")
  expect_identical(md$random[[1L]]@package, "adlaplaceGrf")
  expect_equal(ncol(md$random[[1L]]@theta_map), 2L)

  det_idx <- grep("_det$", names(md$parameters))
  expect_length(det_idx, 1L)
  det <- md$parameters[[det_idx]]
  expect_identical(det@ad_fun, "random_fem_det_2")
  expect_identical(det@ad_kind, "parameters")
  expect_identical(det@package, "adlaplaceGrf")
  expect_identical(det@precision, prec)
  expect_equal(ncol(det@theta_map), 2L)
})

test_that("matern_est returns mean/sd and optional sims on SpatRaster", {
  skip_if_not_installed("terra")
  set.seed(4)
  n <- 16L
  geom <- cbind(runif(n, 0.1, 0.9), runif(n, 0.1, 0.9))
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  dat <- data.frame(y = rnorm(n), z = rnorm(n))
  dat$geometry <- geom
  fit <- adlaplace::adlaplace(
    y ~ z + matern(geometry, knots = knots_list),
    data = dat,
    control = list(maxit = 40L)
  )
  eval_grid <- terra::rast(
    nrows = 5, ncols = 5, xmin = 0, xmax = 1, ymin = 0, ymax = 1
  )

  out0 <- matern_est(fit, eval_grid, n = 0L)
  expect_s4_class(out0, "SpatRaster")
  expect_equal(terra::nlyr(out0), 2L)
  expect_equal(names(out0), c("mean", "sd"))
  expect_true(all(terra::values(out0$sd) >= 0, na.rm = TRUE))

  out2 <- matern_est(fit, eval_grid, n = 2L)
  expect_equal(terra::nlyr(out2), 4L)
  expect_equal(names(out2), c("mean", "sd", "sim1", "sim2"))
  expect_equal(terra::values(out2$mean), terra::values(out0$mean))
})
