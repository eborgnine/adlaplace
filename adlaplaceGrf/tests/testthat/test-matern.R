test_that("crds works for matrix and data.frame without terra", {
  m <- cbind(1:3, 4:6)
  expect_equal(crds(m)[, 1], 1:3)
  df <- data.frame(x = 1:2, y = 3:4)
  expect_equal(crds(df)[, 2], 3:4)
})

test_that("matern term builds design and list precision", {
  pts <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4)))
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  term <- matern(x = pts, knots = knots_list, degree = 2L)
  expect_s4_class(term, "matern")
  expect_identical(term@ad_fun, "random_fem_2")
  dat <- data.frame(dummy = seq_len(nrow(pts)))
  A <- design(term, dat)
  expect_equal(nrow(A), nrow(pts))
  prec <- precision(term, dat)
  expect_true(is.list(prec))
  expect_equal(prec$alpha, 2L)
  th <- theta_info(term)
  expect_equal(nrow(th), 2L)
  expect_true(all(th$transform))
})

test_that("matern degree 3 selects random_fem_3", {
  pts <- cbind(runif(9), runif(9))
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  term <- matern(x = pts, knots = knots_list, degree = 3L)
  expect_identical(term@ad_fun, "random_fem_3")
  expect_equal(term@alpha, 3L)
})

test_that("crds rejects SpatVector polygons when terra is available", {
  skip_if_not_installed("terra")
  v <- terra::vect(cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)), type = "polygons")
  expect_error(crds(v), "points only")
})

test_that("matern formula term works via :: and model_data", {
  set.seed(2)
  n <- 16L
  pts <- cbind(runif(n), runif(n))
  colnames(pts) <- c("x", "y")
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  dat <- data.frame(y = rnorm(n), z = rnorm(n))
  f <- y ~ z + adlaplaceGrf::matern(x = pts, knots = knots_list, degree = 2L)
  environment(f) <- environment()
  terms <- adlaplace::collect_terms(f)
  mat <- terms[[grep("matern", names(terms))[1L]]]
  expect_s4_class(mat, "matern")
  expect_identical(mat@ad_fun, "random_fem_2")

  md <- adlaplace::model_data(f, data = dat)
  expect_length(md$random, 1L)
  prec <- md$random[[1L]]@precision
  expect_true(is.list(prec))
  expect_true(all(c("Q_p", "Q_i", "C_x", "chol", "alpha") %in% names(prec)))
  expect_identical(md$random[[1L]]@ad_fun, "random_fem_2")
  expect_equal(ncol(md$random[[1L]]@theta_map), 2L)
})
