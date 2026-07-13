test_that("grf_bspline returns sparse A, C, G, G2 for degree 2", {
  knots_list <- list(x = seq(0, 1, length.out = 6), y = seq(0, 1, length.out = 6))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- grf_bspline(sites_eval, knots_list, degree = 2L)

  expect_s4_class(fem$A, "dgCMatrix")
  expect_s4_class(fem$C, "dgCMatrix")
  expect_s4_class(fem$G, "dgCMatrix")
  expect_s4_class(fem$G2, "dgCMatrix")
  expect_null(fem$G3)
  expect_equal(nrow(fem$A), nrow(sites_eval))
  expect_equal(ncol(fem$A), prod(fem$n_basis))
  expect_equal(dim(fem$C), rep(prod(fem$n_basis), 2))

  # Partition of unity in the interior for B-splines
  rs <- Matrix::rowSums(fem$A)
  expect_true(all(abs(rs - 1) < 1e-8))
})

test_that("grf_bspline accepts positional knots_list", {
  sites_list <- list(x = seq(-0.2, 1.2, by = 0.2), y = seq(-0.2, 1.2, by = 0.2))
  sites_eval <- do.call(expand.grid, sites_list)
  knots_list <- list(x = seq(-0.2, 1.2, length.out = 5), y = seq(-0.2, 1.2, length.out = 5))
  fem <- grf_bspline(sites_eval, knots_list, degree = 2L)
  # 5 knot lines -> 3 interior; open vector length = 3 + 2*(degree+1) = 9; n_basis = 9 - 3 = 6
  expect_equal(unname(fem$n_basis), c(6L, 6L))
  expect_equal(length(unique(fem$knots$x)), 5L)
})

test_that("grf_bspline returns G3 for degree 3", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- grf_bspline(sites_eval, knots_list, degree = 3L)
  expect_s4_class(fem$G3, "dgCMatrix")
  expect_equal(dim(fem$G3), dim(fem$C))
})

test_that("fem_precision Q2 is SPD on a small grid", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- grf_bspline(sites_eval, knots_list, degree = 2L)
  Q <- fem_precision(kappa = 2, tau = 1, fem$C, fem$G, fem$G2, alpha = 2L)
  ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 1e-10))
})

test_that("fem_chol_pattern returns perm and L1", {
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- grf_bspline(sites_eval, knots_list, degree = 2L)
  payload <- fem_precision_payload(fem, alpha = 2L)
  expect_true(!is.null(payload$chol$perm))
  expect_equal(length(payload$chol$perm), nrow(fem$C))
  expect_true(inherits(payload$chol$L1, "Matrix"))
})
