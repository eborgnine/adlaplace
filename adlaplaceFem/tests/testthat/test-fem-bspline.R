test_that("fem_bspline returns sparse A, C, G, G2 for degree 2", {
  knots_list <- list(x = seq(0, 1, length.out = 6), y = seq(0, 1, length.out = 6))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(sites_eval, knots_list, degree = 2L)

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

test_that("fem_bspline accepts positional knots_list", {
  sites_list <- list(x = seq(-0.2, 1.2, by = 0.2), y = seq(-0.2, 1.2, by = 0.2))
  sites_eval <- do.call(expand.grid, sites_list)
  knots_list <- list(x = seq(-0.2, 1.2, length.out = 5), y = seq(-0.2, 1.2, length.out = 5))
  fem <- fem_bspline(sites_eval, knots_list, degree = 2L)
  # 5 knot lines -> 3 interior; open vector length = 3 + 2*(degree+1) = 9; n_basis = 9 - 3 = 6
  expect_equal(unname(fem$n_basis), c(6L, 6L))
  expect_equal(length(unique(fem$knots$x)), 5L)
})

test_that("fem_bspline returns G3 for degree 3", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(sites_eval, knots_list, degree = 3L)
  expect_s4_class(fem$G3, "dgCMatrix")
  expect_equal(dim(fem$G3), dim(fem$C))
})

test_that("fem_precision Q2 is SPD on a small grid", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(sites_eval, knots_list, degree = 2L)
  Q <- fem_precision(kappa = 2, tau = 1, fem$C, fem$G, fem$G2, alpha = 2L)
  ev <- eigen(as.matrix(Q), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 1e-10))
})

test_that("fem_chol_pattern returns perm and L1", {
  knots_list <- list(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4))
  sites_eval <- do.call(expand.grid, knots_list)
  fem <- fem_bspline(sites_eval, knots_list, degree = 2L)
  payload <- fem_precision_payload(fem, alpha = 2L)
  expect_true(!is.null(payload$chol$perm))
  expect_equal(length(payload$chol$perm), nrow(fem$C))
  expect_true(inherits(payload$chol$L1, "Matrix"))
})

test_that("align_gram_to_pattern skips empty CSC columns", {
  # Upper-tri CSC with an empty middle column (p has a repeated pointer)
  p <- as.integer(c(0, 2, 2, 3))
  i <- as.integer(c(0, 1, 2))
  M <- Matrix::Diagonal(3)
  x <- adlaplaceFem:::align_gram_to_pattern(M, p, i, n = 3L)
  expect_equal(length(x), 3L)
  expect_true(all(is.finite(x)))
})

# Literal transcription of "read M at each pattern position", the definition
# align_gram_to_pattern() is a vectorized form of.
reference_align <- function(M, p, i, n) {
  M <- methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
  x <- numeric(length(i))
  for (col in seq_len(n) - 1L) {
    from <- p[col + 1L]
    to <- p[col + 2L] - 1L
    if (from > to) {
      next
    }
    for (pos in from:to) {
      x[pos + 1L] <- M[i[pos + 1L] + 1L, col + 1L]
    }
  }
  x
}

test_that("align_gram_to_pattern matches a per-entry reference read", {
  set.seed(42)
  n <- 25L
  A <- Matrix::rsparsematrix(n, n, density = 0.15)
  M <- methods::as(
    methods::as(Matrix::forceSymmetric(A), "generalMatrix"),
    "CsparseMatrix"
  )
  pat <- adlaplaceFem:::upper_csc_pattern(M)
  expect_equal(
    adlaplaceFem:::align_gram_to_pattern(M, pat$p, pat$i, pat$n),
    reference_align(M, pat$p, pat$i, pat$n)
  )

  # Pattern strictly larger than M: extra positions must come back as 0.
  wider <- adlaplaceFem:::upper_csc_pattern(M + Matrix::Diagonal(n))
  expect_equal(
    adlaplaceFem:::align_gram_to_pattern(M, wider$p, wider$i, wider$n),
    reference_align(M, wider$p, wider$i, wider$n)
  )

  # Pattern with an empty column, so M has nonzeros outside it: those are
  # ignored rather than erroring on an NA subscript.
  drop_col <- 4L
  keep <- (pat$i + 1L) != drop_col
  sub_i <- pat$i[keep]
  counts <- diff(pat$p)
  pcol <- rep.int(seq_along(counts), counts)[keep]
  sub_p <- as.integer(c(0, cumsum(tabulate(pcol, nbins = pat$n))))
  expect_equal(
    adlaplaceFem:::align_gram_to_pattern(M, sub_p, sub_i, pat$n),
    reference_align(M, sub_p, sub_i, pat$n)
  )
})
