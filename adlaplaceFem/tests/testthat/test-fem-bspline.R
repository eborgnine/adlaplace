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

sparse_max_abs_diff <- function(a, b) {
  d <- methods::as(a - b, "dgCMatrix")@x
  if (!length(d)) 0 else max(abs(d))
}

# Row r of the tensor design is kronecker(By[r, ], Bx[r, ]), with x-index fastest.
tensor_kronecker_ref <- function(x, y, knots_x, knots_y, degree) {
  Bx <- as.matrix(adlaplaceFem:::bspline_eval(knots_x, x, degree, 0L))
  By <- as.matrix(adlaplaceFem:::bspline_eval(knots_y, y, degree, 0L))
  nx <- ncol(Bx)
  ref <- matrix(0, length(x), nx * ncol(By))
  for (r in seq_along(x)) {
    ref[r, ] <- as.numeric(kronecker(By[r, ], Bx[r, ]))
  }
  methods::as(Matrix::Matrix(ref, sparse = TRUE), "dgCMatrix")
}

test_that("tensor_design matches a per-cell Kronecker product", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 4))
  x <- c(0.1, 0.5, 1, -0.2, 0.8, 0.3)
  y <- c(0.2, 0.2, 0.9, 0.4, 1.2, 0)
  for (degree in c(2L, 3L)) {
    fem <- fem_bspline(data.frame(x = x, y = y), knots_list, degree = degree)
    A <- adlaplaceFem:::tensor_design(x, y, fem$knots$x, fem$knots$y, degree)
    ref <- tensor_kronecker_ref(x, y, fem$knots$x, fem$knots$y, degree)
    expect_equal(dim(A), dim(ref))
    expect_equal(sparse_max_abs_diff(A, ref), 0)
  }
})

test_that("serial design blocks match one-shot evaluation, including row order", {
  knots_list <- list(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5))
  x <- seq(0.05, 0.95, by = 0.1)
  y <- seq(0.95, 0.05, by = -0.1)
  fem <- fem_bspline(data.frame(x = x, y = y), knots_list, degree = 2L)
  A1 <- adlaplaceFem:::fem_design_xy(fem, x, y)
  A3 <- adlaplaceFem:::fem_design_blocks(fem, x, y, n_blocks = 3L)
  expect_equal(dim(A3), dim(A1))
  expect_equal(sparse_max_abs_diff(A3, A1), 0)
  # First block is only the leading coordinates; a reversed rbind would miss them.
  n1 <- length(parallel::splitIndices(length(x), 3L)[[1L]])
  expect_equal(
    sparse_max_abs_diff(A3[seq_len(n1), ], A1[seq_len(n1), ]),
    0
  )
})

test_that("hierarchical design blocks match a per-cell fine basis times S", {
  kn <- hb_knots(
    list(xmin = 0, xmax = 1, ymin = 0, ymax = 1, resolution = 0.25),
    c(0.25, 0.75, 0.25, 0.75),
    fact = 2
  )
  pts <- expand.grid(
    x = seq(0.05, 0.95, by = 0.15),
    y = seq(0.1, 0.9, by = 0.2)
  )
  x <- pts$x
  y <- pts$y
  fem <- fem_bspline(pts, kn, degree = 2L)
  expect_true(!is.null(fem$S))
  kf <- fem$knots_finest
  A_fine <- tensor_kronecker_ref(x, y, kf$x, kf$y, fem$degree)
  A_ref <- methods::as(A_fine %*% fem$S, "dgCMatrix")
  A <- adlaplaceFem:::fem_design_xy(fem, x, y)
  A_blocks <- adlaplaceFem:::fem_design_blocks(fem, x, y, n_blocks = 4L)
  expect_equal(sparse_max_abs_diff(A, A_ref), 0)
  expect_equal(sparse_max_abs_diff(A_blocks, A), 0)
})
