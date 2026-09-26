test_that("fem_chol_pattern is accepted by adlaplace::takahashi_davis", {
  n <- 6L
  Q <- Matrix::bandSparse(
    n,
    k = c(-1L, 0L, 1L),
    diagonals = list(rep(-0.4, n - 1L), rep(2, n), rep(-0.4, n - 1L))
  )
  Q <- methods::as(Q, "dgCMatrix")
  S <- adlaplace::takahashi_davis(Q)
  S2 <- adlaplace::takahashi_davis(Q, chol = fem_chol_pattern(Q))
  expect_equal(as.matrix(S2), as.matrix(S), tolerance = 1e-12)
})
