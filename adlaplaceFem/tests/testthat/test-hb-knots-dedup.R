test_that("hb_knots collapses near-duplicate sibling breakpoints", {
  # Sibling level-2 boxes with overlapping y ranges but different ymin so
  # cell-centre arithmetic can produce the same nominal lines with ULP noise
  # when unique() is exact-equality only. Level-1 is a partial refinement
  # (not the full outer) so the hierarchical active set stays well-formed.
  outer <- list(
    xmin = 0, xmax = 32, ymin = 0, ymax = 24, resolution = c(4, 4)
  )
  kn <- hb_knots(
    outer = outer,
    inner = list(
      c(4, 28, 4, 20),
      list(
        c(4, 16, 6, 18),
        c(12, 28, 4, 20)
      )
    ),
    fact = 2L
  )
  for (lev in seq_along(kn)) {
    res <- min(kn[[lev]]$rasters[[1L]]$resolution)
    for (axis in c("x", "y")) {
      spans <- diff(kn[[lev]]$knots[[axis]])
      expect_gt(min(spans), 0.4 * res)
    }
  }
  hb <- hb_basis(kn, degree = 3L)
  fem <- adlaplaceFem:::fem_bspline_hb(numeric(0), numeric(0), hb, degree = 3L)
  pay <- fem_precision_payload(fem, alpha = 3L)
  dbg <- fem_logdet_debug(pay, range = 8, sd = 0.5)
  expect_true(is.finite(dbg$logdet))
  expect_equal(dbg$n_bad, 0L)
})

test_that("hb_basis errors on hand-built near-duplicate knot lines", {
  outer <- list(
    xmin = 0, xmax = 8, ymin = 0, ymax = 8, resolution = c(2, 2)
  )
  kn <- hb_knots(outer = outer, inner = NULL)
  # Inject a near-duplicate y line (1e-10 offset) that survives unique().
  y <- kn[[1L]]$knots$y
  kn[[1L]]$knots$y <- sort(c(y, y[2L] + 1e-10))
  expect_error(
    hb_basis(kn, degree = 2L),
    "degenerate knot span"
  )
})
