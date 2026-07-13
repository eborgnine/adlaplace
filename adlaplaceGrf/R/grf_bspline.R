#' Tensor-product B-spline FEM ingredients for a 2D Matérn GRF
#'
#' Builds the design matrix of tensor-product B-splines at evaluation
#' coordinates together with sparse mass, stiffness, and higher-order Gram
#' matrices used in the SPDE FEM precision for operator orders α = 2 and 3.
#'
#' Coefficient ordering is column-major (`vec`) over the `n_x` by `n_y` array of
#' tensor basis indices: basis `(i, j)` maps to column `i + (j - 1) * n_x`.
#'
#' @param coords Evaluation locations. Methods are defined for `list`
#'   (`list(x, y)`), and (if **terra** is installed) `SpatRaster` / `SpatVector`.
#' @param degree B-spline degree; must be `>= 2`. Default `2` (α = 2). Use
#'   `degree = 3` when assembling `G3` for α = 3.
#' @param knots Knot specification: `NULL` (regular open knots from coordinate
#'   ranges), `list(x, y)` of knot vectors, or (terra) a `SpatRaster` whose
#'   grid defines knot lines.
#' @param n_interior Number of interior knots per axis when `knots = NULL`
#'   (scalar or length-2).
#' @param ... Passed to methods (unused for `list`).
#'
#' @return A list with:
#' \describe{
#'   \item{A}{`dgCMatrix` design, `n` by `n_x * n_y`}
#'   \item{C, G, G2}{sparse Grams}
#'   \item{G3}{sparse Gram or `NULL` if `degree < 3`}
#'   \item{degree, knots, n_basis}{metadata (`n_basis = c(n_x, n_y)`)}
#' }
#'
#' @export
setGeneric(
  "grf_bspline",
  function(coords, degree = 2L, knots = NULL, ...) {
    standardGeneric("grf_bspline")
  }
)

#' @rdname grf_bspline
#' @export
setMethod(
  "grf_bspline",
  signature(coords = "list"),
  function(coords, degree = 2L, knots = NULL, n_interior = 5L, ...) {
    if (is.null(coords$x) || is.null(coords$y)) {
      stop("coords list must have numeric components x and y")
    }
    x <- as.numeric(coords$x)
    y <- as.numeric(coords$y)
    if (length(x) != length(y)) {
      stop("coords$x and coords$y must have the same length")
    }
    if (any(!is.finite(x)) || any(!is.finite(y))) {
      stop("coords must be finite")
    }
    grf_bspline_xy(x, y, degree = degree, knots = knots, n_interior = n_interior)
  }
)

#' Core assembly from numeric x, y and list knots
#' @keywords internal
grf_bspline_xy <- function(x, y, degree = 2L, knots = NULL, n_interior = 5L) {
  degree <- as.integer(degree)
  if (degree < 2L) {
    stop("degree must be >= 2")
  }
  kn <- resolve_knots_list(
    knots,
    xlim = range(x),
    ylim = range(y),
    degree = degree,
    n_interior = n_interior
  )
  nx <- n_basis_knots(kn$x, degree)
  ny <- n_basis_knots(kn$y, degree)

  Cx <- gram_1d(kn$x, degree, 0L, 0L)
  Cy <- gram_1d(kn$y, degree, 0L, 0L)
  Gx <- gram_1d(kn$x, degree, 1L, 1L)
  Gy <- gram_1d(kn$y, degree, 1L, 1L)
  G2x <- gram_1d(kn$x, degree, 2L, 2L)
  G2y <- gram_1d(kn$y, degree, 2L, 2L)

  C <- methods::as(Matrix::kronecker(Cy, Cx), "dgCMatrix")
  G <- methods::as(
    Matrix::drop0(Matrix::kronecker(Cy, Gx) + Matrix::kronecker(Gy, Cx)),
    "dgCMatrix"
  )
  G2 <- methods::as(
    Matrix::drop0(
      Matrix::kronecker(Cy, G2x) +
        Matrix::kronecker(G2y, Cx) +
        2 * Matrix::kronecker(Gy, Gx)
    ),
    "dgCMatrix"
  )

  G3 <- NULL
  if (degree >= 3L) {
    Mx00 <- Cx
    My00 <- Cy
    Mx11 <- Gx
    My11 <- Gy
    Mx22 <- G2x
    My22 <- G2y
    Mx33 <- gram_1d(kn$x, degree, 3L, 3L)
    My33 <- gram_1d(kn$y, degree, 3L, 3L)
    Mx31 <- gram_1d(kn$x, degree, 3L, 1L)
    Mx13 <- gram_1d(kn$x, degree, 1L, 3L)
    Mx20 <- gram_1d(kn$x, degree, 2L, 0L)
    Mx02 <- gram_1d(kn$x, degree, 0L, 2L)
    My31 <- gram_1d(kn$y, degree, 3L, 1L)
    My13 <- gram_1d(kn$y, degree, 1L, 3L)
    My20 <- gram_1d(kn$y, degree, 2L, 0L)
    My02 <- gram_1d(kn$y, degree, 0L, 2L)

    G3 <- Matrix::kronecker(My00, Mx33) +
      Matrix::kronecker(My02, Mx31) +
      Matrix::kronecker(My20, Mx13) +
      Matrix::kronecker(My22, Mx11) +
      Matrix::kronecker(My11, Mx22) +
      Matrix::kronecker(My13, Mx20) +
      Matrix::kronecker(My31, Mx02) +
      Matrix::kronecker(My33, Mx00)
    G3 <- methods::as(Matrix::drop0(G3), "dgCMatrix")
  }

  A <- tensor_design(x, y, kn$x, kn$y, degree)

  list(
    A = A,
    C = C,
    G = G,
    G2 = G2,
    G3 = G3,
    degree = degree,
    knots = kn,
    n_basis = c(x = nx, y = ny)
  )
}

#' Sparse design matrix for tensor-product B-splines at (x, y)
#' @keywords internal
tensor_design <- function(x, y, knots_x, knots_y, degree) {
  Bx <- bspline_eval(knots_x, x, degree, 0L)
  By <- bspline_eval(knots_y, y, degree, 0L)
  nx <- ncol(Bx)
  n <- length(x)
  Bx_t <- Matrix::summary(methods::as(Bx, "TsparseMatrix"))
  By_t <- Matrix::summary(methods::as(By, "TsparseMatrix"))

  max_nnz <- as.integer(n * (degree + 1L)^2)
  rows <- integer(max_nnz)
  cols <- integer(max_nnz)
  vals <- numeric(max_nnz)
  pos <- 0L

  split_x <- split(seq_len(nrow(Bx_t)), Bx_t$i)
  split_y <- split(seq_len(nrow(By_t)), By_t$i)
  common_rows <- intersect(names(split_x), names(split_y))
  for (rn in common_rows) {
    ix <- split_x[[rn]]
    iy <- split_y[[rn]]
    r <- as.integer(rn)
    for (a in ix) {
      i <- Bx_t$j[a]
      bx <- Bx_t$x[a]
      for (b in iy) {
        pos <- pos + 1L
        if (pos > length(rows)) {
          rows <- c(rows, integer(length(rows)))
          cols <- c(cols, integer(length(cols)))
          vals <- c(vals, numeric(length(vals)))
        }
        rows[pos] <- r
        cols[pos] <- i + (By_t$j[b] - 1L) * nx
        vals[pos] <- bx * By_t$x[b]
      }
    }
  }
  if (pos == 0L) {
    return(methods::as(Matrix::Matrix(0, n, nx * ncol(By), sparse = TRUE), "dgCMatrix"))
  }
  methods::as(
    Matrix::sparseMatrix(
      i = rows[seq_len(pos)],
      j = cols[seq_len(pos)],
      x = vals[seq_len(pos)],
      dims = c(n, nx * ncol(By))
    ),
    "dgCMatrix"
  )
}
