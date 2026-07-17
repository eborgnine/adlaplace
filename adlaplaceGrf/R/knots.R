#' Resolve axis knot lines from a terra SpatRaster
#'
#' Canonical SpatRaster-to-knots rule used by [matern()] and
#' [grf_bspline()] terra methods. Returns unique increasing knot-line
#' positions on each axis for [axis_to_open_knots()]: raster extent
#' endpoints (`xmin`, `xmax`, `ymin`, `ymax`) are always included;
#' cell-center coordinates strictly inside the extent become interior knot
#' lines.
#'
#' @param r A terra `SpatRaster` defining the knot grid.
#' @param degree B-spline degree (unused; reserved for future grid rules).
#' @return `list(x = ..., y = ...)` of numeric knot-line positions.
#' @keywords internal
#' @noRd
knots_from_spatraster <- function(r, degree) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "terra is required for SpatRaster knots; install.packages(\"terra\")",
      call. = FALSE
    )
  }
  ext <- terra::ext(r)
  xs <- terra::xFromCol(r, seq_len(terra::ncol(r)))
  ys <- terra::yFromRow(r, seq_len(terra::nrow(r)))
  list(
    x = c(ext$xmin, xs[xs > ext$xmin & xs < ext$xmax], ext$xmax),
    y = c(ext$ymin, ys[ys > ext$ymin & ys < ext$ymax], ext$ymax)
  )
}

#' Build an open B-spline knot vector on an interval
#'
#' Boundary knots are repeated `degree + 1` times (open knot vector).
#' Interior knots are equally spaced unless supplied.
#'
#' @param range Numeric length-2 vector `c(min, max)`.
#' @param degree B-spline degree (order = degree + 1).
#' @param n_interior Number of interior knots (ignored if `interior` given).
#' @param interior Optional numeric vector of interior knots strictly inside
#'   `range`.
#' @return Numeric nondecreasing knot vector.
#' @keywords internal
#' @noRd
open_knot_vector <- function(range, degree, n_interior = 5L, interior = NULL) {
  range <- range(as.numeric(range), na.rm = TRUE)
  if (length(range) != 2L || !is.finite(range[1]) || !is.finite(range[2])) {
    stop("range must be two finite numbers")
  }
  if (range[2] <= range[1]) {
    stop("range must satisfy max > min")
  }
  degree <- as.integer(degree)
  if (degree < 2L) {
    stop("degree must be >= 2 for higher-order FEM Grams")
  }
  if (is.null(interior)) {
    n_interior <- as.integer(n_interior)
    if (n_interior < 0L) {
      stop("n_interior must be non-negative")
    }
    if (n_interior == 0L) {
      interior <- numeric(0)
    } else {
      interior <- seq(range[1], range[2], length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
    }
  } else {
    interior <- sort(as.numeric(interior))
    if (any(interior <= range[1] | interior >= range[2])) {
      stop("interior knots must lie strictly inside range")
    }
  }
  c(rep(range[1], degree + 1L), interior, rep(range[2], degree + 1L))
}

#' Expand unique axis breakpoints to an open knot vector
#'
#' @param t Strictly increasing (or unique) knot-line positions including the
#'   domain endpoints, e.g. `seq(-0.2, 1.2, by = 0.1)`.
#' @param degree B-spline degree.
#' @return Open knot vector with endpoints repeated `degree + 1` times.
#' @keywords internal
#' @noRd
axis_to_open_knots <- function(t, degree) {
  t <- sort(unique(as.numeric(t)))
  if (length(t) < 2L || any(!is.finite(t))) {
    stop("each knots axis must have at least two finite positions")
  }
  interior <- if (length(t) > 2L) t[-c(1L, length(t))] else numeric(0)
  open_knot_vector(range(t), degree, interior = interior)
}

#' Normalize knots to list(x, y) open knot vectors
#'
#' User form is a `knots_list`-style object:
#' `list(x = seq(...), y = seq(...))` of unique, increasing knot-line
#' positions on each axis (domain endpoints included). Those are expanded to
#' open B-spline knot vectors. Already-open vectors (with repeated endpoints)
#' are left as-is.
#'
#' @param knots `list(x, y)` of axis knot positions / open vectors.
#' @param degree B-spline degree.
#' @return `list(x = ..., y = ...)` of open knot vectors.
#' @keywords internal
#' @noRd
resolve_knots_list <- function(knots, degree) {
  if (!is.list(knots) || is.null(knots$x) || is.null(knots$y)) {
    stop("knots must be list(x = ..., y = ...), e.g. list(x = seq(...), y = seq(...))")
  }
  kx <- as.numeric(knots$x)
  ky <- as.numeric(knots$y)
  # Unique axis positions (knots_list style) -> open knot vectors.
  # Repeated endpoints mean the user already passed an open vector.
  if (length(unique(kx)) == length(kx)) {
    kx <- axis_to_open_knots(kx, degree)
  }
  if (length(unique(ky)) == length(ky)) {
    ky <- axis_to_open_knots(ky, degree)
  }
  list(x = kx, y = ky)
}

#' Number of B-spline basis functions for an open knot vector
#' @keywords internal
#' @noRd
n_basis_knots <- function(knots, degree) {
  as.integer(length(knots) - (degree + 1L))
}
