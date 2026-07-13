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

#' Normalize knots to list(x, y) open knot vectors
#'
#' @param knots `NULL`, `list(x, y)`, or (via callers) already-resolved lists.
#' @param xlim,ylim Domain ranges used when `knots` is `NULL`.
#' @param degree B-spline degree.
#' @param n_interior Default number of interior knots per axis when `knots` is
#'   `NULL` (scalar or length-2).
#' @return `list(x = ..., y = ...)` of open knot vectors.
#' @keywords internal
resolve_knots_list <- function(knots, xlim, ylim, degree, n_interior = 5L) {
  if (is.null(knots)) {
    ni <- rep_len(as.integer(n_interior), 2L)
    return(list(
      x = open_knot_vector(xlim, degree, n_interior = ni[1]),
      y = open_knot_vector(ylim, degree, n_interior = ni[2])
    ))
  }
  if (!is.list(knots) || is.null(knots$x) || is.null(knots$y)) {
    stop("knots must be NULL or list(x = ..., y = ...)")
  }
  kx <- as.numeric(knots$x)
  ky <- as.numeric(knots$y)
  # If user passed interior-only or unique knots, expand to open vectors.
  if (length(unique(kx)) == length(kx)) {
    kx <- open_knot_vector(range(kx), degree, interior = if (length(kx) > 2L) kx[-c(1L, length(kx))] else NULL)
  }
  if (length(unique(ky)) == length(ky)) {
    ky <- open_knot_vector(range(ky), degree, interior = if (length(ky) > 2L) ky[-c(1L, length(ky))] else NULL)
  }
  list(x = kx, y = ky)
}

#' Number of B-spline basis functions for an open knot vector
#' @keywords internal
n_basis_knots <- function(knots, degree) {
  as.integer(length(knots) - (degree + 1L))
}
