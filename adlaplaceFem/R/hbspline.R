#' Hierarchical tensor-product B-spline knot specifications
#'
#' Builds nested B-spline spaces from a list of refinement levels. Each level is
#' a terra `SpatRaster` or a list of rasters defining rectangular refinement
#' regions at the same level.
#'
#' @param knots A list of levels: `list(raster0, raster1, list(raster2a, raster2b), ...)`.
#'   Level 0 defines the coarse padded domain; later levels add knot lines inside
#'   their extents. Legacy `list(x = ..., y = ...)` is not hierarchical.
#' @param degree B-spline degree (must be `>= 2`).
#' @return An `"hb_knots"` object used internally by [fem_bspline()] and [matern()].
#' @export
hb_knots <- function(knots, degree = 2L) {
  degree <- as.integer(degree)[1L]
  if (degree < 2L) {
    stop("degree must be >= 2", call. = FALSE)
  }
  if (inherits(knots, "hb_knots")) {
    if (!identical(knots$degree, degree)) {
      stop("hb_knots degree mismatch", call. = FALSE)
    }
    return(knots)
  }
  levels_raw <- normalize_hb_levels(knots)
  n_levels <- length(levels_raw)
  if (n_levels < 1L) {
    stop("hierarchical knots require at least one level", call. = FALSE)
  }

  level_data <- vector("list", n_levels)
  bx_prev <- NULL
  by_prev <- NULL

  for (lev in seq_len(n_levels)) {
    raw <- levels_raw[[lev]]
    regions <- lapply(raw$rasters, raster_extent_list)
    if (lev == 1L) {
      bp <- axis_breakpoints_from_rasters(raw$rasters)
      bx <- bp$x
      by <- bp$y
      omega <- regions
    } else {
      omega <- regions
      for (reg in regions) {
        if (!rect_contained_in_union(reg, level_data[[lev - 1L]]$omega)) {
          stop(
            "refinement extent at level ", lev - 1L,
            " must lie inside the parent domain",
            call. = FALSE
          )
        }
      }
      regions <- lapply(regions, function(reg) {
        snap_extent_to_breaks(reg, bx_prev, by_prev)
      })
      new_x <- interior_breakpoints_from_rasters(raw$rasters, axis = "x")
      new_y <- interior_breakpoints_from_rasters(raw$rasters, axis = "y")
      bx <- sort(unique(c(bx_prev, new_x)))
      by <- sort(unique(c(by_prev, new_y)))
    }
    open_x <- axis_to_open_knots(bx, degree)
    open_y <- axis_to_open_knots(by, degree)
    level_data[[lev]] <- list(
      breakpoints = list(x = bx, y = by),
      knots = list(x = open_x, y = open_y),
      regions = regions,
      omega = regions,
      n_basis = c(
        x = n_basis_knots(open_x, degree),
        y = n_basis_knots(open_y, degree)
      )
    )
    bx_prev <- bx
    by_prev <- by
  }

  structure(
    list(
      levels = level_data,
      degree = degree,
      n_levels = n_levels
    ),
    class = c("hb_knots", "list")
  )
}

#' @describeIn hb_knots Per-level active basis counts and total hierarchical dof.
#' @export
hb_summary <- function(hb) {
  if (!inherits(hb, "hb_knots")) {
    stop("hb must be an hb_knots object", call. = FALSE)
  }
  active <- hb_active_sets(hb, hb$degree)
  rows <- lapply(seq_along(active), function(lev) {
    data.frame(
      level = lev - 1L,
      n_active = length(active[[lev]]$idx),
      n_x = hb$levels[[lev]]$n_basis["x"],
      n_y = hb$levels[[lev]]$n_basis["y"],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$total_active <- sum(out$n_active)
  out
}

#' Detect hierarchical knot list form
#' @keywords internal
#' @noRd
is_hierarchical_knots <- function(knots) {
  if (!is.list(knots) || inherits(knots, "hb_knots")) {
    return(inherits(knots, "hb_knots"))
  }
  if (!is.null(knots$x) && !is.null(knots$y)) {
    return(FALSE)
  }
  if (length(knots) < 1L) {
    return(FALSE)
  }
  all(vapply(knots, is_hb_level_element, logical(1L)))
}

#' @keywords internal
#' @noRd
is_hb_level_element <- function(x) {
  if (inherits(x, "SpatRaster")) {
    return(TRUE)
  }
  is.list(x) && length(x) >= 1L &&
    all(vapply(x, function(z) inherits(z, "SpatRaster"), logical(1L)))
}

#' @keywords internal
#' @noRd
normalize_hb_levels <- function(knots) {
  if (!is.list(knots)) {
    stop("hierarchical knots must be a list of levels", call. = FALSE)
  }
  lapply(knots, function(level) {
    if (inherits(level, "SpatRaster")) {
      list(rasters = list(level))
    } else if (is.list(level) && all(vapply(level, function(z) {
      inherits(z, "SpatRaster")
    }, logical(1L)))) {
      list(rasters = level)
    } else {
      stop(
        "each hierarchical level must be a SpatRaster or list of SpatRasters",
        call. = FALSE
      )
    }
  })
}

#' @keywords internal
#' @noRd
raster_extent_list <- function(r) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("terra is required for hierarchical SpatRaster knots", call. = FALSE)
  }
  ext <- terra::ext(r)
  list(xmin = ext$xmin, xmax = ext$xmax, ymin = ext$ymin, ymax = ext$ymax)
}

#' @keywords internal
#' @noRd
axis_breakpoints_from_rasters <- function(rasters) {
  xs <- numeric(0)
  ys <- numeric(0)
  for (r in rasters) {
    kl <- knots_from_spatraster(r, degree = 2L)
    xs <- c(xs, kl$x)
    ys <- c(ys, kl$y)
  }
  list(x = sort(unique(xs)), y = sort(unique(ys)))
}

#' @keywords internal
#' @noRd
interior_breakpoints_from_rasters <- function(rasters, axis = c("x", "y")) {
  axis <- match.arg(axis)
  vals <- numeric(0)
  for (r in rasters) {
    ext <- terra::ext(r)
    kl <- knots_from_spatraster(r, degree = 2L)
    lines <- kl[[axis]]
    inside <- if (axis == "x") {
      lines > ext$xmin & lines < ext$xmax
    } else {
      lines > ext$ymin & lines < ext$ymax
    }
    vals <- c(vals, lines[inside])
  }
  sort(unique(vals))
}

#' @keywords internal
#' @noRd
snap_extent_to_breaks <- function(reg, bx, by, tol = 1e-6) {
  snap1 <- function(v, br) {
    br <- sort(unique(as.numeric(br)))
    v <- as.numeric(v)
    c(
      max(br[br <= v[1] + tol]),
      min(br[br >= v[2] - tol])
    )
  }
  xs <- snap1(c(reg$xmin, reg$xmax), bx)
  ys <- snap1(c(reg$ymin, reg$ymax), by)
  list(xmin = xs[1], xmax = xs[2], ymin = ys[1], ymax = ys[2])
}

#' @keywords internal
#' @noRd
rect_contained_in_union <- function(reg, omega) {
  any(vapply(omega, function(o) {
    reg$xmin >= o$xmin - 1e-6 && reg$xmax <= o$xmax + 1e-6 &&
      reg$ymin >= o$ymin - 1e-6 && reg$ymax <= o$ymax + 1e-6
  }, logical(1L)))
}

#' @keywords internal
#' @noRd
box_contained_in_union <- function(xmin, xmax, ymin, ymax, omega, tol = 1e-6) {
  any(vapply(omega, function(o) {
    xmin >= o$xmin - tol && xmax <= o$xmax + tol &&
      ymin >= o$ymin - tol && ymax <= o$ymax + tol
  }, logical(1L)))
}

#' Greville abscissae for an open B-spline knot vector
#' @keywords internal
#' @noRd
greville_abscissae <- function(knots, degree) {
  degree <- as.integer(degree)
  n <- n_basis_knots(knots, degree)
  if (n < 1L) {
    return(numeric(0))
  }
  vapply(seq_len(n), function(i) {
    mean(knots[(i + 1L):(i + degree)], na.rm = TRUE)
  }, numeric(1))
}

#' 1D refinement matrix: coarse basis functions in the next finer basis
#' @keywords internal
#' @noRd
refine_matrix_1d <- function(knots_coarse, knots_fine, degree) {
  degree <- as.integer(degree)
  n_coarse <- n_basis_knots(knots_coarse, degree)
  n_fine <- n_basis_knots(knots_fine, degree)
  if (n_coarse == n_fine && identical(knots_coarse, knots_fine)) {
    return(Matrix::Diagonal(n = n_fine))
  }
  grev <- greville_abscissae(knots_coarse, degree)
  uk <- sort(unique(as.numeric(knots_coarse)))
  # Exact nested-space refinement via collocation on all coarse breakpoint spans
  x_eval <- unique(c(
    grev,
    uk,
    seq(uk[1], uk[length(uk)], length.out = max(4L * n_coarse, 20L))
  ))
  Bfine <- as.matrix(bspline_eval(knots_fine, x_eval, degree, 0L))
  Bcoarse <- as.matrix(bspline_eval(knots_coarse, x_eval, degree, 0L))
  # Bcoarse = Bfine %*% R  =>  R = solve least squares
  R <- Matrix::solve(Matrix::crossprod(Bfine), Matrix::crossprod(Bfine, Bcoarse))
  methods::as(R, "dgCMatrix")
}

#' Support interval [left, right] for 1D B-spline basis index
#' @keywords internal
#' @noRd
bspline_support_1d <- function(knots, degree, index) {
  degree <- as.integer(degree)
  i <- as.integer(index)
  c(left = knots[i + 1L], right = knots[i + degree + 1L])
}

#' Active hierarchical basis indices per level
#' @keywords internal
#' @noRd
hb_active_sets <- function(hb, degree) {
  degree <- as.integer(degree)
  n_levels <- hb$n_levels
  omega_next <- list(list(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf))
  # placeholder; level L uses empty omega_{L+1}
  active <- vector("list", n_levels)

  for (lev in seq_len(n_levels)) {
    kn <- hb$levels[[lev]]$knots
    nx <- hb$levels[[lev]]$n_basis["x"]
    ny <- hb$levels[[lev]]$n_basis["y"]
    omega_l <- hb$levels[[lev]]$omega
    omega_lp1 <- if (lev < n_levels) {
      hb$levels[[lev + 1L]]$omega
    } else {
      list()
    }

    idx <- integer(0)
    pairs <- vector("list", 0)
    for (j in seq_len(ny)) {
      sy <- bspline_support_1d(kn$y, degree, j)
      for (i in seq_len(nx)) {
        sx <- bspline_support_1d(kn$x, degree, i)
        in_l <- box_contained_in_union(
          sx["left"], sx["right"], sy["left"], sy["right"], omega_l
        )
        in_lp1 <- length(omega_lp1) > 0L && box_contained_in_union(
          sx["left"], sx["right"], sy["left"], sy["right"], omega_lp1
        )
        if (lev == 1L) {
          # level 0 uses full domain; active when not contained in level 1 region
          if (!in_lp1 || length(omega_lp1) == 0L) {
            idx <- c(idx, i + (j - 1L) * nx)
            pairs[[length(pairs) + 1L]] <- c(i, j)
          }
        } else if (in_l && !in_lp1) {
          idx <- c(idx, i + (j - 1L) * nx)
          pairs[[length(pairs) + 1L]] <- c(i, j)
        }
      }
    }
    active[[lev]] <- list(idx = idx, ij = pairs, n_x = nx, n_y = ny)
  }
  active
}

#' Cumulative 1D refinement from level l to finest
#' @keywords internal
#' @noRd
hb_cumulative_refine_1d <- function(hb, axis = c("x", "y")) {
  axis <- match.arg(axis)
  n_levels <- hb$n_levels
  degree <- hb$degree
  mats <- vector("list", n_levels - 1L)
  for (lev in seq_len(n_levels - 1L)) {
    coarse <- hb$levels[[lev]]$knots[[axis]]
    fine <- hb$levels[[lev + 1L]]$knots[[axis]]
    mats[[lev]] <- refine_matrix_1d(coarse, fine, degree)
  }
  cum <- vector("list", n_levels)
  cum[[n_levels]] <- Matrix::Diagonal(
    n = hb$levels[[n_levels]]$n_basis[[axis]]
  )
  if (n_levels > 1L) {
    for (lev in seq(n_levels - 1L, 1L)) {
      cum[[lev]] <- cum[[lev + 1L]] %*% mats[[lev]]
    }
  }
  cum
}

#' Sparse basis map S: hierarchical coefficients to finest-level basis
#' @param truncate If `TRUE` (default), build truncated (THB) basis functions for
#'   a partition of unity; if `FALSE`, use the untruncated hierarchical embedding.
#' @keywords internal
#' @noRd
hb_basis_map <- function(hb, degree, truncate = TRUE) {
  degree <- as.integer(degree)
  active <- hb_active_sets(hb, degree)
  n_levels <- hb$n_levels
  nx_f <- hb$levels[[n_levels]]$n_basis["x"]
  ny_f <- hb$levels[[n_levels]]$n_basis["y"]

  cum_x <- hb_cumulative_refine_1d(hb, "x")
  cum_y <- hb_cumulative_refine_1d(hb, "y")
  refine_x <- hb_pairwise_refine_1d(hb, "x")
  refine_y <- hb_pairwise_refine_1d(hb, "y")

  cols <- list()
  for (lev in seq_len(n_levels)) {
    if (!length(active[[lev]]$idx)) {
      next
    }
    nx <- active[[lev]]$n_x
    ny <- active[[lev]]$n_y
    for (k in seq_along(active[[lev]]$idx)) {
      ij <- active[[lev]]$ij[[k]]
      i <- ij[1L]
      j <- ij[2L]
      v <- sparse_unit(nx * ny, i + (j - 1L) * nx)
      if (lev < n_levels && truncate) {
        Tx <- refine_x[[lev]]
        Ty <- refine_y[[lev]]
        T <- Matrix::kronecker(Ty, Tx)
        w <- T %*% v
        mask <- active_mask_vector(active[[lev + 1L]]$idx, nrow(T))
        w_sub <- w * mask
        Rx <- cum_x[[lev + 1L]]
        Ry <- cum_y[[lev + 1L]]
        T2 <- Matrix::kronecker(Ry, Rx)
        col <- T2 %*% (w - w_sub)
      } else if (lev < n_levels) {
        Rx <- cum_x[[lev]]
        Ry <- cum_y[[lev]]
        col <- Matrix::kronecker(Rx[, i, drop = FALSE], Ry[, j, drop = FALSE])
      } else {
        col <- sparse_unit(nx_f * ny_f, i + (j - 1L) * nx_f)
      }
      cols[[length(cols) + 1L]] <- col
    }
  }
  if (!length(cols)) {
    stop("hierarchical basis has no active functions", call. = FALSE)
  }
  S <- Reduce(Matrix::cbind2, cols)
  methods::as(S, "dgCMatrix")
}

#' Pairwise 1D refinement matrices between consecutive levels
#' @keywords internal
#' @noRd
hb_pairwise_refine_1d <- function(hb, axis = c("x", "y")) {
  axis <- match.arg(axis)
  n_levels <- hb$n_levels
  degree <- hb$degree
  mats <- vector("list", n_levels - 1L)
  for (lev in seq_len(n_levels - 1L)) {
    coarse <- hb$levels[[lev]]$knots[[axis]]
    fine <- hb$levels[[lev + 1L]]$knots[[axis]]
    mats[[lev]] <- refine_matrix_1d(coarse, fine, degree)
  }
  mats
}

#' Sparse unit vector
#' @keywords internal
#' @noRd
sparse_unit <- function(n, index) {
  methods::as(
    Matrix::sparseMatrix(
      i = as.integer(index),
      j = 1L,
      x = 1,
      dims = c(as.integer(n), 1L)
    ),
    "dgCMatrix"
  )
}

#' Mask vector with 1 at active tensor indices
#' @keywords internal
#' @noRd
active_mask_vector <- function(active_idx, n) {
  m <- numeric(n)
  if (length(active_idx)) {
    m[active_idx] <- 1
  }
  m
}

#' Map a sparse Gram from finest level to hierarchical basis
#' @keywords internal
#' @noRd
hb_project_gram <- function(M, S) {
  methods::as(Matrix::t(S) %*% M %*% S, "dgCMatrix")
}
