#' Hierarchical tensor-product B-spline knot specifications
#'
#' Builds nested knot-line locations from a coarse terra raster and optional
#' inner refinement extents. Knot **lines** do not depend on B-spline degree;
#' pass the result to [hb_basis()] or [matern()] to form the THB space.
#'
#' Each element of `inner` is one refinement level: a single extent, or a list
#' of extents for sibling boxes at that level. `list(ext)` and `list(list(ext))`
#' are the same one-level, one-box specification. `list(ext1, ext2)` is two
#' levels; siblings are `list(list(ext1, ext2))`. Extents are snapped to `outer`
#' cell edges and clipped to the parent level.
#'
#' @param outer A terra `SpatRaster` giving the padded domain (level 0) and
#'   level-0 cell size.
#' @param inner `NULL` (outer only), one extent, or a list of levels. An extent
#'   is a `SpatExtent`, `c(xmin, xmax, ymin, ymax)`, or a `SpatRaster` /
#'   `SpatVector` (extent only; resolution is ignored).
#' @param fact Integer `>= 2`, recycled to `length(inner)`. `fact[i]` is the
#'   cell-size ratio from the previous level, so level-`i` resolution is
#'   `res(outer) / prod(fact[1:i])`.
#' @return An `"hb_knots"` object (breakpoints, snapped regions, resolution).
#' @export
hb_knots <- function(outer, inner = NULL, fact = 2L) {
  if (inherits(outer, "hb_knots")) {
    return(outer)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "terra is required for hb_knots(); install.packages(\"terra\")",
      call. = FALSE
    )
  }
  if (!inherits(outer, "SpatRaster")) {
    stop("outer must be a terra SpatRaster", call. = FALSE)
  }

  inner_levels <- normalize_inner_levels(inner)
  n_refine <- length(inner_levels)
  fact <- as.integer(fact)
  if (n_refine == 0L) {
    fact <- integer(0)
  } else {
    if (length(fact) == 1L) {
      fact <- rep(fact, n_refine)
    }
    if (length(fact) != n_refine) {
      stop("fact must have length 1 or length(inner)", call. = FALSE)
    }
    if (any(is.na(fact)) || any(fact < 2L)) {
      stop("each fact value must be an integer >= 2", call. = FALSE)
    }
  }

  n_levels <- n_refine + 1L
  edges <- outer_cell_edges(outer)
  res0 <- c(terra::xres(outer), terra::yres(outer))
  bp0 <- knots_from_spatraster(outer, degree = 2L)
  omega0 <- list(raster_extent_list(outer))

  level_data <- vector("list", n_levels)
  level_data[[1L]] <- list(
    breakpoints = list(x = bp0$x, y = bp0$y),
    regions = omega0,
    omega = omega0,
    resolution = res0
  )

  bx_prev <- bp0$x
  by_prev <- bp0$y
  res_prev <- res0

  for (i in seq_len(n_refine)) {
    regions <- lapply(inner_levels[[i]], function(reg) {
      snapped <- snap_extent_to_breaks(reg, edges$x, edges$y)
      clipped <- clip_region_to_omega(snapped, level_data[[i]]$omega)
      if (is.null(clipped)) {
        stop(
          "refinement extent at inner level ", i,
          " must lie inside the parent domain",
          call. = FALSE
        )
      }
      clipped
    })
    res_i <- res_prev / as.numeric(fact[i])
    rasters <- lapply(regions, function(reg) region_raster(reg, res_i))
    new_x <- interior_breakpoints_from_rasters(rasters, axis = "x")
    new_y <- interior_breakpoints_from_rasters(rasters, axis = "y")
    bx <- sort(unique(c(bx_prev, new_x)))
    by <- sort(unique(c(by_prev, new_y)))
    level_data[[i + 1L]] <- list(
      breakpoints = list(x = bx, y = by),
      regions = regions,
      omega = regions,
      resolution = res_i
    )
    bx_prev <- bx
    by_prev <- by
    res_prev <- res_i
  }

  structure(
    list(
      levels = level_data,
      n_levels = n_levels,
      fact = unname(fact)
    ),
    class = c("hb_knots", "list")
  )
}

#' Truncated hierarchical B-spline space
#'
#' Expands [hb_knots()] breakpoint lines to open knot vectors and builds the
#' THB embedding (active sets and sparse map `S`) for a B-spline `degree`.
#'
#' @param knots An `"hb_knots"` object from [hb_knots()].
#' @param degree B-spline degree (must be `>= 2`).
#' @return An `"hb_basis"` object used by [fem_bspline()] and [matern()].
#' @export
hb_basis <- function(knots, degree = 2L) {
  degree <- as.integer(degree)[1L]
  if (degree < 2L) {
    stop("degree must be >= 2", call. = FALSE)
  }
  if (inherits(knots, "hb_basis")) {
    if (!identical(knots$degree, degree)) {
      stop("hb_basis degree mismatch", call. = FALSE)
    }
    return(knots)
  }
  if (!inherits(knots, "hb_knots")) {
    stop("knots must be an hb_knots object from hb_knots()", call. = FALSE)
  }

  levels <- lapply(knots$levels, function(lev) {
    open_x <- axis_to_open_knots(lev$breakpoints$x, degree)
    open_y <- axis_to_open_knots(lev$breakpoints$y, degree)
    lev$knots <- list(x = open_x, y = open_y)
    lev$n_basis <- c(
      x = n_basis_knots(open_x, degree),
      y = n_basis_knots(open_y, degree)
    )
    lev
  })

  hb <- structure(
    list(
      levels = levels,
      degree = degree,
      n_levels = knots$n_levels,
      fact = knots$fact,
      knots = knots
    ),
    class = c("hb_basis", "list")
  )
  hb$active <- hb_active_sets(hb, degree)
  hb$S <- hb_basis_map(hb, degree)
  hb
}

#' @describeIn hb_basis Per-level active basis counts and total hierarchical dof.
#' @param hb An `"hb_basis"` object from [hb_basis()].
#' @export
hb_summary <- function(hb) {
  if (!inherits(hb, "hb_basis")) {
    stop("hb must be an hb_basis object from hb_basis()", call. = FALSE)
  }
  active <- hb$active
  if (is.null(active)) {
    active <- hb_active_sets(hb, hb$degree)
  }
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

#' @keywords internal
#' @noRd
is_extent_like <- function(x) {
  inherits(x, "SpatExtent") ||
    inherits(x, "SpatRaster") ||
    inherits(x, "SpatVector") ||
    (is.numeric(x) && length(x) == 4L && is.null(dim(x)))
}

#' @keywords internal
#' @noRd
coerce_extent_region <- function(x) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("terra is required for hb_knots()", call. = FALSE)
  }
  e <- if (inherits(x, "SpatExtent")) {
    x
  } else if (inherits(x, "SpatRaster") || inherits(x, "SpatVector")) {
    terra::ext(x)
  } else if (is.numeric(x) && length(x) == 4L) {
    terra::ext(x[1], x[2], x[3], x[4])
  } else {
    stop(
      "each inner region must be a SpatExtent, length-4 numeric, ",
      "SpatRaster, or SpatVector",
      call. = FALSE
    )
  }
  list(xmin = e$xmin, xmax = e$xmax, ymin = e$ymin, ymax = e$ymax)
}

#' @keywords internal
#' @noRd
normalize_inner_levels <- function(inner) {
  if (is.null(inner)) {
    return(list())
  }
  if (is_extent_like(inner)) {
    return(list(list(coerce_extent_region(inner))))
  }
  if (!is.list(inner) || !length(inner)) {
    stop("inner must be NULL, an extent, or a list of levels", call. = FALSE)
  }
  lapply(inner, function(level) {
    if (is_extent_like(level)) {
      list(coerce_extent_region(level))
    } else if (is.list(level) && length(level) >= 1L) {
      lapply(level, coerce_extent_region)
    } else {
      stop(
        "each inner level must be an extent or a list of extents",
        call. = FALSE
      )
    }
  })
}

#' @keywords internal
#' @noRd
outer_cell_edges <- function(r) {
  e <- terra::ext(r)
  list(
    x = as.numeric(e$xmin) + seq(0, terra::ncol(r)) * terra::xres(r),
    y = as.numeric(e$ymin) + seq(0, terra::nrow(r)) * terra::yres(r)
  )
}

#' @keywords internal
#' @noRd
region_raster <- function(reg, resolution) {
  dx <- resolution[1]
  dy <- if (length(resolution) > 1L) resolution[2] else resolution[1]
  nc <- max(1L, as.integer(round((reg$xmax - reg$xmin) / dx)))
  nr <- max(1L, as.integer(round((reg$ymax - reg$ymin) / dy)))
  terra::rast(
    terra::ext(reg$xmin, reg$xmax, reg$ymin, reg$ymax),
    ncols = nc,
    nrows = nr
  )
}

#' Detect the retired list-of-SpatRaster hierarchical form
#' @keywords internal
#' @noRd
is_legacy_raster_levels <- function(knots) {
  if (!is.list(knots) || inherits(knots, "hb_knots") ||
      inherits(knots, "hb_basis")) {
    return(FALSE)
  }
  if (!is.null(knots$x) && !is.null(knots$y)) {
    return(FALSE)
  }
  if (length(knots) < 1L) {
    return(FALSE)
  }
  all(vapply(knots, function(x) {
    inherits(x, "SpatRaster") ||
      (is.list(x) && length(x) >= 1L &&
         all(vapply(x, function(z) inherits(z, "SpatRaster"), logical(1L))))
  }, logical(1L)))
}

#' @keywords internal
#' @noRd
stop_legacy_raster_levels <- function() {
  stop(
    "hierarchical knots must be created with hb_knots(outer, inner, fact)",
    call. = FALSE
  )
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
    lo <- br[br <= v[1] + tol]
    hi <- br[br >= v[2] - tol]
    if (!length(lo) || !length(hi)) {
      stop(
        "inner extent cannot be snapped to outer cell boundaries",
        call. = FALSE
      )
    }
    c(max(lo), min(hi))
  }
  xs <- snap1(c(reg$xmin, reg$xmax), bx)
  ys <- snap1(c(reg$ymin, reg$ymax), by)
  if (xs[2] <= xs[1] || ys[2] <= ys[1]) {
    stop("snapped inner extent has zero area", call. = FALSE)
  }
  list(xmin = xs[1], xmax = xs[2], ymin = ys[1], ymax = ys[2])
}

#' @keywords internal
#' @noRd
intersect_rect <- function(a, b) {
  list(
    xmin = max(a$xmin, b$xmin),
    xmax = min(a$xmax, b$xmax),
    ymin = max(a$ymin, b$ymin),
    ymax = min(a$ymax, b$ymax)
  )
}

#' Intersect a snapped box with the parent region union.
#' @keywords internal
#' @noRd
clip_region_to_omega <- function(reg, omega) {
  hits <- lapply(omega, function(o) intersect_rect(reg, o))
  ok <- vapply(hits, function(h) {
    h$xmax > h$xmin + 1e-9 && h$ymax > h$ymin + 1e-9
  }, logical(1L))
  if (!any(ok)) {
    return(NULL)
  }
  hits <- hits[ok]
  if (length(hits) == 1L) {
    return(hits[[1]])
  }
  area <- vapply(hits, function(h) {
    (h$xmax - h$xmin) * (h$ymax - h$ymin)
  }, numeric(1))
  hits[[which.max(area)]]
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
