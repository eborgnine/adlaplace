#' Hierarchical tensor-product B-spline knot specifications
#'
#' Builds nested knot-line locations from a coarse raster box and optional
#' inner refinement extents. Knot **lines** do not depend on B-spline degree;
#' pass the result to [hb_basis()] or [matern()] to form the THB space.
#'
#' Each element of `inner` is one refinement level: a single extent, or a list
#' of extents for sibling boxes at that level. `list(ext)` and `list(list(ext))`
#' are the same one-level, one-box specification. `list(ext1, ext2)` is two
#' levels; siblings are `list(list(ext1, ext2))`. Extents are snapped to the
#' nested cell grid at that level (spacing `outer$resolution / prod(fact[1:i])`)
#' and clipped to the parent level.
#'
#' @param outer Level-0 box: `list(xmin, xmax, ymin, ymax, resolution)`, or a
#'   terra `SpatRaster` (converted to that list). `resolution` is `dx` or
#'   `c(dx, dy)`.
#' @param inner `NULL` (outer only), one extent, or a list of levels. An extent
#'   is `c(xmin, xmax, ymin, ymax)`, a terra `SpatExtent`, or
#'   `list(xmin, xmax, ymin, ymax)`.
#' @param fact Integer `>= 2`, recycled to `length(inner)`. `fact[i]` is the
#'   cell-size ratio from the previous level, so level-`i` resolution is
#'   `outer$resolution / prod(fact[1:i])`.
#' @return An `"hb_knots"` list, one element per level, each with `knots` and
#'   `rasters` (`xmin`, `xmax`, `ymin`, `ymax`, `resolution`).
#' @export
hb_knots <- function(outer, inner = NULL, fact = 2L) {
  if (inherits(outer, "hb_knots")) {
    return(outer)
  }
  spec0 <- coerce_outer_spec(outer)

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
  bp0 <- knots_from_spec(spec0)
  outer_box <- list(
    xmin = spec0$xmin,
    xmax = spec0$xmax,
    ymin = spec0$ymin,
    ymax = spec0$ymax
  )

  level_data <- vector("list", n_levels)
  level_data[[1L]] <- list(
    knots = list(x = unname(bp0$x), y = unname(bp0$y)),
    rasters = list(spec0)
  )

  bx_prev <- bp0$x
  by_prev <- bp0$y
  res_prev <- spec0$resolution

  for (i in seq_len(n_refine)) {
    res_i <- res_prev / as.numeric(fact[i])
    edges_i <- spec_cell_edges(raster_spec(outer_box, res_i))
    boxes <- lapply(inner_levels[[i]], function(reg) {
      snapped <- snap_extent_to_breaks(reg, edges_i$x, edges_i$y)
      clipped <- clip_region_to_omega(snapped, level_data[[i]]$rasters)
      if (is.null(clipped)) {
        stop(
          "refinement extent at inner level ", i,
          " must lie inside the parent domain",
          call. = FALSE
        )
      }
      clipped
    })
    specs <- lapply(boxes, function(reg) raster_spec(reg, res_i))
    new_x <- interior_breakpoints_from_specs(specs, axis = "x")
    new_y <- interior_breakpoints_from_specs(specs, axis = "y")
    bx <- sort(unique(c(bx_prev, new_x)))
    by <- sort(unique(c(by_prev, new_y)))
    level_data[[i + 1L]] <- list(
      knots = list(x = unname(bx), y = unname(by)),
      rasters = specs
    )
    bx_prev <- bx
    by_prev <- by
    res_prev <- res_i
  }

  structure(level_data, class = c("hb_knots", "list"))
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

  levels <- lapply(unclass(knots), function(lev) {
    open_x <- axis_to_open_knots(lev$knots$x, degree)
    open_y <- axis_to_open_knots(lev$knots$y, degree)
    lev$open_knots <- list(x = open_x, y = open_y)
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
      knots = knots
    ),
    class = c("hb_basis", "list")
  )
  hb$active <- hb_active_sets(hb, degree)
  hb$S <- hb_basis_map(hb, degree, active = hb$active)
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
is_named_extent <- function(x) {
  is.list(x) && all(c("xmin", "xmax", "ymin", "ymax") %in% names(x))
}

#' @keywords internal
#' @noRd
is_raster_spec <- function(x) {
  is_named_extent(x) &&
    "resolution" %in% names(x) &&
    is.numeric(x$resolution) &&
    length(x$resolution) %in% c(1L, 2L)
}

#' @keywords internal
#' @noRd
is_extent_like <- function(x) {
  inherits(x, "SpatExtent") ||
    inherits(x, "SpatRaster") ||
    inherits(x, "SpatVector") ||
    is_named_extent(x) ||
    (is.numeric(x) && length(x) == 4L && is.null(dim(x)))
}

#' @keywords internal
#' @noRd
extent_from_numeric4 <- function(x) {
  nms <- names(x)
  if (!is.null(nms) && all(c("xmin", "xmax", "ymin", "ymax") %in% nms)) {
    list(
      xmin = unname(x[["xmin"]]),
      xmax = unname(x[["xmax"]]),
      ymin = unname(x[["ymin"]]),
      ymax = unname(x[["ymax"]])
    )
  } else {
    list(
      xmin = unname(x[[1L]]),
      xmax = unname(x[[2L]]),
      ymin = unname(x[[3L]]),
      ymax = unname(x[[4L]])
    )
  }
}

#' @keywords internal
#' @noRd
validate_extent_region <- function(reg) {
  if (!is.finite(reg$xmin) || !is.finite(reg$xmax) ||
      !is.finite(reg$ymin) || !is.finite(reg$ymax) ||
      reg$xmax <= reg$xmin || reg$ymax <= reg$ymin) {
    stop(
      "each extent must have finite xmin < xmax and ymin < ymax",
      call. = FALSE
    )
  }
  reg
}

#' @keywords internal
#' @noRd
coerce_extent_region <- function(x) {
  if (is.numeric(x) && length(x) == 4L && is.null(dim(x))) {
    return(validate_extent_region(extent_from_numeric4(x)))
  }
  if (is_named_extent(x)) {
    return(validate_extent_region(list(
      xmin = unname(x$xmin),
      xmax = unname(x$xmax),
      ymin = unname(x$ymin),
      ymax = unname(x$ymax)
    )))
  }
  if (inherits(x, "SpatExtent")) {
    return(validate_extent_region(list(
      xmin = x$xmin,
      xmax = x$xmax,
      ymin = x$ymin,
      ymax = x$ymax
    )))
  }
  if (inherits(x, "SpatRaster") || inherits(x, "SpatVector")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(
        "terra is required to read SpatRaster or SpatVector extents",
        call. = FALSE
      )
    }
    e <- terra::ext(x)
    return(validate_extent_region(list(
      xmin = e$xmin,
      xmax = e$xmax,
      ymin = e$ymin,
      ymax = e$ymax
    )))
  }
  stop(
    "each inner region must be c(xmin, xmax, ymin, ymax), a SpatExtent, ",
    "or list(xmin, xmax, ymin, ymax)",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
coerce_outer_spec <- function(outer) {
  if (is_raster_spec(outer)) {
    spec <- raster_spec(
      list(
        xmin = outer$xmin,
        xmax = outer$xmax,
        ymin = outer$ymin,
        ymax = outer$ymax
      ),
      as.numeric(outer$resolution)
    )
    return(validate_raster_spec(spec))
  }
  if (inherits(outer, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("terra is required when outer is a SpatRaster", call. = FALSE)
    }
    e <- terra::ext(outer)
    spec <- raster_spec(
      list(xmin = e$xmin, xmax = e$xmax, ymin = e$ymin, ymax = e$ymax),
      c(terra::xres(outer), terra::yres(outer))
    )
    return(validate_raster_spec(spec))
  }
  stop(
    "outer must be list(xmin, xmax, ymin, ymax, resolution) or a SpatRaster",
    call. = FALSE
  )
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
raster_spec <- function(reg, resolution) {
  resolution <- unname(as.numeric(resolution))
  if (length(resolution) == 1L) {
    resolution <- c(resolution, resolution)
  }
  if (length(resolution) != 2L) {
    stop("resolution must be a positive numeric scalar or c(dx, dy)", call. = FALSE)
  }
  list(
    xmin = unname(reg$xmin),
    xmax = unname(reg$xmax),
    ymin = unname(reg$ymin),
    ymax = unname(reg$ymax),
    resolution = resolution
  )
}

#' @keywords internal
#' @noRd
validate_raster_spec <- function(spec) {
  spec <- raster_spec(spec, spec$resolution)
  validate_extent_region(spec)
  if (any(!is.finite(spec$resolution)) || any(spec$resolution <= 0)) {
    stop("resolution must be a positive numeric scalar or c(dx, dy)", call. = FALSE)
  }
  spec
}

#' @keywords internal
#' @noRd
spec_dims <- function(spec) {
  dx <- spec$resolution[1L]
  dy <- if (length(spec$resolution) > 1L) spec$resolution[2L] else spec$resolution[1L]
  list(
    dx = dx,
    dy = dy,
    ncol = max(1L, as.integer(round((spec$xmax - spec$xmin) / dx))),
    nrow = max(1L, as.integer(round((spec$ymax - spec$ymin) / dy)))
  )
}

#' @keywords internal
#' @noRd
spec_cell_edges <- function(spec) {
  d <- spec_dims(spec)
  list(
    x = spec$xmin + seq(0, d$ncol) * d$dx,
    y = spec$ymin + seq(0, d$nrow) * d$dy
  )
}

#' Unique axis knot lines from a raster spec (extent endpoints + cell centers).
#' @keywords internal
#' @noRd
knots_from_spec <- function(spec) {
  d <- spec_dims(spec)
  xs <- spec$xmin + (seq_len(d$ncol) - 0.5) * d$dx
  ys <- spec$ymin + (seq_len(d$nrow) - 0.5) * d$dy
  list(
    x = c(spec$xmin, xs[xs > spec$xmin & xs < spec$xmax], spec$xmax),
    y = c(spec$ymin, ys[ys > spec$ymin & ys < spec$ymax], spec$ymax)
  )
}

#' Level list from hb_knots (the object itself) or hb_basis$levels
#' @keywords internal
#' @noRd
hb_levels <- function(hb) {
  if (inherits(hb, "hb_basis")) {
    return(hb$levels)
  }
  unclass(hb)
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
interior_breakpoints_from_specs <- function(specs, axis = c("x", "y")) {
  axis <- match.arg(axis)
  vals <- numeric(0)
  for (spec in specs) {
    kl <- knots_from_spec(spec)
    lines <- kl[[axis]]
    inside <- if (axis == "x") {
      lines > spec$xmin & lines < spec$xmax
    } else {
      lines > spec$ymin & lines < spec$ymax
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
        "inner extent cannot be snapped to the nested cell grid",
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
  uk_c <- sort(unique(as.numeric(knots_coarse)))
  uk_f <- sort(unique(as.numeric(knots_fine)))
  # Sample every fine span; a coarse-only grid misses locally refined
  # functions and makes Bfine^T Bfine singular (e.g. fact >= 4).
  n_per <- degree + 3L
  n_span <- max(length(uk_f) - 1L, 1L)
  span_pts <- unlist(lapply(seq_len(n_span), function(s) {
    a <- uk_f[s]
    b <- uk_f[s + 1L]
    if (!is.finite(a) || !is.finite(b) || b <= a) {
      return(numeric(0))
    }
    seq(a, b, length.out = n_per + 2L)[-c(1L, n_per + 2L)]
  }), use.names = FALSE)
  x_eval <- sort(unique(c(
    greville_abscissae(knots_coarse, degree),
    greville_abscissae(knots_fine, degree),
    uk_c,
    uk_f,
    span_pts
  )))
  Bfine <- as.matrix(bspline_eval(knots_fine, x_eval, degree, 0L))
  Bcoarse <- as.matrix(bspline_eval(knots_coarse, x_eval, degree, 0L))
  # Bcoarse = Bfine %*% R; QR avoids squaring the condition number
  R <- qr.coef(qr(Bfine, LAPACK = TRUE), Bcoarse)
  if (anyNA(R)) {
    btB <- crossprod(Bfine)
    lam <- 1e-12 * max(1, mean(diag(btB)))
    R <- solve(btB + diag(lam, ncol(Bfine)), crossprod(Bfine, Bcoarse))
  }
  # Structural zeros come back as roundoff (~1e-17). Left in place they become
  # explicit nonzeros, making R -- and hence S -- dense in all but name.
  R[abs(R) < 1e-13 * max(1, max(abs(R)))] <- 0
  Matrix::drop0(methods::as(R, "dgCMatrix"))
}

#' Support bounds for every 1D basis index on an axis
#' @keywords internal
#' @noRd
bspline_support_bounds <- function(knots, degree, n) {
  degree <- as.integer(degree)
  i <- seq_len(n)
  list(left = knots[i + 1L], right = knots[i + degree + 1L])
}

#' Logical `n_x` by `n_y` matrix: support box inside the union of `omega`
#' @keywords internal
#' @noRd
supports_in_union <- function(sx, sy, omega, tol = 1e-6) {
  out <- matrix(FALSE, length(sx$left), length(sy$left))
  for (o in omega) {
    inx <- sx$left >= o$xmin - tol & sx$right <= o$xmax + tol
    iny <- sy$left >= o$ymin - tol & sy$right <= o$ymax + tol
    if (any(inx) && any(iny)) {
      out <- out | outer(inx, iny, "&")
    }
  }
  out
}

#' Active hierarchical basis indices per level
#'
#' Level 0 keeps every function not already covered by level 1; finer levels
#' keep functions whose support lies in that level's regions but not in the
#' next finer regions. Indices are column-major (`i + (j - 1) * n_x`).
#' @keywords internal
#' @noRd
hb_active_sets <- function(hb, degree) {
  degree <- as.integer(degree)
  levs <- hb_levels(hb)
  n_levels <- length(levs)
  active <- vector("list", n_levels)

  for (lev in seq_len(n_levels)) {
    kn <- levs[[lev]]$open_knots
    nx <- as.integer(levs[[lev]]$n_basis[["x"]])
    ny <- as.integer(levs[[lev]]$n_basis[["y"]])
    sx <- bspline_support_bounds(kn$x, degree, nx)
    sy <- bspline_support_bounds(kn$y, degree, ny)
    omega_lp1 <- if (lev < n_levels) levs[[lev + 1L]]$rasters else list()

    in_lp1 <- if (length(omega_lp1)) {
      supports_in_union(sx, sy, omega_lp1)
    } else {
      matrix(FALSE, nx, ny)
    }
    keep <- if (lev == 1L) {
      !in_lp1
    } else {
      supports_in_union(sx, sy, levs[[lev]]$rasters) & !in_lp1
    }

    idx <- which(keep)
    i <- ((idx - 1L) %% nx) + 1L
    j <- ((idx - 1L) %/% nx) + 1L
    active[[lev]] <- list(
      idx = idx,
      i = i,
      j = j,
      ij = .mapply(c, list(i, j), NULL),
      n_x = nx,
      n_y = ny
    )
  }
  active
}

#' Cumulative refinement (level l to finest) from pairwise matrices
#'
#' `mats[[lev]]` maps level `lev` to level `lev + 1`; the returned list holds
#' the product from each level up to the finest.
#' @keywords internal
#' @noRd
hb_cumulative_from_pairwise <- function(mats, n_finest) {
  n_levels <- length(mats) + 1L
  cum <- vector("list", n_levels)
  cum[[n_levels]] <- Matrix::Diagonal(n = n_finest)
  if (n_levels > 1L) {
    for (lev in seq(n_levels - 1L, 1L)) {
      cum[[lev]] <- cum[[lev + 1L]] %*% mats[[lev]]
    }
  }
  cum
}

#' Cumulative 1D refinement from level l to finest
#' @keywords internal
#' @noRd
hb_cumulative_refine_1d <- function(hb, axis = c("x", "y")) {
  axis <- match.arg(axis)
  levs <- hb_levels(hb)
  hb_cumulative_from_pairwise(
    hb_pairwise_refine_1d(hb, axis),
    levs[[length(levs)]]$n_basis[[axis]]
  )
}

#' Sparse basis map S: hierarchical coefficients to finest-level basis
#' @param truncate If `TRUE` (default), build truncated (THB) basis functions for
#'   a partition of unity; if `FALSE`, use the untruncated hierarchical embedding.
#' @keywords internal
#' @noRd
hb_basis_map <- function(hb, degree, truncate = TRUE, active = NULL) {
  degree <- as.integer(degree)
  if (is.null(active)) {
    active <- hb_active_sets(hb, degree)
  }
  levs <- hb_levels(hb)
  n_levels <- length(levs)
  nx_f <- as.integer(levs[[n_levels]]$n_basis[["x"]])
  ny_f <- as.integer(levs[[n_levels]]$n_basis[["y"]])

  refine_x <- hb_pairwise_refine_1d(hb, "x")
  refine_y <- hb_pairwise_refine_1d(hb, "y")
  cum_x <- hb_cumulative_from_pairwise(refine_x, nx_f)
  cum_y <- hb_cumulative_from_pairwise(refine_y, ny_f)

  n_col <- sum(vapply(active, function(a) length(a$idx), integer(1L)))
  if (!n_col) {
    stop("hierarchical basis has no active functions", call. = FALSE)
  }
  rows <- vector("list", n_col)
  vals <- vector("list", n_col)
  lens <- integer(n_col)
  col_id <- 0L

  for (lev in seq_len(n_levels)) {
    act <- active[[lev]]
    n_act <- length(act$idx)
    if (!n_act) {
      next
    }

    # Finest level: the hierarchical function *is* a finest tensor function.
    if (lev == n_levels) {
      for (k in seq_len(n_act)) {
        col_id <- col_id + 1L
        rows[[col_id]] <- act$idx[k]
        vals[[col_id]] <- 1
        lens[col_id] <- 1L
      }
      next
    }

    if (truncate) {
      # (Ty %x% Tx) vec(X) = vec(Tx X t(Ty)) with x fastest, and X = E_ij, so
      # the level-(l+1) coefficients are the outer product of two sparse
      # columns. Truncation zeroes the entries active at level l+1, then
      # cum_{l+1} lifts the remainder to the finest level. No Kronecker needed.
      Tx <- as_dgc_matrix(refine_x[[lev]])
      Ty <- as_dgc_matrix(refine_y[[lev]])
      Rx <- cum_x[[lev + 1L]]
      Ryt <- Matrix::t(cum_y[[lev + 1L]])
      nx1 <- as.integer(levs[[lev + 1L]]$n_basis[["x"]])
      ny1 <- as.integer(levs[[lev + 1L]]$n_basis[["y"]])
      truncated <- logical(nx1 * ny1)
      if (length(active[[lev + 1L]]$idx)) {
        truncated[active[[lev + 1L]]$idx] <- TRUE
      }
      for (k in seq_len(n_act)) {
        cx <- dgc_column(Tx, act$i[k])
        cy <- dgc_column(Ty, act$j[k])
        col_id <- col_id + 1L
        if (!length(cx$i) || !length(cy$i)) {
          rows[[col_id]] <- integer(0)
          vals[[col_id]] <- numeric(0)
          next
        }
        nxr <- length(cx$i)
        nyr <- length(cy$i)
        wi <- rep.int(cx$i, nyr)
        wj <- rep(cy$i, each = nxr)
        wv <- rep.int(cx$x, nyr) * rep(cy$x, each = nxr)
        keep <- !truncated[wi + (wj - 1L) * nx1]
        if (!any(keep)) {
          rows[[col_id]] <- integer(0)
          vals[[col_id]] <- numeric(0)
          next
        }
        W <- Matrix::sparseMatrix(
          i = wi[keep], j = wj[keep], x = wv[keep],
          dims = c(nx1, ny1)
        )
        out <- dgc_nonzeros(Rx %*% W %*% Ryt, nx_f)
        rows[[col_id]] <- out$idx
        vals[[col_id]] <- out$x
        lens[col_id] <- length(out$idx)
      }
      next
    }

    # Untruncated embedding: vec(Rx[, i] %o% Ry[, j]), x fastest.
    Rx <- as_dgc_matrix(cum_x[[lev]])
    Ry <- as_dgc_matrix(cum_y[[lev]])
    for (k in seq_len(n_act)) {
      cx <- dgc_column(Rx, act$i[k])
      cy <- dgc_column(Ry, act$j[k])
      col_id <- col_id + 1L
      nxr <- length(cx$i)
      nyr <- length(cy$i)
      if (!nxr || !nyr) {
        rows[[col_id]] <- integer(0)
        vals[[col_id]] <- numeric(0)
        next
      }
      rows[[col_id]] <- rep.int(cx$i, nyr) +
        (rep(cy$i, each = nxr) - 1L) * nx_f
      vals[[col_id]] <- rep.int(cx$x, nyr) * rep(cy$x, each = nxr)
      lens[col_id] <- nxr * nyr
    }
  }

  Matrix::sparseMatrix(
    i = unlist(rows, use.names = FALSE),
    j = rep.int(seq_len(n_col), lens),
    x = unlist(vals, use.names = FALSE),
    dims = c(nx_f * ny_f, n_col)
  )
}

#' Coerce any Matrix to dgCMatrix
#' @keywords internal
#' @noRd
as_dgc_matrix <- function(M) {
  if (methods::is(M, "dgCMatrix")) {
    return(M)
  }
  methods::as(
    methods::as(methods::as(M, "dMatrix"), "generalMatrix"),
    "CsparseMatrix"
  )
}

#' Nonzero rows and values of one dgCMatrix column
#' @keywords internal
#' @noRd
dgc_column <- function(M, j) {
  lo <- M@p[j]
  hi <- M@p[j + 1L]
  if (hi == lo) {
    return(list(i = integer(0), x = numeric(0)))
  }
  sel <- (lo + 1L):hi
  list(i = M@i[sel] + 1L, x = M@x[sel])
}

#' Column-major linear indices and values of all nonzeros in a sparse matrix
#' @keywords internal
#' @noRd
dgc_nonzeros <- function(M, n_row) {
  M <- as_dgc_matrix(M)
  counts <- diff(M@p)
  list(
    idx = (M@i + 1L) + (rep.int(seq_along(counts), counts) - 1L) * n_row,
    x = M@x
  )
}

#' Pairwise 1D refinement matrices between consecutive levels
#' @keywords internal
#' @noRd
hb_pairwise_refine_1d <- function(hb, axis = c("x", "y")) {
  axis <- match.arg(axis)
  levs <- hb_levels(hb)
  n_levels <- length(levs)
  degree <- hb$degree
  mats <- vector("list", n_levels - 1L)
  for (lev in seq_len(n_levels - 1L)) {
    coarse <- levs[[lev]]$open_knots[[axis]]
    fine <- levs[[lev + 1L]]$open_knots[[axis]]
    mats[[lev]] <- refine_matrix_1d(coarse, fine, degree)
  }
  mats
}

#' Map a sparse Gram from finest level to hierarchical basis
#' @keywords internal
#' @noRd
hb_project_gram <- function(M, S) {
  methods::as(Matrix::t(S) %*% M %*% S, "dgCMatrix")
}
