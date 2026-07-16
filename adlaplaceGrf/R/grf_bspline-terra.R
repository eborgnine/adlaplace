#' Resolve knots from a terra SpatRaster (grid lines)
#' @keywords internal
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
  # Cell centers strictly inside the extent become interior knot lines
  list(
    x = c(ext$xmin, xs[xs > ext$xmin & xs < ext$xmax], ext$xmax),
    y = c(ext$ymin, ys[ys > ext$ymin & ys < ext$ymax], ext$ymax)
  )
}

#' Coerce terra coords / knots then call grf_bspline_xy
#' @keywords internal
grf_bspline_terra <- function(coords, knots, degree = 2L) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "terra is required for SpatRaster / SpatVector input; install.packages(\"terra\")",
      call. = FALSE
    )
  }
  if (inherits(coords, "SpatRaster")) {
    xy <- terra::xyFromCell(coords, seq_len(terra::ncell(coords)))
    x <- xy[, 1]
    y <- xy[, 2]
  } else if (inherits(coords, "SpatVector")) {
    xy <- terra::crds(coords)
    x <- xy[, 1]
    y <- xy[, 2]
  } else {
    stop("internal: expected SpatRaster or SpatVector")
  }
  if (inherits(knots, "SpatRaster")) {
    knots <- knots_from_spatraster(knots, degree)
  }
  grf_bspline_xy(x, y, knots = knots, degree = degree)
}

#' Register terra S4 methods if terra is available
#' @keywords internal
register_terra_methods <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  methods::setMethod(
    "grf_bspline",
    methods::signature(coords = "SpatRaster"),
    function(coords, knots, degree = 2L, ...) {
      grf_bspline_terra(coords, knots, degree = degree)
    }
  )
  methods::setMethod(
    "grf_bspline",
    methods::signature(coords = "SpatVector"),
    function(coords, knots, degree = 2L, ...) {
      grf_bspline_terra(coords, knots, degree = degree)
    }
  )
  invisible(TRUE)
}
