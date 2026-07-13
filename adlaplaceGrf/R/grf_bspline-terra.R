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
  # Use cell-edge lines: unique sorted coordinates as interior+boundary seeds
  list(
    x = open_knot_vector(c(ext$xmin, ext$xmax), degree, interior = xs[xs > ext$xmin & xs < ext$xmax]),
    y = open_knot_vector(c(ext$ymin, ext$ymax), degree, interior = ys[ys > ext$ymin & ys < ext$ymax])
  )
}

#' Coerce terra coords / knots then call grf_bspline_xy
#' @keywords internal
grf_bspline_terra <- function(coords, degree, knots, n_interior = 5L) {
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
  grf_bspline_xy(x, y, degree = degree, knots = knots, n_interior = n_interior)
}

#' Register terra S4 methods if terra is available
#' @keywords internal
register_terra_methods <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  # SpatRaster / SpatVector are S4 classes from terra
  methods::setMethod(
    "grf_bspline",
    signature(coords = "SpatRaster"),
    function(coords, degree = 2L, knots = NULL, n_interior = 5L, ...) {
      grf_bspline_terra(coords, degree, knots, n_interior = n_interior)
    }
  )
  methods::setMethod(
    "grf_bspline",
    signature(coords = "SpatVector"),
    function(coords, degree = 2L, knots = NULL, n_interior = 5L, ...) {
      grf_bspline_terra(coords, degree, knots, n_interior = n_interior)
    }
  )
  invisible(TRUE)
}
