#' Coerce terra coords / knots then call fem_bspline_xy
#' @keywords internal
#' @noRd
fem_bspline_terra <- function(coords, knots, degree = 2L) {
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
  fem_bspline_xy(x, y, knots = knots, degree = degree)
}

#' Register terra S4 methods if terra is available
#' @keywords internal
#' @noRd
register_terra_methods <- function() {
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  methods::setMethod(
    "fem_bspline",
    methods::signature(coords = "SpatRaster"),
    function(coords, knots, degree = 2L, ...) {
      fem_bspline_terra(coords, knots, degree = degree)
    }
  )
  methods::setMethod(
    "fem_bspline",
    methods::signature(coords = "SpatVector"),
    function(coords, knots, degree = 2L, ...) {
      fem_bspline_terra(coords, knots, degree = degree)
    }
  )
  invisible(TRUE)
}
