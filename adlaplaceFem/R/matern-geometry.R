#' Resolve knot lines for matern()
#'
#' Accepts `list(x = ..., y = ...)` or a terra `SpatRaster`. For a raster,
#' axis knot lines are the extent endpoints plus interior cell centers.
#'
#' @param knots Knot specification.
#' @return `list(x = ..., y = ...)` of numeric knot-line positions.
#' @keywords internal
#' @noRd
matern_knots <- function(knots) {
  if (inherits(knots, "SpatRaster")) {
    return(knots_from_spatraster(knots, degree = 2L))
  }
  if (!is.list(knots) || is.null(knots$x) || is.null(knots$y)) {
    stop(
      "knots must be list(x = ..., y = ...) or a SpatRaster",
      call. = FALSE
    )
  }
  knots
}

#' Extract 2D coordinates from a data-column geometry
#'
#' Supports an `n` by `2` numeric matrix/data.frame, a list of length-2
#' numerics, WKT text (via terra), or HEX WKB Point strings (via terra).
#'
#' @param x Column contents.
#' @return Numeric matrix with columns `x` and `y`.
#' @keywords internal
#' @noRd
matern_column_xy <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (is.matrix(x)) {
    if (!is.numeric(x) || ncol(x) != 2L) {
      stop(
        "geometry matrix must be numeric with exactly 2 columns",
        call. = FALSE
      )
    }
    out <- cbind(x = as.numeric(x[, 1L]), y = as.numeric(x[, 2L]))
  } else if (is.list(x) && !is.object(x)) {
    if (!length(x) || !all(vapply(x, function(z) {
      is.numeric(z) && length(z) >= 2L
    }, logical(1)))) {
      stop(
        "geometry list-column entries must be numeric with at least 2 values",
        call. = FALSE
      )
    }
    out <- cbind(
      x = vapply(x, function(z) as.numeric(z)[1L], numeric(1)),
      y = vapply(x, function(z) as.numeric(z)[2L], numeric(1))
    )
  } else if (is.character(x)) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(
        "terra is required to parse WKT/HEX geometries; install.packages(\"terra\")",
        call. = FALSE
      )
    }
    is_wkt <- grepl(
      "^(POINT|MULTIPOINT|LINESTRING|POLYGON|MULTI|GEOMETRY)",
      x,
      ignore.case = TRUE
    )
    if (all(is_wkt)) {
      v <- terra::vect(x)
    } else {
      # HEX WKB points (e.g. terra::as.data.frame(..., geom = "HEX"))
      wkb <- lapply(x, function(hex) {
        hex <- gsub("\\s+", "", hex)
        if (nchar(hex) %% 2L != 0L) {
          stop("HEX geometry string has odd length", call. = FALSE)
        }
        nib <- strsplit(hex, "", fixed = TRUE)[[1L]]
        as.raw(as.integer(matrix(strtoi(nib, 16L), ncol = 2L, byrow = TRUE) %*%
          c(16L, 1L)))
      })
      v <- terra::vect(wkb)
    }
    if (!identical(terra::geomtype(v), "points")) {
      stop("matern() currently supports point geometries only", call. = FALSE)
    }
    xy <- terra::crds(v)
    out <- cbind(x = as.numeric(xy[, 1L]), y = as.numeric(xy[, 2L]))
  } else {
    stop(
      "geometry column must be a 2-column matrix/data.frame, list of ",
      "coordinates, WKT text, or HEX WKB points; got ",
      paste(class(x), collapse = "/"),
      call. = FALSE
    )
  }
  if (any(!is.finite(out))) {
    stop("matern() coordinates must be finite", call. = FALSE)
  }
  out
}
