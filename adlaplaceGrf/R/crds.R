#' Extract point coordinates for spatial model terms
#'
#' Returns an `n` by `2` numeric matrix of `(x, y)` coordinates. Matrix and
#' data.frame inputs do not require **terra**. SpatVector points use
#' `terra::crds()` when terra is installed; non-point geometries error.
#'
#' @param x Point locations: a 2-column matrix/data.frame, or a terra
#'   `SpatVector` of points.
#' @param ... Unused.
#' @return Numeric matrix with columns `x` and `y`.
#' @export
crds <- function(x, ...) {
  UseMethod("crds")
}

#' @rdname crds
#' @export
crds.default <- function(x, ...) {
  stop(
    "crds() does not support class ", paste(class(x), collapse = "/"),
    call. = FALSE
  )
}

#' @rdname crds
#' @export
crds.matrix <- function(x, ...) {
  if (!is.numeric(x) || ncol(x) < 2L) {
    stop("matrix input to crds() must be numeric with at least 2 columns", call. = FALSE)
  }
  out <- cbind(x = as.numeric(x[, 1L]), y = as.numeric(x[, 2L]))
  if (any(!is.finite(out))) {
    stop("crds() coordinates must be finite", call. = FALSE)
  }
  out
}

#' @rdname crds
#' @export
crds.data.frame <- function(x, ...) {
  if (ncol(x) < 2L) {
    stop("data.frame input to crds() needs at least 2 columns", call. = FALSE)
  }
  nm <- names(x)
  if (all(c("x", "y") %in% nm)) {
    out <- cbind(x = as.numeric(x$x), y = as.numeric(x$y))
  } else {
    out <- cbind(x = as.numeric(x[[1L]]), y = as.numeric(x[[2L]]))
  }
  if (any(!is.finite(out))) {
    stop("crds() coordinates must be finite", call. = FALSE)
  }
  out
}

#' @rdname crds
#' @export
crds.SpatVector <- function(x, ...) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "terra is required for SpatVector coordinates; install.packages(\"terra\")",
      call. = FALSE
    )
  }
  geom <- terra::geomtype(x)
  if (!identical(geom, "points")) {
    stop(
      "matern()/crds() currently support SpatVector points only; got '",
      geom, "'",
      call. = FALSE
    )
  }
  xy <- terra::crds(x)
  cbind(x = as.numeric(xy[, 1L]), y = as.numeric(xy[, 2L]))
}
