#' Sum sparse triplet entries with duplicate coordinates
#'
#' @param i,j Integer index vectors (0-based).
#' @param x Numeric values.
#' @return List with components \code{i}, \code{j}, \code{x} (unique coordinates).
#' @keywords internal
#' @export
aggregate <- function(i, j, x) {
  if (!length(x)) {
    return(list(i = integer(), j = integer(), x = numeric()))
  }
  out <- stats::aggregate(
    stats::setNames(list(x), "x"),
    list(i = as.integer(i), j = as.integer(j)),
    sum
  )
  list(i = out$i, j = out$j, x = as.numeric(out$x))
}
