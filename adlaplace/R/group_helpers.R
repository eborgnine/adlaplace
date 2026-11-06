#' @include classes.R
NULL

#' Column name stored in a \code{by} slot (\code{by_group} or character).
#' @keywords internal
by_column_name <- function(term) {
  if (!methods::is(term, "model_term") || !"by" %in% methods::slotNames(term)) {
    return(character(0))
  }
  by <- methods::slot(term, "by")
  if (inherits(by, "by_group")) {
    nm <- by@term
    if (length(nm) == 1L && is.character(nm) && nzchar(nm)) {
      return(nm)
    }
    return(character(0))
  }
  if (is.character(by) && length(by) == 1L && nzchar(by)) {
    return(by)
  }
  character(0)
}

#' Fill or check the grouping levels of a model term
#'
#' Terms with a \code{by} slot of class \code{by_group} get their levels from
#' \code{data} when none were supplied at construction. When levels were
#' supplied, every non-missing value in the data must already be one of them.
#' Terms without a \code{by_group} slot are returned unchanged.
#'
#' Missing values in the grouping column are an error. A factor column with
#' no levels stored on the term uses \code{levels()} (including unused levels);
#' any other column uses \code{unique()} in data order.
#'
#' @param term A model term.
#' @param data A data frame containing the grouping column.
#' @return \code{term}, with \code{@by} levels and labels filled when they
#'   were empty.
#' @export
add_by_levels <- function(term, data) {
  if (!"by" %in% methods::slotNames(term)) {
    return(term)
  }
  if (!inherits(term@by, "by_group")) {
    return(term)
  }

  by_var <- term@by@term
  if (length(by_var) == 0 || !is.character(by_var)) {
    return(term)
  }
  if (!(by_var %in% names(data))) {
    stop("Grouping variable '", by_var, "' not found in data", call. = FALSE)
  }

  x <- data[[by_var]]
  if (anyNA(x)) {
    stop(
      "grouping variable '", by_var, "' has missing values",
      call. = FALSE
    )
  }

  if (length(term@by@levels) == 0) {
    if (is.factor(x)) {
      levels <- levels(x)
      labels <- as.character(levels)
    } else {
      unique_values <- unique(x)
      is_numeric <- is.numeric(unique_values)
      levels <- as.character(unique_values)
      if (is_numeric) {
        labels <- formatC(
          as.numeric(levels),
          width = ceiling(log10(max(abs(as.numeric(levels))))),
          flag = "0"
        )
      } else {
        labels <- as.character(levels)
      }
    }
    term@by <- methods::new(
      "by_group",
      term = term@by@term,
      levels = levels,
      labels = labels
    )
    return(term)
  }

  data_levels <- as.character(unique(x))
  missing_levels <- setdiff(data_levels, term@by@levels)
  if (length(missing_levels) > 0) {
    stop(
      "grouping variable '", by_var, "' has levels not in the term: ",
      paste(missing_levels, collapse = ", "),
      call. = FALSE
    )
  }
  term
}

#' Row index of each observation in a set of group levels
#'
#' @param data A data frame.
#' @param by Grouping column name.
#' @param levels Character vector of allowed levels, in block order.
#' @return Integer vector of length \code{nrow(data)}, with values in
#'   \code{1:length(levels)}.
#' @export
group_index <- function(data, by, levels) {
  by <- as.character(by)[1L]
  if (!length(by) || is.na(by) || !by %in% names(data)) {
    stop("grouping variable '", by, "' not found in data", call. = FALSE)
  }
  x <- data[[by]]
  if (anyNA(x)) {
    stop("grouping variable '", by, "' has missing values", call. = FALSE)
  }
  levels <- as.character(levels)
  group <- match(as.character(x), levels)
  if (anyNA(group)) {
    unknown <- unique(as.character(x)[is.na(group)])
    stop(
      "grouping variable '", by, "' has levels not in the term: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  group
}

#' Scatter design columns into one block per group
#'
#' Column \code{j} of \code{B} for a row in group \code{g} is written to
#' column \code{(g - 1) * ncol(B) + j}. Basis columns stay together inside
#' each group block.
#'
#' @param B Design matrix (base \code{matrix} or any \code{Matrix}).
#' @param group Integer group index per row, from [group_index()].
#' @param n_levels Number of group blocks.
#' @return A \code{dgCMatrix} with \code{ncol(B) * n_levels} columns.
#'   Column names are left unset.
#' @export
scatter_by_group <- function(B, group, n_levels) {
  n_levels <- as.integer(n_levels)[1L]
  if (length(n_levels) != 1L || is.na(n_levels) || n_levels < 1L) {
    stop("n_levels must be a positive integer", call. = FALSE)
  }
  if (!is.matrix(B) && !methods::is(B, "Matrix")) {
    stop("B must be a matrix", call. = FALSE)
  }
  if (length(group) != nrow(B)) {
    stop(
      "length(group) (", length(group), ") must equal nrow(B) (", nrow(B), ")",
      call. = FALSE
    )
  }
  group <- as.integer(group)
  if (anyNA(group) || any(group < 1L | group > n_levels)) {
    stop("group indices must be between 1 and n_levels", call. = FALSE)
  }

  if (!inherits(B, "dgCMatrix")) {
    if (!methods::is(B, "Matrix")) {
      B <- Matrix::Matrix(B, sparse = TRUE)
    }
    if (methods::is(B, "diagonalMatrix")) {
      B <- methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix")
    }
    B <- methods::as(methods::as(B, "generalMatrix"), "CsparseMatrix")
    B <- methods::as(B, "dMatrix")
    B <- methods::as(B, "generalMatrix")
    if (!inherits(B, "dgCMatrix")) {
      B <- methods::as(B, "dgCMatrix")
    }
  }

  nb <- ncol(B)
  counts <- diff(B@p)
  cols <- rep.int(seq_along(counts) - 1L, counts)
  rows <- B@i
  col_out <- cols + (group[rows + 1L] - 1L) * nb
  Matrix::sparseMatrix(
    i = rows + 1L,
    j = col_out + 1L,
    x = B@x,
    dims = c(nrow(B), nb * n_levels),
    repr = "C"
  )
}

#' Gamma labels for a grouped design, basis fastest within each level
#'
#' @param base_labels Character vector of labels for one group's columns.
#' @param by_labels Character vector of group labels, in block order.
#' @return Character vector of length \code{length(base_labels) * length(by_labels)},
#'   each entry \code{paste0(base, "_g", group)}.
#' @export
group_gamma_labels <- function(base_labels, by_labels) {
  base_labels <- as.character(base_labels)
  by_labels <- as.character(by_labels)
  n <- length(base_labels)
  paste0(
    rep(base_labels, times = length(by_labels)),
    "_g",
    rep(by_labels, each = n)
  )
}
