#' By-Group Classes and Functions
#'
#' @description Functions for creating and managing by_group objects for hierarchical
#' model grouping information.
#'
#' @name by_group
#' @rdname by_group
#' @export
NULL

#' @rdname by_group
#' @param term Character, the grouping variable name
#' @param data A data frame. If provided and levels is missing, unique values are
#'   extracted from data[[term]]
#' @param levels Character vector of unique levels for the grouping variable.
#'   If missing and data is provided, unique values are extracted from data[[term]]
#' @return A by_group object
#' @examples
#' # Create with explicit levels
#' bg <- by_group(term = "group", levels = c("A", "B", "C"))
#'
#' # Create from data
#' dat <- data.frame(group = c("A", "B", "A", "C"))
#' bg <- by_group(term = "group", data = dat)
by_group <- function(term, data = NULL, levels = NULL) {
  if (missing(term)) {
    # Return empty by_group object for use in prototypes
    return(methods::new("by_group", term = character(0), levels = character(0), labels = character(0)))
  }

  if (is.null(levels)) {
    if (is.null(data)) {
      return(methods::new("by_group",
        term = term,
        levels = character(0),
        labels = character(0)
      ))
    }
    if (!(term %in% names(data))) {
      stop("Grouping variable '", term, "' not found in data")
    }
    unique_values <- unique(data[[term]])
    is_numeric <- is.numeric(unique_values)
    levels <- as.character(unique_values)
  } else {
    is_numeric <- is.numeric(levels)
  }

  # Construct labels from levels with formatC for numeric
  if (is_numeric) {
    labels <- formatC(
      as.numeric(levels),
      width = ceiling(log10(max(abs(as.numeric(levels))))),
      flag = "0"
    )
  } else {
    labels <- as.character(levels)
  }

  methods::new("by_group",
    term = term,
    levels = levels,
    labels = labels
  )
}

#' @rdname by_group
#' @description Updates the by_group levels in a model term from data.
#' If the term has no by slot, returns it unchanged.
#' If levels are empty, they are populated from data.
#' If levels exist, new values in data are added with a warning.
#'
#' @param term A model term object (inherited from model class)
#' @param data A data frame containing the grouping variable
#' @return The term with updated by_group levels
#' @export
#' @examples
#' # Assuming hrp is an hrpoly term with a by slot
#' hrp <- hrpoly(x = "age", by = "site")
#' dat <- data.frame(age = 1:10, site = rep(c("A", "B"), each = 5))
#' hrp <- add_by_levels(hrp, dat)
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
    stop("Grouping variable '", by_var, "' not found in data")
  }

  data_levels <- as.character(unique(data[[by_var]]))

  if (length(term@by@levels) == 0) {
    # Create new by_group with levels from data
    unique_values <- unique(data[[by_var]])
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
    term@by <- methods::new("by_group",
      term = term@by@term,
      levels = levels,
      labels = labels
    )
  } else {
    # Check if all data levels are in existing levels
    missing_levels <- setdiff(data_levels, term@by@levels)
    if (length(missing_levels) > 0) {
      warning("Adding new levels to by_group: ", paste(missing_levels, collapse = ", "))
      # Determine if existing levels are numeric
      is_numeric <- is.numeric(term@by@levels)
      new_levels <- c(term@by@levels, missing_levels)
      if (is_numeric) {
        new_labels <- formatC(
          as.numeric(new_levels),
          width = ceiling(log10(max(abs(as.numeric(new_levels))))),
          flag = "0"
        )
      } else {
        new_labels <- as.character(new_levels)
      }
      term@by <- methods::new("by_group",
        term = term@by@term,
        levels = new_levels,
        labels = new_labels
      )
    }
  }

  term
}
