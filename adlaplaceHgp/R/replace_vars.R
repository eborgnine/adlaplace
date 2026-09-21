#' @importFrom adlaplace replace_vars
#' @export
setMethod("replace_vars", "hiwp", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  x <- methods::callNextMethod(x, map)
  if (length(x@by@term)) {
    x@by@term <- vapply(
      x@by@term,
      function(nm) {
        if (nm %in% names(map)) as.character(map[[nm]]) else nm
      },
      character(1)
    )
  }
  x
})

#' @importFrom adlaplace replace_vars
#' @export
setMethod("replace_vars", "rsiid", function(x, map) {
  if (is.null(map) || !length(map)) {
    return(x)
  }
  map <- stats::setNames(as.character(map), names(map))
  # Grouping factor lives in @name; exposure multiplier in @mult
  new_name <- if (x@name %in% names(map)) as.character(map[[x@name]]) else x@name
  new_mult <- if (length(x@mult) && x@mult %in% names(map)) {
    as.character(map[[x@mult]])
  } else {
    x@mult
  }
  x@name <- new_name
  x@mult <- new_mult
  x@label <- paste(c(new_name, new_mult, "rsiid"), collapse = "_")
  x@formula <- stats::as.formula(paste0("~ 0 + ", new_name), env = new.env())
  x
})
