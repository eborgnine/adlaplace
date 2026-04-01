#' Create Model Terms
#'
#' @description
#' Creates model terms for use in hierarchical models.
#'
#' @param x Variable name (or expression)
#' @param model Model type (e.g., "iid", "iwp", "fpoly", etc.)
#' @param ... Additional arguments passed to the model constructor
#'
#' @return A model term object
#'
#' @examples
#' # Create an IWP term
#' term <- f(x = 0:10, model = "iwp", ref_value = 5, knots = seq(0, 10, 2))
#'
#' @export
f <- function(x, model = "iid", ...) {
  x_str <- deparse(substitute(x))
  model_fun <- get(model, envir = parent.frame(), mode="function")
  model_fun(x_str, ...)
}
#my_env = new.env(parent = asNamespace("hpolcc"))
#my_env$f = function(...) {
#  adlaplace::f(...)
#}
#eval(parse(text="f(x , model = \"hiwp\")"), envir = my_env)

#' Parse Model Terms from Formula
#'
#' @description
#' Parses a formula and creates model terms using constructors from specified packages.
#'
#' @param formula Model formula containing f() calls
#' @param model_package Character vector of package names to search for model constructors
#' @return List of model term objects
#'
#' @examples
#' # Parse formula using hpolcc models
#' terms <- collect_terms(y ~ f(x, model = "iwp"), model_package = "hpolcc")
#'
#' @export
collect_terms <- function(formula, package = character(0), verbose = FALSE) {
  model_package <- unique(c(package, "adlaplace"))
  term_labels <- attr(stats::terms(formula), "term.labels")

  # Ensure packages are loaded
  for (pkg in model_package) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(paste("Package", pkg, "not available, skipping"))
    }
  }

  terms_1 <- lapply(term_labels, function(lab) {
    if (grepl("f\\(", lab)) {
      term_obj <- NULL

      # Try each package's namespace until we find one that works
      for (pkg in model_package) {
        pkg_env <- tryCatch(
          {
            asNamespace(pkg)
          },
          error = function(e) NULL
        )

        if (!is.null(pkg_env)) {
          try_result <- try(
            {
              eval(parse(text = lab), envir = pkg_env)
            },
            silent = TRUE
          )

          if (!inherits(try_result, "try-error")) {
            term_obj <- try_result
            if (verbose) {
              message(
                paste(
                  "Model term", substr(lab, 1, 20),
                  "... found in package", pkg
                )
              )
            }
            break
          }
        }
      }

      if (is.null(term_obj)) {
        warning(paste("Failed to parse term", lab, "in any of the specified packages"))
        # Fallback to linear term
        return(list(lab))
      }

      return(term_obj)
    } else {
      result <- list(linear(lab))
      names(result) <- result[[1]]@term
      return(result)
    }
  })

  do.call(c, terms_1)
}

#' Create Formula with Empty Environment
#'
#' @description
#' Creates a formula with a new empty environment, useful for avoiding
#' unintended variable capture from the calling environment.
#'
#' @param prefix Optional prefix to add to variables
#' @param x Variable name
#'
#' @return A formula object with empty environment
#'
#' @examples
#' # Create formula with empty environment
#' make_empty_formula("my_", "var")
#'
#' @export
make_empty_formula <- function(prefix = "", x) {
  formula(paste0("~ 0 + ", prefix, x), env = new.env())
}
