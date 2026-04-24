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
  model_fun <- get(model, envir = parent.frame(), mode = "function")
  model_fun(x_str, ...)
}
# my_env = new.env(parent = asNamespace("hpolcc"))
# my_env$f = function(...) {
#  adlaplace::f(...)
# }

#' Parse Model Terms from Formula
#'
#' @description
#' Parses a formula and creates model terms using constructors from specified packages.
#'
#' @param formula Model formula containing f() calls
#' @param package Character vector of package names to search for model constructors
#' @param verbose print extra information
#' @return List of model term objects
#'
#' @examples
#' # Parse formula using hpolcc models
#' terms <- collect_terms(y ~ f(x, model = "iwp"), package = "hpolcc")
#'
#' @export
collect_terms <- function(
  formula, package = character(0), verbose = FALSE
) {
  if (!methods::is(formula, "formula")) {
    warning("formula must be of class formula")
  }
  model_package <- unique(c(package, "adlaplace"))
  term_labels <- attr(stats::terms(formula), "term.labels")

  # Ensure packages are loaded
  pkg_env <- list()
  for (pkg in model_package) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(paste("Package", pkg, "not available, skipping"))
    }
    pkg_env[[pkg]] <- asNamespace(pkg)
  }


  terms_1 <- lapply(term_labels, function(lab) {
    term_obj <- NULL

    # Try each package's namespace until we find one that works

    for (pkg in names(pkg_env)) {
      try_result <- try(
        {
          eval(parse(text = lab), envir = pkg_env[[pkg]])
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

    if (is.null(term_obj)) {
      # term is probably a variable name, but if it has () in it, it's probably a model that wasn't found
      if (grepl("[(]", lab)) {
        warning(paste("Failed to parse term", lab, "in any of the specified packages"))
      }
      term_obj <- linear(lab)
    }
    if (!is.list(term_obj)) {
      term_obj <- list(term_obj)
      names(term_obj) <- lab
    }

    term_obj
  })

  terms_1 <- do.call(c, terms_1)

  # add intercept
  if (
    attr(stats::terms(formula), "intercept") &&
      !any(vapply(terms_1, inherits, logical(1), what = "intercept"))
  ) {
    terms_1$intercept <- intercept()
  }

  outcome_var <- all.vars(formula)[1]
  terms_1$response <- response(outcome_var)

  terms_1
}


make_empty_formula <- function(prefix = "", x) {
  stats::formula(paste0("~ 0 + ", prefix, x), env = new.env())
}
