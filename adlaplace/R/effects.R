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

#' Parse Model Terms from Formula
#'
#' @description
#' Parses a formula and creates model terms using constructors from specified packages.
#'
#' @param formula Model formula containing f() calls (e.g. \code{f(x, model = "iwp")},
#'   \code{linear(x)}). Symbols such as \code{x} refer to column names in
#'   \code{data} passed to \code{model_setup()}, not objects in the calling environment.
#' @param package Character vector of package names to search for model constructors
#' @param verbose print extra information
#' @return List of model term objects
#'
#' @examples
#' terms <- collect_terms(y ~ f(x, model = "iwp"))
#'
#' @keywords internal
coerce_term_call_symbols <- function(expr) {
  if (is.symbol(expr) || is.name(expr)) {
    return(call("linear", as.character(expr)))
  }
  if (!is.call(expr)) {
    return(expr)
  }
  args <- as.list(expr)
  for (i in seq_along(args)) {
    if (i == 1L) {
      next
    }
    a <- args[[i]]
    if (is.symbol(a) || is.name(a)) {
      args[[i]] <- as.character(a)
    } else if (is.call(a)) {
      args[[i]] <- coerce_term_call_symbols(a)
    }
  }
  as.call(args)
}

#' @keywords internal
eval_term_label <- function(lab, pkg_envs) {
  expr <- parse(text = lab, keep.source = FALSE)[[1]]
  expr <- coerce_term_call_symbols(expr)
  for (pkg in names(pkg_envs)) {
    try_result <- try(
      eval(expr, envir = pkg_envs[[pkg]]),
      silent = TRUE
    )
    if (!inherits(try_result, "try-error")) {
      return(try_result)
    }
  }
  NULL
}

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
    term_obj <- eval_term_label(lab, pkg_env)
    if (!is.null(term_obj) && verbose) {
      message(
        "Model term ", substr(lab, 1, 20),
        "... parsed from formula"
      )
    }

    if (is.null(term_obj)) {
      if (grepl("[(]", lab)) {
        stop(
          "Failed to parse term '", lab,
          "' in any of the specified packages",
          call. = FALSE
        )
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

