#' Strip quotes from a term name string
#'
#' @param x Character vector of term names.
#' @return Character vector with surrounding quotes removed.
#' @keywords internal
strip_term_name <- function(x) {
  x <- as.character(x)
  gsub("^\"|\"$", "", x)
}

#' Parse Model Terms from Formula
#'
#' @description
#' Parses a formula and creates model terms using constructors from specified packages.
#'
#' @param formula Model formula with constructor terms (e.g. \code{iwp(x, ...)},
#'   \code{iid(fac)}, \code{nbinom(y, ...)} on the LHS). Bare symbols such as
#'   \code{x} are coerced to \code{linear(x)}. Symbols refer to column names in
#'   \code{data} passed to \code{model_data()}, not objects in the calling environment.
#' @param package Character vector of package names to search for model constructors
#' @param verbose print extra information
#' @return List of model term objects
#'
#' @export
collect_terms <- function(
  formula, package = character(0), verbose = FALSE
) {
  if (!methods::is(formula, "formula")) {
    warning("formula must be of class formula")
  }
  model_package <- unique(c(package, "adlaplace"))

  # Ensure packages are loaded
  pkg_env <- list()
  for (pkg in model_package) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(paste("Package", pkg, "not available, skipping"))
    }
    pkg_env[[pkg]] <- asNamespace(pkg)
  }

  formula_terms <- stats::terms(formula)
  term_labels <- rownames(attr(formula_terms, "factors")) # attr(stats::terms(formula), "term.labels")
  response_idx <- attr(formula_terms, "response")
  response_label <- if (response_idx > 0L) term_labels[response_idx] else NULL

  terms_1 <- lapply(term_labels, function(lab) {
    # a bare response symbol (e.g. y ~ x) defaults to a Gaussian
    # observation term; bare covariate symbols default to linear()
    if (identical(lab, response_label) && !grepl("[(]", lab)) {
      if (verbose) {
        message("Response '", lab, "' defaulting to gaussian(", lab, ")")
      }
      term_obj <- gaussian(lab)
      term_obj <- list(term_obj)
      names(term_obj) <- lab
      return(term_obj)
    }

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

  terms_1
}

#' @noRd
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
      if (i == 2L) {
        args[[i]] <- as.character(a)
      }
    } else if (is.call(a)) {
      args[[i]] <- coerce_term_call_symbols(a)
    }
  }
  as.call(args)
}

#' @noRd
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
