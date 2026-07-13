#' Per-shard layout for AD density evaluation
#'
#' @param x Character vector of term names.
#' @return Character vector with surrounding quotes removed.
#' @export
strip_term_name <- function(x) {
  x <- as.character(x)
  gsub("^\"|\"$", "", x)
}

#' Parse Model Terms from Formula
#'
#' @description
#' Parses a formula and creates model terms by evaluating constructor calls in
#' the formula environment (so \code{pkg::term(...)} and terms from attached
#' packages on the search path work). Falls back to the \pkg{adlaplace}
#' namespace for built-in constructors.
#'
#' @param formula Model formula with constructor terms (e.g. \code{iwp(x, ...)},
#'   \code{iid(fac)}, \code{nbinom(y, ...)} on the LHS). Bare symbols such as
#'   \code{x} are coerced to \code{linear(x)}. Unnamed first positional symbols
#'   are treated as column names in \code{data}; named arguments (e.g.
#'   \code{x = sites} for spatial terms) are evaluated as objects in the
#'   formula environment. The named argument \code{size} is also treated as a
#'   column name when passed as a symbol (\code{binomial(y, size = N)}).
#' @param verbose print extra information
#' @return List of model term objects
#'
#' @export
collect_terms <- function(formula, verbose = FALSE) {
  if (!methods::is(formula, "formula")) {
    warning("formula must be of class formula")
  }

  formula_env <- environment(formula)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  formula_terms <- stats::terms(formula)
  term_labels <- rownames(attr(formula_terms, "factors"))
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

    term_obj <- eval_term_label(lab, formula_env)
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
          "'. Use pkg::term(...) or library(pkg) so constructors are visible.",
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

#' Column-name formals coerced even when named (symbols -> character).
#' @noRd
.column_name_args <- c("size")

#' @noRd
coerce_term_call_symbols <- function(expr) {
  if (is.symbol(expr) || is.name(expr)) {
    return(call("linear", as.character(expr)))
  }
  if (!is.call(expr)) {
    return(expr)
  }
  args <- as.list(expr)
  nm <- names(args)
  for (i in seq_along(args)) {
    if (i == 1L) {
      next
    }
    a <- args[[i]]
    arg_name <- if (is.null(nm)) {
      ""
    } else {
      nm[[i]]
    }
    if (is.null(arg_name) || is.na(arg_name)) {
      arg_name <- ""
    }
    if (is.symbol(a) || is.name(a)) {
      is_unnamed <- !nzchar(arg_name)
      # Only the first positional (unnamed) argument is treated as a column
      # name, so named object args like matern(x = sites) stay objects.
      # Named size=N is also a column name.
      if ((is_unnamed && i == 2L) || arg_name %in% .column_name_args) {
        args[[i]] <- as.character(a)
      }
    } else if (is.call(a)) {
      args[[i]] <- coerce_term_call_symbols(a)
    }
  }
  as.call(args)
}

#' @noRd
eval_term_label <- function(lab, formula_env) {
  expr <- parse(text = lab, keep.source = FALSE)[[1]]
  expr <- coerce_term_call_symbols(expr)
  try_result <- try(
    eval(expr, envir = formula_env),
    silent = TRUE
  )
  if (!inherits(try_result, "try-error")) {
    return(try_result)
  }
  try_result <- try(
    eval(expr, envir = asNamespace("adlaplace")),
    silent = TRUE
  )
  if (!inherits(try_result, "try-error")) {
    return(try_result)
  }
  NULL
}
