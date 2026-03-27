#' Helpers for the `hnlm` function
#' @name hnlm_helpers
#' @rdname hnlm_helpers
#'
#' @description
#' Common setup helpers used by `hnlm()`. These functions collect and
#' augment term specifications, construct design matrices, and prepare
#' gamma/theta metadata.
#'
#' The helper functions include:
#' - `collect_terms()`: set up the terms object used throughout fitting.
#' - `get_extra()`: construct effect-specific quantities.
#' - `get_design()`: construct the design matrix for a given effect.
#' - `get_gamma_setup()`: set up gamma metadata for an effect.
#' - `get_precision()`: construct a precision matrix for a random effect.
#' - `get_theta_setup()`: set up theta metadata for an effect.
#' - `add_fpoly()`: add fixed polynomial terms implied by an effect.
#' - `add_rpoly()`: add random polynomial terms implied by an effect.
NULL

#' @rdname hnlm_helpers
collect_terms <- function(formula) {
  term_labels <- attr(terms(formula), "term.labels")

  terms_1 <- lapply(term_labels, function(lab) {
    if (grepl("f\\(", lab)) {
      eval(parse(text = lab))
    } else {
      linear(lab)
    }
  })
  terms_all = do.call(c, terms_1)
  terms_all
}

#' @rdname hnlm_helpers
get_extra <- function(term, data, cc_matrix) {
  list2env(term, envir = environment())

  if (term$run_as_is) {
    return(term)
  }

  if (term$model == "iid") {
    term$n <- length(unique(data[[term$var]]))
  }

  if (!is.null(term$group_var)) {
    term$groups <- factor(unique(data[[term$group_var]])) |> droplevels()
    if (any(term$groups == "GLOBAL")) {
      stop(
        "Please change value (name) of GLOBAL in the group variable for ",
        term$var,
        " -- ",
        term$model
      )
    }
    term$groups <- factor(
      term$groups,
      levels = c("GLOBAL", levels(term$groups))
    )
    term$ngroups <- length(term$groups)
  }

  if (!is.null(term$knots)) {
    term$nknots <- length(term$knots)
  }

  if (!(term$var %in% names(data))) {
    warning("cant find ", term$var, " in data")
  }
  term$range <- range(data[[term$var]])

  term
}

#' @rdname hnlm_helpers
get_design <- function(term, data) {
  list2env(term, envir = environment())

  if (term$model == "iwp") {
    design(term, data)
  } else if (term$model == "hiwp") {
    design(term, data)
  } else if (term$model == "iid") {
    iidDesign(term, data)
  } else if (term$model == "rpoly") {
    rpolyDesign(term, data)
  } else if (term$model == "fpoly") {
    fpolyDesign(term, data)
  } else if (term$model == "hrpoly") {
    hrpolyDesign(term, data)
  } else if (term$model == "hfpoly") {
    stop("hfpoly design is not implemented")
  } else {
    stop("Unknown term model (", term$model, ")")
  }
}

#' @rdname hnlm_helpers
get_gamma_setup <- function(term) {
  list2env(term, envir = environment())

  if (is.null(term$groups)) {
    term$groups <- NA
  }
  if (is.null(term$knots)) {
    term$knots <- NA
  }

  if (term$model %in% c("fpoly", "hrpoly")) {
    term$order <- seq_len(term$p)
  } else {
    term$order <- NA
  }

  if (term$model %in% c("iwp", "hiwp")) {
    term$basis <- seq(1, len = term$nknots - 1)
    term$order <- NA
  } else {
    term$basis <- NA
  }
  if (term$model == "hiwp" && term$include_global) {
    term$groups <- c(levels(term$groups)[1], as.character(term$groups))
  }

  result <- expand.grid(
    var = term$var,
    model = term$model,
    group = term$groups,
    basis = term$basis,
    order = term$order
  )
  if (!all(is.na(result$basis))) {
    result$basis <- formatC(
      result$basis,
      width = max(ceiling(c(1, log10(result$basis))), na.rm = TRUE),
      flag = "0"
    )
  }
  the_na <- apply(result, 2, function(xx) all(is.na(xx)))
  result$name <- apply(
    result[, !the_na, drop = FALSE],
    1,
    paste,
    collapse = "_"
  )
  result
}

#' @rdname hnlm_helpers
get_precision <- function(term) {
  list2env(term, envir = environment())

  if (term$model == "iwp") {
    iwpPrecision(term)
  } else if (term$model == "hiwp") {
    hiwpPrecision(term)
  } else if (term$model == "iid") {
    iidPrecision(term)
  } else if (term$model == "rpoly") {
    rpolyPrecision(term)
  } else if (term$model == "hrpoly") {
    hrpolyPrecision(term)
  } else {
    stop("Unknown term model (", term$model, ")")
  }
}

#' @rdname hnlm_helpers
get_theta_setup <- function(theta_info, term) {
  if (term$model == "iwp") {
    iwpTheta(theta_info, term)
  } else if (term$model == "hiwp") {
    hiwpTheta(theta_info, term)
  } else if (term$model == "iid") {
    iidTheta(theta_info, term)
  } else if (term$model == "rpoly") {
    rpolyTheta(theta_info, term)
  } else if (term$model == "hrpoly") {
    hrpolyTheta(theta_info, term)
  } else if (term$model == "fpoly") {
    fpolyTheta(theta_info, term)
  } else {
    stop("Unknown term model (", term$model, ")")
  }
}

#' @rdname hnlm_helpers
add_fpoly <- function(term) {
  new_terms <- NULL

  if (!is.null(term$fpoly_p) && term$fpoly_p > 0) {
    new_f <- paste0(
      "~ f(__var__, model = 'fpoly', ref_value = __rv__, p = __p__, ",
      "boundary_is_random = __bound__)"
    ) |>
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$fpoly_p) |>
      gsub(
        pattern = "__bound__",
        replacement = c(term$boundary_is_random, TRUE)[1]
      ) |>
      as.formula()
    new_terms[[length(new_terms) + 1]] <- collect_terms(new_f)[[1]]
  } else {
    warning("fpoly_p must be positive")
  }

  new_terms
}

#' @rdname hnlm_helpers
add_rpoly <- function(term) {
  new_terms <- NULL

  if (!is.null(term$rpoly_p) && term$rpoly_p > 0) {
    new_f <- "~ f(__var__, model = 'rpoly', ref_value = __rv__, p = __p__)" |>
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$rpoly_p) |>
      as.formula()
    new_terms[[length(new_terms) + 1]] <- collect_terms(new_f)[[1]]
  }

  if (!is.null(term$hrpoly_p) && term$hrpoly_p > 0) {
    new_f <- paste0(
      "~ f(__var__, model = 'hrpoly', ref_value = __rv__, p = __p__, ",
      "group_var = __gv__, init = __init__,lower = __lower__,upper = __upper__, include_global = F)"
    ) |>
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$hrpoly_p) |>
      gsub(pattern = "__init__", replacement = deparse(term$init_hrpoly)) |>
      gsub(pattern = "__lower__", replacement = deparse(term$lower_hrpoly)) |>
      gsub(pattern = "__upper__", replacement = deparse(term$upper_hrpoly)) |>
      gsub(pattern = "__gv__", replacement = term$group_var) |>
      as.formula()

    new_terms[[length(new_terms) + 1]] <- collect_terms(new_f)[[1]]
  }

  new_terms
}
