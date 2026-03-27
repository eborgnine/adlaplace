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
#' - `gamma_info()`: set up gamma metadata for an effect.
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
      result = list(linear(lab))
      names(result) = result[[1]]@name
      result

    }
  })
  terms_all = do.call(c, terms_1)


  terms_all
}

