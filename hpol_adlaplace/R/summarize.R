#' Post-fit design matrix constructor
#'
#' @description
#' Replicates the loop in `hnlm()` that constructs the design matrices,
#' without modifying `terms` or creating the associated precision
#' matrices.
#'
#' @param terms The `terms` component from an `hnlm()` fit.
#' @param df A data frame with the columns required to rebuild the design
#'   matrices.
#'
#' @return A list containing the fixed-effect matrix `X` and random-
#'   effect matrix `A`.
get_new_xa <- function(terms, df) {
  x_list <- list()
  a_list <- list()

  for (term in terms) {
    if (!(term$var %in% names(df))) {
      next
    }

    if (term$run_as_is) {
      x_sub <- Matrix::sparse.model.matrix(term$f, df)
      if (is.factor(df[[term$var]])) {
        x_sub <- x_sub[, -1, drop = FALSE]
      }
      x_list <- c(x_list, list(x_sub))
      next
    }

    if (term$model %in% "fpoly") {
      x_sub <- poly(
        df[[term$var]] - term$ref_value,
        raw = TRUE,
        simple = TRUE,
        degree = term$p
      ) |> as("TsparseMatrix")

      if (identical(term$boundary_is_random, TRUE)) {
        colnames(x_sub) <- paste0(
          term$var,
          "_fpoly_",
          seq(from = 1, len = ncol(x_sub))
        )
        a_list <- c(a_list, list(x_sub))
      } else {
        colnames(x_sub) <- paste0(
          term$var,
          seq(from = 1, by = 1, len = ncol(x_sub))
        )
        x_list <- c(x_list, list(x_sub))
      }
      next
    }

    if (term$model %in% "iid") {
      a_sub <- Matrix::Diagonal(nrow(df))
      colnames(a_sub) <- paste0("factor(", term$var, ")", df[[term$var]])
      a_list <- c(a_list, list(a_sub))
      next
    }

    a_sub <- hpolcc:::get_design(term, df)
    a_list <- c(a_list, list(a_sub))
  }

  if (length(a_list)) {
    a_matrix <- do.call(cbind, a_list)
    a_matrix <- as(a_matrix, "CsparseMatrix")
    a_matrix <- a_matrix[, which(diff(a_matrix@p) > 0), drop = FALSE]
  } else {
    a_matrix <- Matrix::Matrix(0, nrow = nrow(df), ncol = 0, sparse = TRUE)
  }

  if (length(x_list)) {
    x_matrix <- do.call(cbind, x_list)
  } else {
    x_matrix <- matrix(nrow = nrow(df), ncol = 0)
  }

  list(X = x_matrix, A = a_matrix)
}

#' Summarize effects for a fitted model
#'
#' @description
#' Construct a data frame suitable for plotting the fitted effect for a
#' given exposure variable, optionally across groups.
#'
#' @param fit Output of an `hnlm()` call.
#' @param exposure_var Character scalar naming the variable of interest.
#' @param group_var Character scalar naming the grouping variable used to
#'   segment the results. It can be `NA` to indicate no grouping.
#' @param group Character vector giving the specific group values of
#'   interest.
#' @param values Numeric vector giving the exposure values at which to
#'   evaluate the effect.
#' @param ref_values Named list of reference values used to populate the
#'   remaining variables in the design matrix construction.
#'
#' @return A data frame ready for plotting.
#' @export
get_effect <- function(
  fit,
  exposure_var,
  group_var,
  group,
  values,
  ref_values
) {
  vars <- c(
    as.character(fit$formula)[2],
    unique(sapply(fit$terms, "[[", "var")),
    unique(unlist(sapply(fit$terms, "[[", "group_var"))),
    fit$cc_design$time_var,
    fit$cc_design$strat_var
  )
  vars <- vars[!is.na(vars)]

  group_var[is.na(group_var)] <- "__NONE__"

  df <- data.frame(row.names = seq_len(length(values) * length(group)))
  for (v in vars) {
    df[[v]] <- if (is.null(ref_values[[v]])) 0 else ref_values[[v]]
  }
  df[[group_var]] <- rep(group, each = length(values)) |> unlist()
  df[[exposure_var]] <- rep(values, times = max(1, length(group)))

  xa <- get_new_xa(fit$terms, df)
  pars <- fit$obj$env$last.par.best
  beta <- pars[names(pars) == "beta"]
  gamma <- pars[names(pars) == "gamma"]

  data.frame(
    variable = exposure_var,
    var_value = values,
    effect_value = as.numeric(xa$X %*% beta + xa$A %*% gamma),
    group = df[[group_var]]
  )
}

format_result <- function(obj) {
  fit_list <- list(
    random = list(
      hessian = try(obj$env$spHess(par = obj$env$last.par.best, random = TRUE)),
      est = obj$env$last.par.best[grep("gamma", names(obj$env$last.par.best))]
    ),
    param = list(
      est = obj$env$last.par.best[
        grep("gamma", names(obj$env$last.par.best), invert = TRUE)
      ]
    )
  )

  names(fit_list$random$est) <- colnames(obj$env$.data$A)
  if (!inherits(fit_list$random$hessian, "try-error")) {
    colnames(fit_list$random$hessian) <- colnames(obj$env$.data$A)
    rownames(fit_list$random$hessian) <- colnames(obj$env$.data$A)
  }

  names(fit_list$param$est)[grep("beta", names(fit_list$param$est))] <-
    colnames(obj$env$.data$X)

  fit_list
}
