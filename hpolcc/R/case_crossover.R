#' Build a case-crossover design specification
#'
#' @param ... Named design options (\code{strat_vars}, \code{time_var},
#'   \code{time_lag}, \code{time_size}, \code{scheme}).
#' @return A list with case-crossover design settings.
#' @keywords internal
ccDesign <- function(...) {
  cc_design <- list(
    strat_vars = NULL,
    time_var = NULL,
    time_lag = 7,
    time_size = 4,
    scheme = "time stratified"
  )

  params <- list(...)
  cc_design[names(params)] <- params

  cc_design
}

#' Map observations to case-crossover strata
#'
#' @param cc_design A \code{ccDesign} list.
#' @param data Data table with stratification variables.
#' @param outcome Optional outcome variable name for stratum filtering.
#' @return Sparse matrix mapping observations to retained strata.
#' @keywords internal
setStrata <- function(cc_design, data, outcome = NULL) {
  if (is.null(cc_design$time_var) &&
    is.null(cc_design$strat_vars)) {
    stop("Provide stratification (or time) variables.")
  }
  if (!is.null(cc_design$time_var) &&
    cc_design$scheme != "time stratified") {
    stop("Only time stratified scheme is implemented...")
  }

  if (is.null(cc_design$strat_vars)) {
    data$interaction <- as.numeric(factor(data$strat_vars[1]))
  } else {
    old_cols <- setdiff(names(data), "interaction")
    data <- data[,
      c(.SD, list(interaction = .GRP)),
      by = eval(cc_design$strat_vars),
      .SDcols = old_cols
    ]
    table_interaction <- mean(table(data$interaction) <= 1)
    if (table_interaction > 0.5) {
      warning("more than half of strata only have one observation, might be missing the zeros in dataset?")
    }
  }

  if (!is.null(outcome)) {
    summaryDT <- data[, .(
      Noutcome = sum(get(outcome[1]), na.rm = TRUE),
      Ndays = .N
    ), by = interaction]
  } else {
    summaryDT <- data[, .(
      Noutcome = 1,
      Ndays = .N
    ), by = interaction]
  }
  summaryDT$toKeep <- summaryDT$Noutcome > 0 & summaryDT$Ndays > 1
  summaryDT <- summaryDT[summaryDT$toKeep, ]

  data$interaction_sub <- match(data$interaction, summaryDT$interaction)
  whichToKeep <- which(!is.na(data$interaction_sub))

  cc_matrix_sub <- cbind(i = whichToKeep, j = data$interaction_sub[whichToKeep])

  Matrix::sparseMatrix(
    i = cc_matrix_sub[, "i"], j = cc_matrix_sub[, "j"],
    dims = c(nrow(data), max(cc_matrix_sub[, "j"]))
  )
}
