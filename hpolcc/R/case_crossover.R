setStrata <- function(cc_design, data, outcome = NULL) {
  if (is.null(cc_design$time_var) &
    is.null(cc_design$strat_vars)) {
    stop("Provide statification (or time) variables.")
  }
  if (!is.null(cc_design$time_var) &
    cc_design$scheme != "time stratified") {
    stop("Only time stratified scheme is implemented...")
  }
  # Now, the rows cc_days is run over with each possible case day within it.
  # To allow bidirectional designs, include (shadow) parameter that
  # force the TMB template to use only the first day as case day, and
  # create a row for each case day.

  #  theEnv = list2env(cc_design, envir = environment())

  if (is.null(cc_design$strat_vars)) {
    data$interaction <- as.numeric(factor(data$strat_vars[1]))
  } else {
    old_cols <- setdiff(names(data), "interaction")
    data <- data[,
      c(.SD, list(interaction = .GRP)),
      by = eval(cc_design$strat_vars),
      .SDcols = old_cols
    ]
    # check
    table_interaction <- mean(table(data$interaction) <= 1)
    if (table_interaction > 0.5) {
      warning("more than half of strata only have one observation, might be missing the zeros in dataset?")
    }
  }

  if (!is.null(outcome)) {
    # Summing by interaction for the outcome variable
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


  cc_sparse <- Matrix::sparseMatrix(
    i = cc_matrix_sub[, "i"], j = cc_matrix_sub[, "j"],
    dims = c(nrow(data), max(cc_matrix_sub[, "j"]))
  )

  return(cc_sparse)
}


ccDesign <- function(...) {
  # default
  cc_design <- list(
    strat_vars = NULL,
    time_var = NULL,
    time_lag = 7,
    time_size = 4,
    scheme = "time stratified"
  )

  params <- list(...)
  cc_design[names(params)] <- params

  return(cc_design)
}


removeUnusedStrata <- function(x, Sstrata, outcome, NeventCutoff = 0) {
  data.table::setDT(x)

  if (!all(Sstrata %in% names(x))) {
    stop("Some strata columns not found in data")
  }
  if (!outcome[1] %in% names(x)) {
    stop("Outcome variable not found in data")
  }

  # Calculate strata statistics using efficient data.table operations
  NinStrata <- x[, .(
    Ndays = .N,
    Nevents = sum(.SD[[1]], na.rm = TRUE) # Explicit NA handling), by = Sstrata, .SDcols = outcome[1]]
  ), by = Sstrata, .SDcols = outcome[1]]

  # Filter using a single pass with compound condition
  keep_strata <- NinStrata[Ndays > 1L &
    Nevents > NeventCutoff, ..Sstrata]

  # Use an anti-join pattern for better memory efficiency
  xout <- x[keep_strata, on = Sstrata, nomatch = NULL]

  # Add metadata more efficiently
  data.table::set(xout, j = "y", value = xout[[outcome[1]]])
  data.table::setattr(xout, "strata", Sstrata)
  data.table::setattr(xout, "outcome", outcome[1])


  return(xout)
}
