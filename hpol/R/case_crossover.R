#' Stratify the observations for usage in a case-crossover model
#'
#' @description This function creates a stratification of the data based on the specified stratification variables
#' and the time-related parameters, if applicable.
#'
#' @param cc_design A list containing the design parameters for the stratification. This includes:
#'   \itemize{
#'     \item strat_vars: character or NULL (default). Variables with which to stratify the data.
#'     \item time_var: character vector or NULL (default). Variable giving the timestamps.
#'     \item time_lag: integer (default is 7). Lag used to (further) stratify the data (only used when !is.null(time_var)).
#'     \item time_size: integer (default is 4). Maximum size of strata (only used when !is.null(time_var)).
#'     \item scheme: character (default is "time stratified"). Scheme used to (further) stratify the data (using lag, and only when !is.null(time_var)).
#'   }
#' @param data A data frame containing the data to be stratified.
#'
#' @return A matrix where each row represents a strata, and each column represents a day (or unit) for that strata.
#'
#' @details The function performs stratification based on the specified `strat_vars` and `time_var`. If a time variable is provided,
#' additional stratification is done using the time lag and size parameters. The result is a matrix where each row represents a
#' group of data points corresponding to a particular stratification level, with zeroes representing empty cases.
#' Note that it handles only `time stratified` for now. If ever extended, beware that some other stratification
#' methods create overlapping strata; this may require modifying other parts of the package.
#'
#' @import data.table

setStrata <- function(cc_design, data, outcome = NULL) {
  if (is.null(cc_design$time_var) &
      is.null(cc_design$strat_vars))
    stop("Provide statification (or time) variables.")
  if (!is.null(cc_design$time_var) &
      cc_design$scheme != "time stratified")
    stop("Only time stratified scheme is implemented...")
  # Now, the rows cc_days is run over with each possible case day within it.
  # To allow bidirectional designs, include (shadow) parameter that
  # force the TMB template to use only the first day as case day, and
  # create a row for each case day.
  
  #  theEnv = list2env(cc_design, envir = environment())
  
    if (is.null(cc_design$strat_vars)) {
      data$interaction = as.numeric(factor(data$strat_vars[1]))
    } else {
    old_cols <- setdiff(names(data), "interaction")
    data <- data[,
                 c(.SD, list(interaction = .GRP)),
                 by = eval(cc_design$strat_vars),
                 .SDcols = old_cols
      ]
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
summaryDT = summaryDT[summaryDT$toKeep,]

data$interaction_sub = match(data$interaction, summaryDT$interaction)
whichToKeep = which(!is.na(data$interaction_sub))

cc_matrix_sub = cbind(i = whichToKeep, j = data$interaction_sub[whichToKeep])



    cc_sparse = Matrix::sparseMatrix(
      i = cc_matrix_sub[,'i'], j = cc_matrix_sub[,'j'], 
      dims = c(nrow(data), max(cc_matrix_sub[,'j']))
    )
  
  return(cc_sparse)
}


#' Create a list of stratification parameters
#'
#' @description This function creates a list of design parameters used in stratification. The default parameters can be
#' customized by passing arguments to the function.
#'
#' @param ... Additional arguments to customize the stratification parameters. These arguments are used to update the
#' default parameters.
#'
#' @return A list containing the stratification design parameters. The list includes:
#'   \itemize{
#'     \item strat_vars: character or NULL. Variables to stratify the data.
#'     \item time_var: character or NULL. Time variable.
#'     \item time_lag: integer. Lag used to stratify data (default is 7).
#'     \item time_size: integer. Size of strata (default is 4).
#'     \item scheme: character. Stratification scheme (default, and sole implemented scheme is "time stratified").
#'   }
#'
#' @export
ccDesign <- function(...) {
  # default
  cc_design <- list(
    strat_vars = NULL,
    time_var = NULL,
    time_lag = 7,
    time_size = 4,
    scheme = "time stratified"
  )
  
  params = list(...)
  cc_design[names(params)] <- params
  
  return(cc_design)
}

#' @export
removeUnusedStrata = function(x, Sstrata, outcome, NeventCutoff = 0) {
  setDT(x)
  
  if (!all(Sstrata %in% names(x)))
    stop("Some strata columns not found in data")
  if (!outcome[1] %in% names(x))
    stop("Outcome variable not found in data")
  
  # Calculate strata statistics using efficient data.table operations
  NinStrata <- x[, .(
    Ndays = .N,
    Nevents = sum(.SD[[1]], na.rm = TRUE)  # Explicit NA handling), by = Sstrata, .SDcols = outcome[1]]
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



# JUNK THAT MAY BE USEFUL IN GENERALIZING THE PACKAGE...
# if(design$scheme == "unidirectional"){
#   control_days <- purrr::map(-(design$n_control:1)*design$lag, ~ case_day + .x) |> Reduce(f="cbind")
#   if(design$n_control == 1) control_days <- as.matrix(control_days)
# }else if(design$scheme == "bidirectional"){
#   if(design$n_control %% 2 == 0){a <- design$n_control/2; a <- design$lag*(-a:a)[-(a+1)]}
#   else{a <- (design$n_control+1)/2; a <- (-a:a)[-c(a+1,2*a+1)]}
#   control_days <- purrr::map(a, ~ case_day + .x) |> Reduce(f="cbind")
# }else if(design$scheme == "time stratified"){
#   case_day_id <- match(case_day, time)
#   if(design$stratum_rule == "sequential"){
#     t0 <- min(time)
#     # do something with model$design$stratum_var --- data[,model$design$stratum_var]
#     # stop("error)
#     # id for the stratum (window_id, dow_id)
#     id <- paste(floor((time - t0)/(design$lag * (design$n_control+1))),
#                 (time - t0) %% design$lag, sep = "-")
#     # id <- paste(floor((time - t0)/(design$lag * (design$n_control+1))),
#     #             (time - t0) %% design$lag,
#     #             stratum_var, sep = "-")
#
#   }else if(design$stratum_rule == "month"){
#     id <- paste(format(data[, model$time_index], "%Y-%m"), time %% design$lag, sep=".")
#
#   }else stop("The stratum rule", design$stratum_rule, "is not implemented.")
#
#   # stata (case and control days togeteher)
#   stratum <- split(time, id)
