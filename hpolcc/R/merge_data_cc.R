#' Merge exposure data with event data for case-crossover models
#'
#' @description
#' Build a case-crossover dataset by:
#' \enumerate{
#'   \item deriving a time-strata variable from the supplied time column,
#'   \item aggregating event counts within strata and time,
#'   \item expanding exposure rows to the strata observed in the event data, and
#'   \item joining aggregated counts back onto the exposure rows.
#' }
#'
#' @param x A `data.table` containing event rows.
#' @param exposure A `data.table` containing exposure rows.
#' @param strata_variables Character vector of stratification variables present
#'   in `x`, excluding the derived time-strata column.
#' @param time_variable Name of the date/time column shared by `x` and
#'   `exposure`.
#' @param strata_fun Function used to derive the time-strata variable from
#'   `time_variable`, defaults to a function that formats dates as "Year-Month-DayOfWeek".
#' @param time_strata_name Name of the derived time-strata column.
#' @param count_name Name of the aggregated event-count column.
#'
#' @return A `data.table` containing the expanded exposure rows and aggregated
#'   event counts. Missing event counts are filled with zero.
#'
#' @export
merge_data_cc <- function(
  x,
  exposure,
  strata_variables,
  time_variable = "date",
  strata_fun = function(x) format(x, "%Y-%b-%a"),
  time_strata_name = "timeStrata",
  count_name = "count"
) {
  stopifnot(data.table::is.data.table(x))
  stopifnot(data.table::is.data.table(exposure))
  stopifnot(time_variable %in% names(x))
  stopifnot(time_variable %in% names(exposure))

  strata_variables <- unique(c(time_strata_name, strata_variables))
  exposure_join_vars <- intersect(strata_variables, names(exposure))

  x <- data.table::copy(x)
  exposure <- data.table::copy(exposure)

  x[, (time_strata_name) := strata_fun(get(time_variable))]
  exposure[, (time_strata_name) := strata_fun(get(time_variable))]

  stopifnot(all(strata_variables %in% names(x)))
  stopifnot(all(exposure_join_vars %in% names(exposure)))
  stopifnot(all(exposure_join_vars %in% strata_variables))

  if (anyNA(x[, ..strata_variables])) {
    stop("x has missing values in stratification columns")
  }
  if (anyNA(exposure[, ..exposure_join_vars])) {
    stop("exposure has missing values in exposure join columns")
  }

  x_agg <- x[
    , setNames(list(.N), count_name),
    by = c(time_variable, strata_variables)
  ]

  keys <- unique(x[, ..strata_variables])

  exposure_all <- exposure[
    keys,
    on = exposure_join_vars,
    nomatch = 0L
  ]

  data_cc <- merge(
    exposure_all,
    x_agg,
    by = c(time_variable, strata_variables),
    all.x = TRUE
  )
  data_cc[is.na(get(count_name)), (count_name) := 0L]

  num_per_strata <- data_cc[, .N, by = strata_variables]
  if (any(num_per_strata$N == 1L)) {
    warning("There are strata with only one record in the merged data.")
  }

  data_cc
}
