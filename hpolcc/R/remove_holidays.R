#' Remove Holidays from Data
#'
#' This function removes holidays from a dataset based on the specified type.
#'
#' @param data A data frame containing a 'date' column.
#' @param type A character string specifying the type of holiday removal. Options are "rm_none", "rm_all".
#' @return A data frame with holidays removed according to the specified type.
#' @export
removeHolidays <- function(data, type = "rm_all") {
  if (type == "rm_none") {
    return(data)
  }

  date_var = grep("^[Dd]ate$", names(data), value=TRUE)
  if (length(date_var) ==0 ) {
    warning("no date variable")
  }
  years <- as.integer(format(data$date, "%Y"))

  year_range <- min(years):max(years)

  # First Wed not in the same week as New Year
  new_year_wday <- as.POSIXlt(paste(year_range, "01", "01", sep = "-"))$wday
  begin_day <- paste(year_range, "01", 12 - new_year_wday, sep = "-")

  # Last Tuesday not in the same week of Christmas
  christmas_wday <- as.POSIXlt(paste(year_range, "12", "25", sep = "-"))$wday
  end_day <- paste(year_range, "12", 25 - (5:11)[christmas_wday], sep = "-")

  keep_within <- cbind(begin_day, end_day)
  to_keep <- unlist(lapply(seq_len(nrow(keep_within)), function(i) {
    which(data$date %between% keep_within[i, ])
  }))

  data <- data[to_keep, ]

  if (type == "rm_all") {
    hol_to_rm <- c(
      timeDate::GoodFriday,
      timeDate::EasterSunday,
      timeDate::Easter,
      timeDate::EasterMonday,
      timeDate::CAVictoriaDay,
      timeDate::CALabourDay,
      timeDate::CAThanksgivingDay,
      timeDate::CACivicProvincialHoliday,
      timeDate::CACanadaDay,
      timeDate::CaRemembranceDay
    )
    holiday_index <- timeDate::holiday(
      min(years):max(years),
      hol_to_rm
    ) |>
      as.Date()
    data <- data[!(data$date %in% holiday_index), ]
  }

  data
}