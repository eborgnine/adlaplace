# requires timeDate and lubridate package!
removeHolidays <- function(data, type){

  if(type == "rm_none") return(data)

  if(!("year" %in% names(data))){
    cat("Computing <year> variable from data$date. Hopefully it exists\n")
    years <- year(data$date)
  }else{
    years <- data$year
  }
  year_range <- min(years):max(years)

  # First Wed not in the same week as New Year
  new_year_wday <- lubridate::wday(paste(year_range,"01","01", sep="-"))
  begin_day <- paste(year_range,"01",12 - new_year_wday, sep="-")

  # Last tuesday not in the same week of Christmas
  christmas_wday <- lubridate::wday(paste(year_range,"12","25", sep="-"))
  end_day <- paste(year_range, "12", 25 - (5:11)[christmas_wday], sep="-")

  # intervals
  keep_within <- cbind(begin_day,end_day)
  to_keep <- unlist(lapply(1:nrow(keep_within), function(i){
    which(data$date %between% keep_within[i,])
  }))

  data <- data[to_keep,]

  if(type == "rm_all"){
    # remove other holidays (need timeDate package)
    hol_to_rm <- c(timeDate::GoodFriday, 
                   timeDate::EasterSunday, 
                   timeDate::Easter, 
                   timeDate::EasterMonday, 
                   timeDate::CAVictoriaDay, 
                   timeDate::CALabourDay, 
                   timeDate::CAThanksgivingDay, 
                   timeDate::CACivicProvincialHoliday, 
                   timeDate::CACanadaDay, 
                   timeDate::CaRemembranceDay)
    holiday_index <- timeDate::holiday(min(years):max(years), hol_to_rm) |> as.Date()
    data <- data[!(data$date %in% holiday_index),]
  }

  return(data)
}


