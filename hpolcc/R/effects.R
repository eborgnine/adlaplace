
get_by_levels = function(term, data) {
  if(length(term@by)) {
  if(!length(term@by_levels)) {
    unique_values <- unique(data[[term@by]])
    if(is.numeric(unique_values)) {
      unique_values_string <- 
          formatC(unique_values, width = ceiling(log10(max(unique_values))), flag = "0")
    } else {
      unique_values_string = as.character(unique_values)
    }
    term@by_levels = unique_values
    term@by_labels = unique_values_string
  }
  }
  term
}