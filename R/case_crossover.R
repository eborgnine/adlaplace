#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param cc_design A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param formula A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param data A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details cc_design should contain (at least some of) the following elements:
#' strat_vars: character or NULL (default). Variables with which to stratify the data.
#' time_var: character vector or NULL (default).  Variable giving the timestamps.
#' time_lag: integer (default is 7). Lag used to (further) stratify the data (only used when !is.null(time_var)).
#' time_size: integer (default is 4). Maximum size of strata (only used when !is.null(time_var)).
#' scheme: character (default is "time stratified"). Scheme used to (further) stratify the data (using lag, and only when !is.null(time_var)).
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
setStrata <- function(cc_design, data){
  
  if(is.null(cc_design$time_var) & is.null(cc_design$strat_vars)) stop("Provide statification (or time) variables.")
  if(!is.null(cc_design$time_var) & cc_design$scheme != "time stratified") stop("Only time stratified scheme is implemented...")
  # Now, the rows cc_days is run over with each possible case day within it.
  # To allow bidirectional designs, include (shadow) parameter that
  # force the TMB template to use only the first day as case day, and
  # create a row for each case day.
  
  list2env(cc_design, envir = environment())

  strat_split <- if(is.null(cc_design$strat_vars)){
    list(1:nrow(data))
  }else{
    split(1:nrow(data), interaction(data[strat_vars]), drop = T)
  }
  
  # if no more stratification do do, return right away
  if(is.null(cc_design$time_var)){
    max_len <- sapply(strat_split, length)
    cc_matrix <- matrix(unlist(strat_split), length(strat_split), max_len, byrow=T)
    
    # filter out strata of size 1
    cc_matrix <- cc_matrix[apply(cc_matrix > 0, 1, any),]
    return(cc_matrix)
  }
  
  # if not, also use time_var to further stratify 
  if(cc_design$time_lag %% 1 != 0) stop("Please provide a valid (integer) time_lag.")
  strat_split <- lapply(strat_split, \(ss){
    split(ss, as.integer(data[ss,][[time_var]]) %% time_lag, drop = F)
  }) |> unlist(recursive = F)
    
    
  if(!is.null(time_size)){
    
    strat_split <- lapply(strat_split, \(ss){
      
      times <- data[ss,][[time_var]]
      l <- length(ss)
      i <- 1
      max_lag <- time_lag*time_size
      new_ss <- list()
      while(i < l){
        maybe <- i + 1:min(time_size-1, l-i)
        new_s <- c(times[i], times[maybe][times[maybe]-times[i] < max_lag])
        if(length(new_s) != 1) new_ss[[length(new_ss) + 1]] <- new_s
        i <- i+length(new_s)
      }
      
      return(new_ss)
    }) |> unlist(recursive = F)

    # for each case day, enumerates control days (0 means empty)
    cc_matrix <- sapply(strat_split, \(ss) c(ss, rep(0, time_size-length(ss)))) |> t()
  }
  
  # filter out strata of size 1
  cc_matrix <- cc_matrix[apply(cc_matrix > 0, 1, any),]
  return(cc_matrix)
}  

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




ccDesign <- function(...){
  
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
