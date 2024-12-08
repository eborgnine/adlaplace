#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param x A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param y A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details Provide any additional details about the function, such as edge cases, assumptions, or implementation notes.
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
iwp <- function(x, p = 2, 
                ref_value, knots, range = NULL, 
                rpoly_p = p-1, fpoly_p = 0,
                theta_hyper = .01) {
  l <- list(var = deparse(substitute(x)), 
             type = "iwp", 
             p = p, 
             ref_value = ref_value,
             knots = knots,
            range = range,
            rpoly_p = rpoly_p, 
            fpoly_p = fpoly_p,
            theta_hyper = theta_hyper,
            run_as_is = F)
  
  if(!is.null(range)){
    if(ref_value < range[1] | ref_value > range[2])
      stop("ref_value not in the range of the data for", var, " (bs).")
  }

  return(l)
}

iwpDesign <- function(term, data){
  list2env(term, envir = environment())
  
  ref_pos <- which(knots == ref_value)
  
  # should not happen
  if(length(ref_pos) == 0) stop("ref_value of", var, "cannot be found in the corresponding knots vector. \n")
  if(range[1] < knots[1] & range[2] > rev(knots)[1]) warning("knots for ", var, " do not span its range. Continuing anyway. \n")
  
  methods::as(local_poly(knots = knots-ref_value, refined_x = data[[var]]-ref_value, p = p), "dgTMatrix")
}

iwpPrecision <- function(term){
  as(compute_weights_precision(knots=term$knots), "dgTMatrix")
}

iwpTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  
  if(is.null(term$theta_id)) theta_id <- 0L
  if(theta_id < 0) theta_id <- 0L
  if(length(theta_id) != 1) stop("iwpTheta ", var, " ", type)
  
  if(is.null(term$theta_hyper)) theta_hyper <- .001
  if(length(theta_hyper) != 1) stop("iwpTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- 9
  # if(length(theta_init) != 1) stop("iwpTheta ", var, " ", type)
  theta_init <- 9

  list(id = theta_id, hyper = theta_hyper, init = theta_init)
}


#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param x A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param y A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details Provide any additional details about the function, such as edge cases, assumptions, or implementation notes.
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
hiwp <- function(x, p = 2, ref_value, knots, range, 
                 group_var, 
                 rpoly_p = p-1, fpoly_p = 0,
                 theta_hyper = .01) {
  l <- list(var = deparse(substitute(x)), 
       type = "hiwp", 
       p = p,
       ref_value = ref_value,
       knots = knots, 
       range = range,
       group_var = deparse(substitute(group_var)),
       rpoly_p = rpoly_p, 
       fpoly_p = fpoly_p,
       theta_hyper = theta_hyper,
       run_as_is = F)
  
  if(!is.null(range)){
    if(ref_value < range[1] | ref_value > range[2])
      stop("ref_value not in the range of the data for", var, " (bs).")
  }

  return(l)
}

hiwpDesign <- function(){
  
}

hiwpPrecision <- function(){
  
}

hiwpTheta <- function(theta_info, term){
  stop("hiwpTheta TO DO")
  
  list2env(term, envir = environment())
  
  lg <- length(group_levels)
  lt <- length(theta_id)
  
  if(is.null(theta_id)) return(rep(0, l+1))
  if(!(lt %in% c(1,2,l+1))) stop("theta_id must be NULL or of length 1, 2, or same as group_levels (in hiwp for ", var, ")")
  
  if(any(theta_id > 0) & any(theta_id < 0)) stop("mixed theta_id (<0 and >0) for hiwp ", var, ". Not yet implemented.")
  if(lt == l+1) return(theta_id)
  
  if(lt == 1){
    m <- ifelse(theta_id < 0, min(theta_info$id) - 1, theta_id)
    return(rep(m, l+1))
  }else{
    if(theta[1] == theta[2]) stop("length(theta)=2 but they are equal(?)")
    if(theta_id[1] < 0){
      m <- min(theta_info$id) - 1:2
      return(c(m[1], rep(m[2], l)))
    }else{
      return(c(theta[1], rep(theta[2], l)))
    }
  }
}




#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param x A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param y A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details Provide any additional details about the function, such as edge cases, assumptions, or implementation notes.
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
rpoly <- function(x, p = 2, ref_value, theta_hyper = .01) {
  l <- list(var = deparse(substitute(x)), 
            type = "rpoly", 
            p = p, 
            ref_value = ref_value,
            theta_hyper = theta_hyper,
            run_as_is = F)
  return(l)
}

rpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, degree = p)
  D[,1:ncol(D),drop=F]
}

rpolyPrecision <- function(term){
  as(matrix(1), "dgTMatrix")
}

rpolyTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  
  if(is.null(term$theta_id)) theta_id <- rep(0,p)
  if(length(theta_id) != p) stop("rpolyTheta ", var, " ", type)
  
  if(is.null(theta_hyper)) theta_hyper <- rep(.001, p)
  if(length(theta_hyper) != p) stop("rpolyTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- rep(9, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", type)
  theta_init <- rep(9, p)

  list(id = theta_id, hyper = theta_hyper, init = theta_init)
}
