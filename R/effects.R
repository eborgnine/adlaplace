.my_theta_init <- 8

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
                rpoly_p = 0, fpoly_p = 1) {
  l <- list(var = deparse(substitute(x)), 
            type = "iwp", 
            p = p, 
            ref_value = ref_value,
            knots = knots,
            range = range,
            rpoly_p = rpoly_p, 
            fpoly_p = fpoly_p,
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
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  
  if(is.null(term$theta_id)) theta_id <- m + 1
  if(theta_id < 0) theta_id <- m + 1
  if(length(theta_id) != 1) stop("iwpTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- .my_theta_init
  # if(length(theta_init) != 1) stop("iwpTheta ", var, " ", type)
  theta_init <- .my_theta_init
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)),
       id = theta_id, init = theta_init)
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
hiwp <- function(x, p = 2, ref_value, knots, range = NULL, 
                 group_var,
                 include_global = T,
                 hrpoly_p = p-1, hfpoly_p = 0, # with include_global = F
                 rpoly_p = 0, fpoly_p = include_global) {
  l <- list(var = deparse(substitute(x)), 
            type = "hiwp", 
            p = p,
            ref_value = ref_value,
            knots = knots, 
            range = range,
            group_var = deparse(substitute(group_var)),
            include_global = include_global,
            hrpoly_p = hrpoly_p, 
            hfpoly_p = hfpoly_p,
            rpoly_p = rpoly_p, 
            fpoly_p = fpoly_p,
            run_as_is = F)
  
  if(!is.null(range)){
    if(ref_value < range[1] | ref_value > range[2])
      stop("ref_value not in the range of the data for", var, " (bs).")
  }
  
  return(l)
}

hiwpDesign <- function(term, data){
  list2env(term, envir = environment())
  
  A0 <- iwpDesign(term, data)
  id_split <- split(1:nrow(data), factor(data[[group_var]], levels = groups), drop = F)
  if(include_global) id_split <- c(list(1:nrow(data)), id_split)
  
  Afinal <- Matrix(0, nrow=nrow(data), ncol=ncol(A0)*length(id_split)) |> as("dgTMatrix")
  for(k in seq_along(id_split)){
    Afinal[id_split[[k]], (k-1)*ncol(A0) + 1:ncol(A0)] <- A0[id_split[[k]],]
  } 
    
  Afinal
}

hiwpPrecision <- function(term){
  list2env(term, envir = environment())
  replicate(include_global+ngroups, iwpPrecision(term)) |> .bdiag()
}

hiwpTheta <- function(theta_info, term){
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  list2env(term, envir = environment())
  ig <- include_global

  if(is.null(term$theta_id)){
    theta_id <- m + c(include_global + rep(1, ngroups))
    if(include_global) theta_id <- c(m + 1, theta_id)
  } 
  if(length(theta_id) != include_global+ngroups) stop("hiwpTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", type)
  theta_init <- rep(.my_theta_init, include_global+ngroups)
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)), id = theta_id, init = theta_init)
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
fpoly <- function(x, p = 2, ref_value = 0) {
  l <- list(var = deparse(substitute(x)), 
            type = "fpoly", 
            p = p, 
            ref_value = ref_value,
            run_as_is = F)
  return(l)
}

fpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, degree = p)
  D[,1:ncol(D),drop=F]
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
rpoly <- function(x, p = 2, ref_value) {
  l <- list(var = deparse(substitute(x)), 
            type = "rpoly", 
            p = p, 
            ref_value = ref_value,
            run_as_is = F)
  return(l)
}

rpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, raw=T, degree = p)
  D[,1:ncol(D),drop=F]
}

rpolyPrecision <- function(term){
  as(matrix(1), "dgTMatrix")
}

rpolyTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  
  if(is.null(term$theta_id)) theta_id <- m + 1:p
  if(length(theta_id) != p) stop("rpolyTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", type)
  theta_init <- rep(.my_theta_init, p)
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)), id = theta_id, init = theta_init)
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
hrpoly <- function(x, p = 1, ref_value, group_var, include_global = T) {
  l <- list(var = deparse(substitute(x)), 
            type = "hrpoly", 
            p = p, 
            ref_value = ref_value,
            group_var = deparse(substitute(group_var)),
            include_global = include_global,
            run_as_is = F)
  return(l)
}

hrpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  ig <- include_global
  
  id_split <- split(1:nrow(data), factor(data[[group_var]], levels = groups), drop = F)
  pp <- (ig + length(id_split)) * p
  mm <- c(0, sapply(id_split, length)) |> cumsum()
  
  A0 <- rpolyDesign(term, data)
  Afinal <- Matrix(0, nrow=nrow(data), ncol=pp) |> as("dgTMatrix")
  if(ig) Afinal[,1:p] <- A0
  for(k in seq_along(id_split)) 
    Afinal[id_split[[k]], (ig+k-1)*p + 1:p] <- A0[id_split[[k]],]
  
  Afinal
}

hrpolyPrecision <- function(term){
  list2env(term, envir = environment())
  pp <- (include_global+ngroups)*p
  sparseMatrix(i=1:pp, j=1:pp, rep = "T") |> as("dgTMatrix")
}

hrpolyTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  if(is.null(term$theta_id)){
    theta_id <- m + c(p*include_global + rep(1:p, ngroups))
    if(include_global) theta_id <- c(m + 1:p, theta_id)
  } 
  if(length(theta_id) != (include_global+ngroups)*p) stop("hrpolyTheta ", var, " ", type)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", type)
  theta_init <- rep(.my_theta_init, (include_global+ngroups)*p)
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)), id = theta_id, init = theta_init)
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
bs <- function(x, o = 2, ref_value, range = NULL, nknots = NULL, run_as_is = T, ...) {
  
  stop("bs() is not implemented.")
  
  l <- c(list(var = deparse(substitute(x)), 
              type = "bs", 
              o = o, 
              ref_value = ref_value,
              nknots = nknots,
              range = range,
              run_as_is = run_as_is,
              prefix = "splines::"),
         list(...))
  
  if(!is.null(range)){
    if(ref_value < range[1] | ref_value > range[2])
      stop("ref_value not in the range of the data for", var, " (bs).")
  }
  
  return(l)
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
od <- function(x) {
  list(var = deparse(substitute(x)), 
       type = "od",
       run_as_is = F
  )
}

odDesign <- function(term, data){
  list2env(term, envir = environment())
  # sparseMatrix(x=1, i=1:nrow(data), j=1:nrow(data), rep = "T")[,-to_remove] |> as("dgTMatrix")
  sparseMatrix(x=1, i=1:nrow(data), j=1:nrow(data), rep = "T") |> as("dgTMatrix")
}

odPrecision <- function(term){
  list2env(term, envir = environment())
  # m <- n - length(to_remove)
  # sparseMatrix(x=1, i=1:m, j=1:m, rep = "T") |> as("dgTMatrix")
  sparseMatrix(x=1, i=1:n, j=1:n, rep = "T") |> as("dgTMatrix")
}

odTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  
  if(is.null(term$theta_id)) theta_id <- m + 1
  if(theta_id < 0) theta_id <- m + 1
  if(length(theta_id) != 1) stop("odTheta ", var, " ", type)
  
  theta_init <- .my_theta_init
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)), id = theta_id, init = theta_init)
}
