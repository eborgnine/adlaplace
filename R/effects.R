#' Functions to specify terms of the formula (and associated utility functions)
#'
#' @description This suite of functions implements various moels (fixed and random effects) polynomials, integrated Wiener processes, hierarchical extensions.
#'
#' @param x A numeric vector representing the data to which the model is applied.
#' @param p A numeric value specifying the polynomial degree (default: varies by function).
#' @param ref_value A numeric value specifying the reference point for the basis functions.
#' @param knots A numeric vector specifying the knot locations for spline-based models.
#' @param range A numeric vector of length 2 specifying the range of the data (optional).
#' @param group_var A factor or grouping variable used to define hierarchical structures (optional).
#' @param include_global A logical value indicating whether to include a global component in hierarchical models (default: `TRUE`).
#' @param rpoly_p A numeric value specifying the degree for (random effects) polynomial components to add to the model as a separate term.
#' @param fpoly_p A numeric value specifying the degree for (fixed effects) polynomial components to add to the model as a separate term.
#' @param hrpoly_p A numeric value specifying the degree for hierarchical (random effects) polynomial components to add to the model as a separate term.
#' @param hfpoly_p A numeric value specifying the degree for hierarchical (fixed effects) polynomial components to add to the model as a separate term.
#' @param theta_info A list containing information about parameter indexing and initialization for hierarchical models.
#' @param term A list representing the term structure for the design matrix or precision calculations.
#' @param data A data frame containing the data to which the model terms are applied.
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A list or matrix, depending on the function, representing the design matrix, precision matrix, or parameter information for the specified model term.
#'
#' @details These functions provide polynomial and spline-based models:
#' - **Raw Polynomial Functions (rpoly, hrpoly)**: Construct design matrices for raw polynomial terms (random effects polynomial terms), with hierarchical extensions available for group-specific modeling.
#' - **Fitted Polynomial Functions (fpoly, hfpoly)**: Generate fitted polynomial terms (fixed effects polynomial terms) that adapt to specific data features, also with hierarchical extensions.
#' - **Integrated Wiener Process (iwp, hiwp)**: Develop design and precision matrices for integrated Wiener process components, with options for hierarchical group structures.
#' - **Spline-based Models (splines, knots)**: Support for spline terms, allowing for flexible non-linear modeling with specified knot locations.
#' - **Utility Functions**:
#'   - **Design Functions (`*Design`)**: Construct design matrices tailored to specific modeling terms.
#'   - **Precision Matrix Functions (`*Precision`)**: Construct precision matrices corresponding to specific modeling terms.
#'   - **Theta Functions (`*Theta`)**: Handle the variance parameters associated with random effects.
#'
#' The utility functions are specifically used within the `hm` function and are not be exported.
#'
#' @examples
#' These were generated with chatGPT and have not been reviewed. See vignette for a human generated example.
#' 
#' # Example usage for iwp
#' term <- iwp(x = 1:10, ref_value = 5, knots = c(0, 5, 10))
#' design_matrix <- iwpDesign(term, data = data.frame(x = 1:10))
#' 
#' # Example usage for hierarchical model (hiwp)
#' term <- hiwp(x = 1:10, ref_value = 5, knots = c(0, 5, 10), group_var = factor(rep(1:2, each = 5)))
#' design_matrix <- hiwpDesign(term, data = data.frame(x = 1:10, group_var = factor(rep(1:2, each = 5))))
#'
#' # Example usage for raw polynomial models
#' term <- rpoly(x = 1:10, ref_value = 5, p = 2)
#' design_matrix <- rpolyDesign(term, data = data.frame(x = 1:10))
#'
#' # Example usage for hierarchical raw polynomial (hrpoly)
#' term <- hrpoly(x = 1:10, ref_value = 5, group_var = factor(rep(1:2, each = 5)))
#' design_matrix <- hrpolyDesign(term, data = data.frame(x = 1:10, group_var = factor(rep(1:2, each = 5))))
#'
#' # Example usage for fitted polynomial models
#' term <- fpoly(x = 1:10, ref_value = 5, p = 3)
#' design_matrix <- fpolyDesign(term, data = data.frame(x = 1:10))
#'
#' # Example usage for hierarchical fitted polynomial (hfpoly)
#' term <- hfpoly(x = 1:10, ref_value = 5, group_var = factor(rep(1:2, each = 5)))
#' design_matrix <- hfpolyDesign(term, data = data.frame(x = 1:10, group_var = factor(rep(1:2, each = 5))))
#' @rdname effects_and_utilities 
NULL

.my_theta_init <- 8


#' @rdname effects_and_utilities 
#' @export
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

#' @rdname effects_and_utilities 
iwpDesign <- function(term, data){
  list2env(term, envir = environment())
  
  ref_pos <- which(knots == ref_value)
  
  # should not happen
  if(length(ref_pos) == 0) stop("ref_value of", var, "cannot be found in the corresponding knots vector. \n")
  if(range[1] < knots[1] & range[2] > rev(knots)[1]) warning("knots for ", var, " do not span its range. Continuing anyway. \n")
  
  methods::as(local_poly(knots = knots-ref_value, refined_x = data[[var]]-ref_value, p = p), "dgTMatrix")
}

#' @rdname effects_and_utilities 
iwpPrecision <- function(term){
  as(compute_weights_precision(knots=term$knots), "dgTMatrix")
}

#' @rdname effects_and_utilities 
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





#' @rdname effects_and_utilities 
#' @export
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

#' @rdname effects_and_utilities 
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

#' @rdname effects_and_utilities 
hiwpPrecision <- function(term){
  list2env(term, envir = environment())
  replicate(include_global+ngroups, iwpPrecision(term)) |> .bdiag()
}

#' @rdname effects_and_utilities 
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




#' @rdname effects_and_utilities 
#' @export
fpoly <- function(x, p = 2, ref_value = 0) {
  l <- list(var = deparse(substitute(x)), 
            type = "fpoly", 
            p = p, 
            ref_value = ref_value,
            run_as_is = F)
  return(l)
}

#' @rdname effects_and_utilities 
fpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, degree = p)
  D[,1:ncol(D),drop=F]
}




#' @rdname effects_and_utilities 
#' @export
rpoly <- function(x, p = 2, ref_value) {
  l <- list(var = deparse(substitute(x)), 
            type = "rpoly", 
            p = p, 
            ref_value = ref_value,
            run_as_is = F)
  return(l)
}

#' @rdname effects_and_utilities 
rpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, raw=T, degree = p)
  D[,1:ncol(D),drop=F]
}

#' @rdname effects_and_utilities 
rpolyPrecision <- function(term){
  as(matrix(1), "dgTMatrix")
}

#' @rdname effects_and_utilities 
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




#' @rdname effects_and_utilities 
#' @export
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

#' @rdname effects_and_utilities 
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

#' @rdname effects_and_utilities 
hrpolyPrecision <- function(term){
  list2env(term, envir = environment())
  pp <- (include_global+ngroups)*p
  sparseMatrix(i=1:pp, j=1:pp, rep = "T") |> as("dgTMatrix")
}

#' @rdname effects_and_utilities 
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







#' @rdname effects_and_utilities 
#' @export
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





#' @rdname effects_and_utilities 
#' @export
od <- function(x) {
  list(var = deparse(substitute(x)), 
       type = "od",
       run_as_is = F
  )
}

#' @rdname effects_and_utilities 
odDesign <- function(term, data){
  list2env(term, envir = environment())
  # sparseMatrix(x=1, i=1:nrow(data), j=1:nrow(data), rep = "T")[,-to_remove] |> as("dgTMatrix")
  sparseMatrix(x=1, i=1:nrow(data), j=1:nrow(data), rep = "T") |> as("dgTMatrix")
}

#' @rdname effects_and_utilities 
odPrecision <- function(term){
  list2env(term, envir = environment())
  # m <- n - length(to_remove)
  # sparseMatrix(x=1, i=1:m, j=1:m, rep = "T") |> as("dgTMatrix")
  sparseMatrix(x=1, i=1:n, j=1:n, rep = "T") |> as("dgTMatrix")
}

#' @rdname effects_and_utilities 
odTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$id), 0, max(theta_info$id))
  
  if(is.null(term$theta_id)) theta_id <- m + 1
  if(theta_id < 0) theta_id <- m + 1
  if(length(theta_id) != 1) stop("odTheta ", var, " ", type)
  
  theta_init <- .my_theta_init
  
  list(var = rep(var, length(theta_id)), type = rep(type, length(theta_id)), id = theta_id, init = theta_init)
}
