#' Functions to Specify Terms of the Formula (and Associated Utility Functions)
#' @rdname effects_and_utilities
#' @description This suite of functions implements various models (fixed and random effects), polynomials, integrated Wiener processes, and hierarchical extensions. 
#' The primary interface for these models is the `f()` function, which allows specification of model terms via the `model` argument. 
#'
#' @param x A numeric vector representing the data to which the model is applied.
#' @param model A character string specifying the model of model term. Options include:
#'   - `"fpoly"`: Fixed polynomial effect terms.
#'   - `"rpoly"`: Random polynomial effect terms.
#'   - `"hfpoly"`: Fixed hierarchical polynomial effect terms.
#'   - `"hrpoly"`: Random hierarchical polynomial effect terms.
#'   - `"iwp"`: Integrated Wiener process terms (random).
#'   - `"hiwp"`: Hierarchical integrated Wiener process terms (random).
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
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A list or matrix, depending on the function, representing the design matrix, precision matrix, or parameter information for the specified model term.
#'
#' @details 
#' The `f()` function serves as the primary entry point for specifying model terms. Internally, it dispatches to specific functions based on the `model` argument:
#' - `"fpoly"`: Calls the `fpoly()` function to specify fixed polynomial effect terms.
#' - `"rpoly"`: Calls the `rpoly()` function to specify random polynomial effect terms.
#' - `"hfpoly"`: Calls the `hfpoly()` function to specify hierarchical fixed polynomial effect terms.
#' - `"hrpoly"`: Calls the `hrpoly()` function to specify hierarchical random polynomial effect terms.
#' - `"iwp"`: Calls the `iwp()` function to construct integrated Wiener process terms.
#' - `"hiwp"`: Calls the `hiwp()` function for hierarchical integrated Wiener process terms.
#'
#' Utility functions:
#' - **Design Functions (`*Design`)**: Construct design matrices tailored to specific modeling terms.
#' - **Precision Matrix Functions (`*Precision`)**: Construct precision matrices corresponding to specific modeling terms.
#' - **Theta Functions (`*Theta`)**: Handle the variance parameters associated with random effects.
#'
#' These utility functions are specifically used within the `hnlm` framework and are not exported.
#'
#' @examples
#' These were generated with chatGPT and have not been reviewed. See vignette for a human generated example.
#' 
#' # Example usage with f()
#' term <- f(x = 0:10, model = "iwp", ref_value = 5, knots = seq(0,10,2))
#' term <- f(x = 0:10, model = "hiwp", ref_value = 5, knots = seq(0,10,2), group_var = "city")
#'
#' @export
f <- function(x, model = c("iwp", "hiwp", "fpoly", "rpoly", "hfpoly", "hrpoly", "iid"), ...) {
  model <- match.arg(model)
  x <- deparse(substitute(x))
  
  switch(model,
         iwp = iwp(x, ...),
         hiwp = hiwp(x, ...),
         fpoly = fpoly(x, ...),
         rpoly = rpoly(x, ...),
         hfpoly = hfpoly(x, ...),
         hrpoly = hrpoly(x, ...),
         iid = iid(x, ...),
         stop("Unknown model")
  )
}



.my_theta_init <- -2

#' @rdname effects_and_utilities 
#' @export
iwp <- function(x, p = 2, 
                ref_value, knots, range = NULL, 
                rpoly_p = 0, fpoly_p = p-1) {
  l <- list(var = x, 
            model = "iwp", 
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
  res=local_poly(knots = knots-ref_value, refined_x = data[[var]]-ref_value, p = p)
  res = methods::as(res, "TsparseMatrix")

  dimnames(res) = list(
    rownames(data), 
    paste(term$var, term$model, 1:ncol(res), sep='_')
  )
  res
}

#' @rdname effects_and_utilities 
iwpPrecision <- function(term){
  as(compute_weights_precision(knots=term$knots), "TsparseMatrix")
}

#' @rdname effects_and_utilities 
iwpTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$map), 0, max(theta_info$map))
  
  if(is.null(term$theta_map)) theta_map <- m + 1
  if(theta_map < 0) theta_map <- m + 1
  if(length(theta_map) != 1) stop("iwpTheta ", var, " ", model)
  
  # if(is.null(theta_init)) theta_init <- .my_theta_init
  # if(length(theta_init) != 1) stop("iwpTheta ", var, " ", model)
  theta_init <- .my_theta_init
  
  # predictive SD
  iqrKnots = diff(quantile(term$knots, c(0.25, 0.75)))
  pMhalf = term$p-1/2
  logPsd = pMhalf * log(iqrKnots)- 0.5*log(2*pMhalf)-lfactorial(term$p-1)
  names(logPsd) = var
  
  list(var = var, model = model,
       name = paste0(var, "_", model),
       name_mapped = paste0(var, "_", model),
       level_mapped = "GLOBAL",
       psd_scale_log = logPsd,
       map = theta_map, init = theta_init)
}





#' @rdname effects_and_utilities 
#' @export
hiwp <- function(x, p = 2, ref_value, knots, range = NULL, 
                 group_var,
                 include_global = T,
                 hrpoly_p = p-1, hfpoly_p = 0, # with include_global = F
                 rpoly_p = 0, fpoly_p = include_global*(p-1)) {
  l <- list(var = x, 
            model = "hiwp", 
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
hiwpDesign <- function(term, data, use_dev_version = F){
  list2env(term, envir = environment())
  
  A0 <- iwpDesign(term, data)
  id_split <- split(1:nrow(data), 
                    factor(data[[term$group_var]], levels = term$groups), 
                    drop = F)
  if(term$include_global) id_split <- c(list(global=1:nrow(data)), id_split)
  
  if(TRUE){
    # fixed
    A0split = mapply(function(AA, xx) {
      res = as(AA[xx, ], "TsparseMatrix")
      cbind(i=xx[1+res@i], j=res@j+1, x=res@x)
    }, 
    xx = id_split, MoreArgs = list(AA=A0), SIMPLIFY=FALSE)
    A0combine = cbind(
      as.data.frame(do.call(rbind, A0split)), 
      split = rep(1:length(A0split), unlist(lapply(A0split, nrow)))
    )
    A0combine[,'j2'] = A0combine[,'j'] + ncol(A0) * (A0combine[,'split']-1)
    Afinal = Matrix::sparseMatrix(i=A0combine$i, j=A0combine$j2, x=A0combine$x,
      dims = c(nrow(data), ncol(A0)*length(id_split)),
      dimnames = list(rownames(data), 
        paste(term$var, term$model, 
          rep(names(id_split), each=ncol(A0)),
          rep(1:ncol(A0), length(id_split)), sep='_')))
  }else{
    
    # the slow way
    Afinal <- Matrix::Matrix(0, nrow=nrow(data), ncol=ncol(A0)*length(id_split)) |> as("TsparseMatrix")
    for(k in seq_along(id_split)) 
      Afinal[id_split[[k]], (k-1)*ncol(A0) + 1:ncol(A0)] <- A0[id_split[[k]],]
    AfinalOld = Afinal
  }
  
  Afinal
}

#' @rdname effects_and_utilities 
hiwpPrecision <- function(term){
  list2env(term, envir = environment())
  Matrix::.bdiag(replicate(include_global+ngroups, iwpPrecision(term)))
}

#' @rdname effects_and_utilities 
hiwpTheta <- function(theta_info, term){
  m <- ifelse(is.null(theta_info$map), 0, max(theta_info$map))
  list2env(term, envir = environment())

  if(is.null(term$theta_map)){
    # technically (p-1) * include_global below right? No, just for slope components.
    theta_map <- m + c(include_global + rep(1, ngroups))
    if(include_global) theta_map <- c(m + 1, theta_map)
  } 
  if(length(theta_map) != include_global+ngroups) stop("hiwpTheta ", var, " ", model)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", model)
  theta_init <- rep(.my_theta_init, include_global+ngroups)
  
  name <- paste(rep(var, ngroups), rep("iwp", ngroups), as.character(term$groups), sep="_")
  if(include_global) name <- c(paste(var, "iwp", "GLOBAL", sep="_"), name)
  
  name_mapped <- paste(var, "iwp", "LOCAL", sep="_")
  level_mapped <- "LOCAL"
  if(include_global){
    name_mapped <- c(paste(var, "iwp", "GLOBAL", sep="_"), name_mapped)
    level_mapped <- c("GLOBAL", level_mapped)
  } 
  
  list(var = rep(var, length(theta_map)), model = rep(model, length(theta_map)), 
       name = name,
       name_mapped = name_mapped,
       level_mapped = level_mapped,
       map = theta_map, init = theta_init)
}




#' @rdname effects_and_utilities 
#' @export
fpoly <- function(x, p = 2, ref_value = 0) {
  l <- list(var = x, 
            model = "fpoly", 
            p = p, 
            ref_value = ref_value,
            run_as_is = F)
  return(l)
}

#' @rdname effects_and_utilities 
fpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, degree = p)
  D <- D[,1:ncol(D),drop=F]
  colnames(D) <- paste0(term$var, c('', seq(from=1, by=1, len=ncol(D)-1)))
  
  D
}




#' @rdname effects_and_utilities 
#' @export
rpoly <- function(x, p = 2, ref_value) {
  l <- list(var = x, 
            model = "rpoly", 
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
  as(matrix(1), "TsparseMatrix")
}

#' @rdname effects_and_utilities 
rpolyTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$map), 0, max(theta_info$map))
  
  if(is.null(term$theta_map)) theta_map <- m + 1:p
  if(length(theta_map) != p) stop("rpolyTheta ", var, " ", model)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", model)
  theta_init <- rep(.my_theta_init, p)
  
  list(var = rep(var, length(theta_map)), model = rep(model, length(theta_map)), 
       name = paste0(var, "_", model, "_", 1:p),
       name_mapped = paste0(var, "_", model, "_", 1:p),
       level_mapped = "GLOBAL",
       map = theta_map, init = theta_init)
}



#' @rdname effects_and_utilities 
#' @export
hrpoly <- function(x, p = 1, ref_value, group_var, include_global = T) {
  l <- list(var = x, 
            model = "hrpoly", 
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
  Afinal <- Matrix::Matrix(0, nrow=nrow(data), ncol=pp) |> as("TsparseMatrix")
  if(ig) Afinal[,1:p] <- A0
  for(k in seq_along(id_split)) 
    Afinal[id_split[[k]], (ig+k-1)*p + 1:p] <- A0[id_split[[k]],]
  colnames(Afinal) = paste(
    term$var, term$model, 
    rep(term$groups, each=term$p), rep(1:term$p, length(term$groups)), sep='_')

  Afinal
}

#' @rdname effects_and_utilities 
hrpolyPrecision <- function(term){
  list2env(term, envir = environment())
  pp <- (include_global+ngroups)*p
  Matrix::sparseMatrix(i=1:pp, j=1:pp, rep = "T") |> as("TsparseMatrix")
}

#' @rdname effects_and_utilities 
hrpolyTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  
  m <- ifelse(is.null(theta_info$map), 0, max(theta_info$map))
  if(is.null(term$theta_map)){
    theta_map <- m + c(p*include_global + rep(1:p, ngroups))
    if(include_global) theta_map <- c(m + 1:p, theta_map)
  } 
  if(length(theta_map) != (include_global+ngroups)*p) stop("hrpolyTheta ", var, " ", model)
  
  # if(is.null(theta_init)) theta_init <- rep(.my_theta_init, p)
  # if(length(theta_init) != p) stop("rpolyTheta ", var, " ", model)
  theta_init <- rep(.my_theta_init, (include_global+ngroups)*p)
  
  name <- paste(rep(var, ngroups), rep("rpoly", ngroups*p), rep(as.character(term$groups), each=p), rep(1:p, times=ngroups), sep="_")
  if(include_global) name <- c(paste(var, "rpoly", "GLOBAL", 1:p, sep="_"), name)

  name_mapped <- paste(var, "rpoly", "LOCAL", 1:p, sep="_")
  level_mapped <- "LOCAL"
  if(include_global){
    name_mapped <- c(paste(var, "rpoly", "GLOBAL", 1:p, sep="_"), name_mapped)
    level_mapped <- c("GLOBAL", level_mapped)
  } 
  
  list(var = rep(var, length(theta_map)), model = rep(model, length(theta_map)), 
       name = name, name_mapped = name_mapped,
       level_mapped = level_mapped,
       map = theta_map, init = theta_init)
}







# #' @rdname effects_and_utilities 
# #' @export
# bs <- function(x, o = 2, ref_value, range = NULL, nknots = NULL, run_as_is = T, ...) {
#   
#   stop("bs() is not implemented.")
#   
#   l <- c(list(var = x, 
#               model = "bs", 
#               o = o, 
#               ref_value = ref_value,
#               nknots = nknots,
#               range = range,
#               run_as_is = run_as_is,
#               prefix = "splines::"),
#          list(...))
#   
#   if(!is.null(range)){
#     if(ref_value < range[1] | ref_value > range[2])
#       stop("ref_value not in the range of the data for", var, " (bs).")
#   }
#   
#   return(l)
# }





#' @rdname effects_and_utilities 
#' @export
iid <- function(x) {
  list(var = x, 
       model = "iid",
       run_as_is = F
  )
}

#' @rdname effects_and_utilities 
iidDesign <- function(term, data){
  list2env(term, envir = environment())
  ff <- paste0("~ 0 + factor(", var, ")") |> formula()
  Matrix::sparse.model.matrix(ff, data) |> as("TsparseMatrix")
}

#' @rdname effects_and_utilities 
iidPrecision <- function(term){
  list2env(term, envir = environment())
  Matrix::sparseMatrix(x=1, i=1:n, j=1:n, rep = "T") |> as("TsparseMatrix")
}

#' @rdname effects_and_utilities 
iidTheta <- function(theta_info, term){
  list2env(term, envir = environment())
  m <- ifelse(is.null(theta_info$map), 0, max(theta_info$map))
  
  if(is.null(term$theta_map)) theta_map <- m + 1
  if(theta_map < 0) theta_map <- m + 1
  if(length(theta_map) != 1) stop("odTheta ", var, " ", model)
  
  theta_init <- .my_theta_init
  
  list(var = rep(var, length(theta_map)), model = rep(model, length(theta_map)), 
       name = paste(var, model, sep="_"), 
       name_mapped = paste(var, model, sep="_"),
       level_mapped = "GLOBAL",
       map = theta_map, init = theta_init)
}
