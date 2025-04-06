#' Post-fit design matrices constructor.
#' @description Replicates the loop in `hlnm` that construct the design matrices, without modifying `terms` 
#' and creating the other objects (i.e., the precision matrices).
#'
#' @param terms The output `terms` of `hlnm`.
#' @param df A data.frame with all the necessary columns (as in the `data` input of `hlnm`).
#'
#' @return A list containing X and A, the fixed and random effects design matrices.
#'
getNewXA <- function(terms, df){
  X <- NULL
  A <- NULL
  k <- 1
  while(k <= length(terms)){
    term <- terms[[k]]
    
    if(! (term$var %in% names(df))) {
      k <- k+1
      next
    }
    
    if(term$run_as_is){
      Xsub <- sparse.model.matrix(term$f, df)
      if(is.factor(df[[term$var]])) {
        Xsub = Xsub[,-1] 
      }
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }
    
    if(term$model %in% "fpoly"){
      Xsub <- poly(df[[term$var]] - term$ref_value, raw = T, simple = T,
                   degree = term$p) |> as("TsparseMatrix")
      colnames(Xsub) <- paste0(term$var,  seq(from=1, by=1, len=ncol(Xsub)))
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }    
    
    if(term$model %in% "iid"){
      Asub <- Matrix(data = 0, nrow = nrow(df), ncol = term$n) |> as("TsparseMatrix")
      colnames(Asub) <- paste0(term$var, seq(from=1, by=1, len=ncol(Asub)))
      A <- cbind(A, Asub)
      k <- k+1
      next
    }    
    
    
    # below takes care of random effects
    # design matrix
    Asub <- getDesign(term, df)
    A <- cbind(A, Asub)

    k <- k+1
  }
  return(list(X = X, A = A))
}


#' Helper to visualize results (of an `hlnm` fit) associated with a given variable.
#' @description Replicates the loop in `hlnm` that construct the design matrices, without modifying `terms` 
#' and creating the other objects (i.e., the precision matrices).
#'
#' @param fit Output of an `hlnm` call.
#' @param exposure_var character indicating the variable of interest.
#' @param group_var character indicating the variable upon which to segment the results (should match the `group_var`
#' used in the hierarchical models of the `exposure_var`). Can be a vector, but this might lead to problems.
#' @param group character vector indicating specific group(s) (possible value(s) of the `group_var` variable) of interest.
#' @param values Numeric values at which to evaluate the `exposure_var` (set to the range spanned by the knots, if applicable).
#' @param ref_values Named list indicating the reference values for each variable that uses one. 
#' It assumes that the same reference value is used whenever a variable appears in multiple terms in the `formula`
#' input of `hlnm`.
#'
#' @return A data.frame ready to use for plotting the results.
#'
#' @export
getEffect <- function(fit, exposure_var, group_var, group, values, ref_values){
  vars <- c(
    as.character(fit$formula)[2],
    sapply(fit$terms, "[[", "var") |> unique(),
    sapply(fit$terms, "[[", "group_var") |> unlist() |> unique(),
    fit$cc_design$time_var,
    fit$cc_design$strat_var
  )
  
  group_var[is.na(group_var)] <- "__NONE__" # hopefully no one uses this...
  
  df <- data.frame(row.names = 1:(length(values)*length(group)))
  for(v in vars) df[[v]] <- ifelse(is.null(ref_values[[v]]), 0, ref_values[[v]])
  df[[group_var]] <- rep(group, each = length(values)) |> unlist()
  df[[exposure_var]] <- rep(values, times = max(1,length(group)))
  
  # # takes care of overdispersion
  # term_models <- c("", sapply(fit$terms, "[[", "model"))
  # for(v in vars[which(term_models == "iid")]) df[[v]] <- factor(1:nrow(df))
  
  list2env(getNewXA(fit$terms, df), envir = environment())
  pars <- fit$obj$env$last.par.best
  beta <- pars[names(pars) == "beta"]
  gamma <- pars[names(pars) == "gamma"]

  #
  
  data.frame(variable = exposure_var, var_value = values, 
             effect_value = as.numeric(X %*% beta + A %*% gamma), 
             group = df[[group_var]])
}



formatResult = function(obj) {
  fitList = list(
    random = list(
      hessian=try(obj$env$spHess(par = obj$env$last.par.best,random=TRUE)),
      est = obj$env$last.par.best[grep("gamma", names(obj$env$last.par.best))]),
    param = list(
      est = obj$env$last.par.best[grep("gamma", names(obj$env$last.par.best), invert=TRUE)]
    )
  )
  
  names(fitList$random$est) = colnames(obj$env$.data$A)
  if(! 'try-error' %in% class(fitList$random$hessian))     
    colnames(fitList$random$hessian) = 
      rownames(fitList$random$hessian) = colnames(obj$env$.data$A)

  names(fitList$param$est)[grep("beta", names(fitList$param$est))] = 
    colnames(obj$env$.data$X)

    fitList
  }
