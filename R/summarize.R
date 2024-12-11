getNewXA <- function(terms, df){
  X <- NULL
  A <- NULL
  k <- 1
  while(k <= length(terms)){
    term <- terms[[k]]
    
    if(term$run_as_is){
      Xsub <- sparse.model.matrix(term$f, df)
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }
    
    if(term$type %in% "fpoly"){
      Xsub <- poly(df[[term$var]] - term$ref_value, raw = T, simple = T) |> as("dgTMatrix")
      X <- cbind(X, Xsub)
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
  
  list2env(getNewXA(fit$terms, df), envir = environment())
  pars <- fit$obj$env$last.par.best
  beta <- pars[names(pars) == "beta"]
  gamma <- pars[names(pars) == "gamma"]

  data.frame(variable = exposure_var, var_value = values, 
             effect_value = as.numeric(X %*% beta + A %*% gamma), 
             group = df[[group_var]])
}
