#' Helpers for the `hnlm` function
#' @name hnlm_helpers
#' @rdname hnlm_helpers
#' 
#' @description This page provides common information for functions that handle the setup and 
#' processing of terms used in the model.
#' These functions include:
#' - `collectTerms()`: Sets up the terms object, which is the core of the model.
#' - `getExtra()`: Constructs additional quantities specific to an effect.
#' - `getDesign()`: Constructs the design matrix for a given effect.
#' - `getGammaSetup()`: Sets up the gamma parameters specific to an effect.
#' - `getPrecision()`: Constructs the precision matrix for a random effect.
#' - `getThetaSetup()`: Sets up theta information for a specific effect.
#' - `addFPoly()`: Adds fixed polynomial effects to the model.
#' - `addRPoly()`: Adds random polynomial effects to the model.
#' 
NULL



# sets up the terms object that is core to the hnlm function.
#' @rdname hnlm_helpers
collectTerms <- function(formula){
  term_labels <- attr(terms(formula), "term.labels")
  prefix <- NULL
  
  terms <- lapply(term_labels, function(lab) {
    # Check if the term is defined using `f()`
    if (grepl("f\\(", lab)) {
      term <- eval(parse(text = lab))
      if (!is.null(term$model)) {
        term$f <- formula(paste0("~ 0 + ", term$prefix, lab))
      } else {
        stop("The term defined with `f()` must specify a model.")
      }
      
    } else {
      term <- list(var = lab, model = "", f = formula(paste0("~ 0 + ", lab)), run_as_is = TRUE)
    }
    
    return(term)
  })
  
  names(terms) <- NULL
  return(terms)
}

  

# construct additional quantites specific to an effect.
#' @rdname hnlm_helpers
getExtra <- function(term, data, cc_matrix){
  list2env(term, envir = environment())
  if(term$run_as_is) return(term)
  
  if(term$model == "iid"){
    term$n <- length(unique(data[[term$var]]))
    # term$to_remove <- cc_matrix[,1] # IF WE GENERALIZE TO SOME OVERLAPPING STRATA, watch out here.
  } 
    
  
  if(!is.null(term$group_var)){
    term$groups <- factor(data[[group_var]] |> unique()) |> droplevels()
    if(any(term$groups == "GLOBAL")) stop("Please change value (name) of GLOBAL in the group variable for ", term$var, " -- ", term$model)
    term$groups <- factor(term$groups, levels = c("GLOBAL", levels(term$groups)))
    
    term$ngroups <- term$groups |> length() 
  }
  
  if(!is.null(term$knots)){
    term$nknots <- length(knots)
  }
  
  
  return(term)
}




# construct the design matrix specific to an effect.
#' @rdname hnlm_helpers
getDesign <- function(term, data){
  list2env(term, envir = environment())
  
  Asub <- 
    if(term$model == "iwp"){
      iwpDesign(term, data)
    }else if(term$model == "hiwp"){
      hiwpDesign(term, data)
    }else if(term$model == "iid"){
      iidDesign(term, data)
    }else if(term$model == "rpoly"){
      rpolyDesign(term, data)
    }else if(term$model == "fpoly"){
      fpolyDesign(term, data)
    }else if(term$model == "hrpoly"){
      hrpolyDesign(term, data)
    }else if(term$model == "hfpoly"){
      hfpolyDesign(term, data)
    }else{
      stop("Unknown term model (", term$model, ")")
    }
  
  return(Asub)
}


#' @rdname hnlm_helpers
getGammaSetup <- function(term){
  list2env(term, envir = environment())
  
  gamma_split <- if(model == "iwp"){
      length(knots)-1 # IS THIS ALWAYS THAT
      
  }else if(model == "hiwp"){
    # if(include_global) return(c(length(knots)-1,ngroups*(length(knots)-1)))
    # return(ngroups*(length(knots)-1))
    rep(length(knots)-1, ngroups+include_global)

  }else if(model == "iid"){
    # return(term$n - length(term$to_remove))
    term$n
    
  }else if(model == "rpoly"){
    rep(1,p)
    
  }else if(model == "hrpoly"){
    rep(1, (include_global + ngroups)*p)
  }

  gamma_var <- rep(var, sum(gamma_split))
  gamma_id <- rep(id, sum(gamma_split))

  gamma_pick <- if(model == "iwp"){
    factor("GLOBAL") # IS THIS ALWAYS THAT
    
  }else if(model == "hiwp"){
    gp <- lapply(seq_along(groups), \(i) rep(as.integer(groups)[i]-1, gamma_split[i+include_global])) |> unlist()
    if(include_global) gp <- c(rep(0, gamma_split[1]), gp)
    paste0(var, "__", gp)
    
  }else if(model == "iid"){
    paste0(var, "__", 0)
    
  }else if(model == "rpoly"){
    paste0(var, "__", 0)
    
  }else if(model == "hrpoly"){
    gp <- lapply(seq_along(groups), \(i) rep(as.integer(groups)[i]-1, gamma_split[i+include_global])) |> unlist()
    if(include_global) gp <- c(rep(0, gamma_split[1]), gp)
    paste0(var, "__", gp)
  }
  
  return(list(var = gamma_var, id = gamma_id, split = gamma_split, pick = gamma_pick))
}


# construct the Penalty/Precision matrix specific to a random effect.
#' @rdname hnlm_helpers
getPrecision <- function(term){
  list2env(term, envir = environment())
  
  Qsub <- 
    if(term$model == "iwp"){
      iwpPrecision(term)
    }else if(term$model == "hiwp"){
      hiwpPrecision(term)
    }else if(term$model == "iid"){
      iidPrecision(term)
    }else if(term$model == "rpoly"){
      rpolyPrecision(term)
    }else if(term$model == "hrpoly"){
      hrpolyPrecision(term)
    }else{
      stop("Unknown term model (", term$model, ")")
    }
  
  return(Qsub)
}


#' @rdname hnlm_helpers
getThetaSetup <- function(theta_info, term){
  
  theta_setup <- 
    if(term$model == "iwp"){
      iwpTheta(theta_info, term)
    }else if(term$model == "hiwp"){
      hiwpTheta(theta_info, term)
    }else if(term$model == "iid"){
      iidTheta(theta_info, term)
    }else if(term$model == "rpoly"){
      rpolyTheta(theta_info, term)
    }else if(term$model == "hrpoly"){
      hrpolyTheta(theta_info, term)
    }else{
      stop("Unknown term model (", term$model, ")")
    }
  return(theta_setup)
}



#' @rdname hnlm_helpers
addFPoly <- function(term){
  
  new_terms <- NULL
  
  # standard poly stuff
  if(!is.null(term$fpoly_p) && term$fpoly_p > 0){
    new_f <- "~ f(__var__, model = 'fpoly', ref_value = __rv__, p = __p__)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$fpoly_p) |>
      as.formula()
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  if(!is.null(term$hfpoly_p) && term$hfpoly_p > 0){
    stop("hfpoly_p not yet implemented")
    # include additional fixed polynomial effects in the model
    new_f <- "~ f(__var__, model = 'hfpoly', ref_value = __rv__, p = __p__, include_global = F)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$hfpoly_p) |>
      as.formula(env)
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  new_terms
}



#' @rdname hnlm_helpers
addRPoly <- function(term){
  
  new_terms <- NULL
  
  if(!is.null(term$rpoly_p) && term$rpoly_p > 0){
    # include additional random polynomial effects in the model
    new_f <- "~ f(__var__, model = 'rpoly', ref_value = __rv__, p = __p__)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$rpoly_p) |>
      as.formula()
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  if(!is.null(term$hrpoly_p) && term$hrpoly_p > 0){
    # include additional random polynomial effects in the model
    new_f <- "~ f(__var__, model = 'hrpoly', ref_value = __rv__, p = __p__, group_var = __gv__, include_global = F)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$hrpoly_p) |> 
      gsub(pattern = "__gv__", replacement = term$group_var) |>
      as.formula()
    
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  new_terms
}


