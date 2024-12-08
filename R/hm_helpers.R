# sets up the terms object that is core to the hm function.
collectTerms <- function(formula){
  term_labels <- attr(terms(formula), "term.labels")
  prefix <- NULL
  terms <- lapply(term_labels, \(lab){
    if (grepl("iwp\\(|hiwp\\(|bs\\(|rpoly\\(|hrpoly\\(|fpoly\\(|hfpoly\\(", lab)) {
      term <- eval(parse(text = lab))
      term$f <- formula(paste0("~ 0 + ", term$prefix, lab))
    }else{
      term <- list(var = lab, f = formula(paste0("~ 0 + ", lab)), run_as_is = T)
    }
    return(term)
  })
  names(terms) <- NULL
  return(terms)
}

  

# construct additional quantites specific to an effect.
getExtra <- function(term, data){
  list2env(term, envir = environment())

  if(!is.null(term$group_var)){
    term$ngroups <- if(is.factor(data[[group_var]])){
      data[[group_var]] |> droplevels() |> nlevels()
    }else{
      data[[group_var]] |> unique() |> length()
    }
  }
  
  if(!is.null(term$knots)){
    terms$nknots <- length(knots)
  }
  
  
  return(term)
}




# construct the design matrix specific to an effect.
getDesign <- function(term, data){
  list2env(term, envir = environment())
  
  Asub <- 
    if(term$type == "iwp"){
      iwpDesign(term, data)
    }else if(term$type == "hiwp"){
      hiwpDesign(term, data)
    }else if(term$type == "od"){
      odDesign(term, data)
    }else if(term$type == "rpoly"){
      rpolyDesign(term, data)
    }else if(term$type == "fpoly"){
      fpolyDesign(term, data)
    }else if(term$type == "hrpoly"){
      hrpolyDesign(term, data)
    }else if(term$type == "hfpoly"){
      hfpolyDesign(term, data)
    }else{
      stop("Unknown term type (", term$type, ")")
    }
  
  return(Asub)
}


getGammaSplit <- function(term){
  list2env(term, envir = environment())
  
    if(type == "iwp"){
      return(length(knots)-1) # IS THIS ALWAYS THAT
      
    }else if(type == "hiwp"){
      # if(include_global) return(c(length(knots)-1,ngroups*(length(knots)-1)))
      # return(ngroups*(length(knots)-1))
      return(rep(length(knots)-1, ngroups+include_global))

    }else if(type == "od"){
      return(term$n)
      
    }else if(type == "rpoly"){
      return(1)
      
    }else if(type == "hrpoly"){
      return(rep(1, (include_global + ngroups)*p))
    }

  return(NULL)
}


# construct the Penalty/Precision matrix specific to a random effect.
getPrecision <- function(term){
  list2env(term, envir = environment())
  
  Qsub <- 
    if(term$type == "iwp"){
      iwpPrecision(term)
    }else if(term$type == "hiwp"){
      hiwpPrecision(term)
    }else if(term$type == "od"){
      odPrecision(term)
    }else if(term$type == "rpoly"){
      rpolyPrecision(term)
    }else if(term$type == "hrpoly"){
      hrpolyPrecision(term)
    }else{
      stop("Unknown term type (", term$type, ")")
    }
  
  return(Qsub)
}




getThetaSetup <- function(theta_info, term){
  
  theta_id <- 
    if(term$type == "iwp"){
      iwpTheta(theta_info, term)
    }else if(term$type == "hiwp"){
      hiwpTheta(theta_info, term)
    }else if(term$type == "od"){
      odTheta(theta_info, term)
    }else if(term$type == "rpoly"){
      rpolyTheta(theta_info, term)
    }else if(term$type == "hrpoly"){
      hrpolyTheta(theta_info, term)
    }else{
      stop("Unknown term type (", term$type, ")")
    }
  return(theta_id)
}



addFPoly <- function(term){
  
  new_terms <- NULL
  
  # standard poly stuff
  if(!is.null(term$fpoly_p) && term$fpoly_p > 0){
    new_f <- "~ fpoly(__var__, ref_value = __rv__, p = __p__)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$fpoly_p) |>
      as.formula()
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  if(!is.null(term$hfpoly_p) && term$hfpoly_p > 0){
    stop("fpoly_p not yet implemented")
    # include additional fixed polynomial effects in the model
    new_f <- "~ fpoly(__var__, ref_value = __rv__, p = __p__, include_global = F)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$hfpoly_p) |>
      as.formula(env)
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  new_terms
}



addRPoly <- function(term){
  
  new_terms <- NULL
  
  if(!is.null(term$rpoly_p) && term$rpoly_p > 0){
    # include additional random polynomial effects in the model
    new_f <- "~ rpoly(__var__, ref_value = __rv__, p = __p__)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$rpoly_p) |>
      as.formula()
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  if(!is.null(term$hrpoly_p) && term$hrpoly_p > 0){
    # include additional random polynomial effects in the model
    new_f <- "~ hrpoly(__var__, ref_value = __rv__, p = __p__, group_var = __gv__, include_global = F)" |> 
      gsub(pattern = "__var__", replacement = term$var) |>
      gsub(pattern = "__rv__", replacement = term$ref_value) |>
      gsub(pattern = "__p__", replacement = term$hrpoly_p) |> 
      gsub(pattern = "__gv__", replacement = term$group_var) |>
      as.formula()
    
    new_terms[[length(new_terms)+1]] <- collectTerms(new_f)[[1]]
  }
  
  new_terms
}


