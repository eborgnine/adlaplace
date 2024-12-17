#' Fit (some) hierarchical non-linear models (that we like)
#'
#' @description This function fits a hierarchical model using the specified formula and data, incorporating case-crossover designs, fixed and random effects, and precision matrices for random effects. It allows for flexible inclusion of stratification variables, time variables, and complex random effect structures.
#'
#' @param formula A formula object specifying the model to be fitted.
#' @param data A data frame containing the variables specified in the formula and any additional variables required for the model.
#' @param cc_design An object specifying the case-crossover design, including stratification and time variables. Defaults to the output of `ccDesign()`.
#' @param weight_var (Optional) A character string specifying the column in the data frame used for weights. If provided, it must exist in `data`.
#' @param tmb_parameters (Optional) A list of initial parameter values for the TMB optimization, including `beta`, `gamma`, and `theta`.
#' @param for_dev Logical; if `TRUE`, the function returns intermediate objects for development purposes. Defaults to `FALSE`.
#'
#' @return A list containing the fitted TMB object, the formula, terms used in the model, the case-crossover design, and information about the gamma and theta parameters.
#' @details The function handles fixed effects, random effects, and their associated precision matrices. It also optimizes the model using TMB with options for additional preprocessing and handling specific random effect structures.
#'
#' @import Matrix
#' @import TMB
#' 
#' @examples
#' # See vignette for basic usage
#' 
#' @useDynLib hpoltest
#' @export
#' 
hnlm <- function(formula, data, cc_design = ccDesign(), weight_var, tmb_parameters = NULL, for_dev = F) {
  
  data <- as.data.frame(data)
  
  # Check inputs
  if (!is(formula, "formula")) stop("formula must be a formula.")
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!missing(weight_var) && !is.character(weight_var) && !(weight_var %in% colnames(data)))
    stop("weight_var must be a character vector.")

  # Order the rows of data appropriately.
  if(is.null(cc_design$strat_vars) & is.null(cc_design$time_var)) stop("Provide a valid stratification (or time) variable.")
  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)
  new_order <- eval(str2lang(paste0("order(", paste0("data[['", strat_time_vars, "']]", collapse = ", "), ")")))
  data <- data[new_order,]

  
  # setup the data for case-crossover
  cc_matrix <- setStrata(cc_design = cc_design, data = data)
  
  # setup of the design matrices and other parameters
  # terms carries all the information throughout
  terms <- collectTerms(formula)
  
  # design matrices
  X <- NULL # fixed effects
  A <- NULL # random effects
  Qs <- list() # precisions for random effects
  
  beta_info <- list()
  gamma_info <- list()
  theta_info <- list()

  # loop
  k <- 1
  while(k <= length(terms)){
    term <- terms[[k]] |> getExtra(data=data, cc_matrix=cc_matrix)
    term$id <- k
    if(!is.factor(data[[term$var]][1]) && is.null(term$range))
      range <- term$range <- range(data[[term$var]])
    
    if(term$run_as_is){
      Xsub <- sparse.model.matrix(term$f, data)
      beta_info$var <- c(beta_info$var, term$var)
      beta_info$pick <- c(beta_info$pick, paste0(term$pick, "__", 0))
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }
    

    if(term$type %in% "fpoly"){
      Xsub <- poly(data[[term$var]] - term$ref_value, raw = T, simple = T) |> as("dgTMatrix")
#      colnames(Xsub) = paste0(term$var, 1:ncol(Xsub))
      beta_info$var <- c(beta_info$var, term$var)
      beta_info$pick <- c(beta_info$pick, paste0(term$pick, "__", 0))
      X <- cbind(X, Xsub)
      k <- k+1
      next
    }    
    
    # below takes care of random effects
    
    # design matrix
    Asub <- getDesign(term, data)
    A <- cbind(A, Asub)

    gamma_setup <- getGammaSetup(term)

    gamma_info$var <- c(gamma_info$var, gamma_setup$var)
    gamma_info$id <- c(gamma_info$id, gamma_setup$id)
    gamma_info$split <- c(gamma_info$split, gamma_setup$split)
    gamma_info$pick <- c(gamma_info$pick, gamma_setup$pick)
    # if(!is.null(term$group_var)) gamma_info$nrep <- c(gamma_info$nrep, length(split(1:nrow(data), data[[term$group_var]])))
    # Note: for iwp 1 knot removed for constraints
    
    # Add fized and random polynomial effects
    terms <- c(terms, addFPoly(term), addRPoly(term))
    
    
    # precision matrix
    Qs[[length(Qs) + 1]] <- getPrecision(term)
    
    # theta parameters
    theta_setup <- getThetaSetup(theta_info, term)
    
    theta_info$var <- c(theta_info$var, theta_setup$var)
    theta_info$type <- c(theta_info$type, theta_setup$type)
    theta_info$id <- c(theta_info$id, theta_setup$id)
    theta_info$init <- c(theta_info$init, theta_setup$init)
    
    # update term with new elements
    terms[[k]] <- term
    k <- k+1
  }
  
  if(for_dev) return(list(X = X, A = A, gamma_split = gamma_info$split, Qs = Qs, theta_info=theta_info, new_order = new_order))

  y <- data[[all.vars(formula)[1]]]
  tmb_data <- list(
    X = X, A = A, y = y,
    # gamma_nreplicate = gamma_info$nreplicate, # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
    Q = Qs |> .bdiag(), 
    gamma_split = gamma_info$split,
    # theta_id = theta_info$id
    cc_matrix = cc_matrix
  )
  
  if(is.null(tmb_parameters)){
    tmb_parameters <- list(
      beta = rep(0, ncol(X)),
      gamma = rep(0, ncol(A)),
      theta = theta_info$init
    )
  }

  # OPTIMIZATION ----
  # # preliminary run fixing the random effects for iwp, hiwp and od
  # # (but not the corresponding random slopes)
  # to_rm_ids <- which(theta_info$type %in% c("od", "iwp", "hiwp"))
  # 
  # if(length(to_rm_ids) > 0){
  #   gamma_split <- gamma_info$split
  #   map0 <- list(theta = factor(rep(NA, length(tmb_parameters$theta))),
  #                gamma = factor(1:length(tmb_parameters$gamma)))
  #   
  #   cs <- cumsum(c(0,gamma_split))
  #   for(id in to_rm_ids) map0$gamma[(cs[id] + 1):cs[id+1]] <- NA
  # 
  #   r <- NULL
  #   if(!all(is.na(map0$gamma))) r <- "gamma"
  #   obj0 <- MakeADFun(data = tmb_data,
  #                     parameters = tmb_parameters,
  #                     map = map0,
  #                     random = r,
  #                     DLL = "hpoltest")
  #   opt0 <- nlminb(start = obj0$par, objective = obj0$fn, gradient = obj0$gr,
  #                  control = list(eval.max=2000, iter.max=2000))
  #   tmb_parameters$beta <- opt0$par[names(opt0$par) == "beta"]
  #   if(!all(is.na(map0$gamma))) tmb_parameters$gamma[!is.na(map0$gamma)] <- obj0$env$last.par.best[names(obj0$env$last.par.best) == "gamma"]
  # }
  
  # Run optimization
  map <- list(theta = factor(theta_info$id))
  r <- NULL
  if(length(tmb_parameters$gamma) > 0) r <- "gamma"
  obj <- MakeADFun(data = tmb_data,
                   parameters = tmb_parameters,
                   random = r,
                   map = map,
                   DLL = "hpoltest")
  
  obj$fn(obj$par)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, 
                control = list(eval.max=2000, iter.max=2000))
  
  # Return the result
  return(list(obj = obj, formula = formula, 
              terms = terms, cc_design = cc_design,
              beta_info = beta_info, gamma_info = gamma_info, theta_info = theta_info,
              new_order))
}
