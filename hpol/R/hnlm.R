#' Fit (some) hierarchical non-linear models (that we like)
#'
#' @description This function fits a hierarchical model using the specified formula and data, incorporating case-crossover designs, fixed and random effects, and precision matrices for random effects. It allows for flexible inclusion of stratification variables, time variables, and complex random effect structures.
#'
#' @param formula A formula object specifying the model to be fitted.
#' @param data A data frame containing the variables specified in the formula and any additional variables required for the model.
#' @param cc_design An object specifying the case-crossover design, including stratification and time variables. Defaults to the output of `ccDesign()`.
#' @param dirichlet if TRUE, fit dirichlet-multinomial model, with time-within-strata gamma random effect, otherwise multinomial.
#' @param weight_var (Optional) A character string specifying the column in the data frame used for weights. If provided, it must exist in `data`.
#' @param tmb_parameters (Optional) A list of initial parameter values for the TMB optimization, including `beta`, `gamma`, and `theta`.
#' @param optim_parameters (Optional) A list of   parameters passed to `nlminb`
#' @param for_dev Logical; if `TRUE`, the function returns intermediate objects for development purposes. Defaults to `FALSE`.
#' @param verbose logical, if `TRUE` occasional information is printed
#'
#'
#' @return A list containing the fitted TMB object, the formula, terms used in the model, the case-crossover design, and information about the gamma and theta parameters.
#' @details The function handles fixed effects, random effects, and their associated precision matrices. It also optimizes the model using TMB with options for additional preprocessing and handling specific random effect structures.
#'
#' @import Matrix
#'
#' @examples
#' # See vignette for basic usage
#'
#' @useDynLib hpolcc
#' @export
hnlm <- function(formula,
                 data,
                 cc_design = ccDesign(),
                 weight_var,
                 dirichelet = FALSE,
                 tmb_parameters = NULL,
                 optim_parameters = list(eval.max = 2000, iter.max = 2000),
                 optimizer = c('nlminb', 'optim'),
                 for_dev = FALSE,
                 verbose = FALSE,
                 ...) {
  data.table::setDT(data)
  
  # Check inputs
  if (!is(formula, "formula"))
    stop("formula must be a formula.")
  #  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!missing(weight_var) &&
      !is.character(weight_var) &&
      !(weight_var %in% colnames(data)))
    stop("weight_var must be a character vector.")
  
  # Order the rows of data appropriately.
  if (is.character(cc_design)) {
    cc_design = ccDesign(strat_vars = cc_design)
  }
  if (is.null(cc_design$strat_vars) &
      is.null(cc_design$time_var))
    stop("Provide a valid stratification (or time) variable.")
  
  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)
  data.table::setorderv(data, strat_time_vars)
  
  #  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)
  #  new_order <- eval(str2lang(paste0("order(", paste0("data[['", strat_time_vars, "']]", collapse = ", "), ")")))
  #  data <- data[new_order,]
  
  
  # setup the data for case-crossover
  if (verbose) {
    cat("setting strata\n")
  }
  
  cc_matrix <- setStrata(
    cc_design = cc_design, data = data, 
    outcome = all.vars(formula)[1])

  if (verbose) cat("collecting terms\n")
  # setup of the design matrices and other parameters
  # terms carries all the information throughout
  terms <- collectTerms(formula)
  
  # design matrices
  Xlist <- list() # <- matrix(nrow=nrow(data), ncol=0) # fixed effects
  Alist <- list()#matrix(nrow=nrow(data), ncol=0) # random effects
  Qs <- list() # precisions for random effects
  
  beta_info <- list()
  gamma_info <- list()
  theta_info <- list()
  
  # loop
  k <- 1
  while (k <= length(terms)) {
    if (verbose)
      cat(k, ' ')
    
    term <- getExtra(terms[[k]], data = data, cc_matrix = cc_matrix)
    term$id <- k
    if (!is.factor(data[[term$var]][1]) &&
        !is.character(data[[term$var]][1]) &&
        is.null(term$range))
      range <- term$range <- range(data[[term$var]])
    
    if (term$run_as_is) {
      Xsub <- Matrix::sparse.model.matrix(term$f, data)
      if (is.factor(data[[term$var]])) {
        Xsub = Xsub[, -1]
      }
      beta_info$var <- c(beta_info$var, term$var)
      beta_info$pick <- c(beta_info$pick, paste0(term$pick, "__", 0))
      Xlist[[k]] <- Xsub
      k <- k + 1
      next
    }
    
    
    if (term$model %in% "fpoly") {
      Xsub <- as(
        poly(
          data[[term$var]] - term$ref_value,
          degree = term$p,
          raw = TRUE,
          simple = TRUE
        ),
        "TsparseMatrix"
      )
      colnames(Xsub) <- paste0(term$var, c('', seq(
        from = 1,
        by = 1,
        len = ncol(Xsub) - 1
      )))
      beta_info$var <- c(beta_info$var, term$var)
      beta_info$pick <- c(beta_info$pick, paste0(term$pick, "__", 0))
      Xlist[[k]] <- Xsub #cbind(X, Xsub)
      k <- k + 1
      next
    }
    # below takes care of random effects
    
    # design matrix
    Asub <- getDesign(term, data)
    Alist[[k]] <- Asub #cbind(A, Asub)
    
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
    Qs[[k]] <- getPrecision(term)
    
    # theta parameters
    theta_setup <- getThetaSetup(theta_info, term)
    
    theta_info$var <- c(theta_info$var, theta_setup$var)
    theta_info$model <- c(theta_info$model, theta_setup$model)
    theta_info$map <- c(theta_info$map, theta_setup$map)
    theta_info$init <- c(theta_info$init, theta_setup$init)
    theta_info$psd_scale_log <- c(theta_info$psd_scale_log, rep_len(c(theta_setup$psd_scale_log, 0), length(theta_info$init)))
    
    # update term with new elements
    terms[[k]] <- term
    k <- k + 1
  } # done k loop
  # final element of theta is the dirichlet SD
  theta_info$var = c(theta_info$var, 'overdisp')
  theta_info$model = c(theta_info$model, '')
  theta_info$psd_scale = exp(theta_info$psd_scale_log)

  if (dirichelet) {
    theta_info$map = c(theta_info$map, max(theta_info$map) + 1)
    dirichletStart = 0.01
    theta_info$init <- c(theta_info$init, dirichletStart)
  } 

  
  uniqueTheta = which(!duplicated(theta_info$map))
  theta_info$name = paste(theta_info$var[uniqueTheta], 
    c(theta_info$model[uniqueTheta]), 
    sep = '_')
  
  if (verbose)
    cat('.\n')
  if (length(Alist)) {
    A = do.call(cbind, Alist) |> as("TsparseMatrix")
  } else {
    A = matrix(nrow = 0, ncol = 0) |> as("TsparseMatrix")
  }
  if (length(Xlist)) {
    X = do.call(cbind, Xlist)
  } else {
    X = matrix(nrow = nrow(data), ncol = 0)
  }

  allMap = 
    rep(theta_info$map[1:length(gamma_info$split)], gamma_info$split)-1
  
  tmb_data <- list(
    X = X,
    A = A,
    y = data[[all.vars(formula)[1]]],
    # gamma_nreplicate = gamma_info$nreplicate, # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
    # Q = do.call(bdiag, Qs[!unlist(lapply(Qs, is.null))]), #Qs |> .bdiag(),
    Q =  Matrix::bdiag(Qs[!sapply(Qs, is.null)]),
    map = allMap,
    psd_scale = theta_info$psd_scale[theta_info$map],
    # theta_id = theta_info$id
    cc_matrix = cc_matrix
  )
  
  config = list(verbose = verbose,  dirichelet = as.integer(dirichelet))
  
  tmb_data = formatHpolData(tmb_data)
  
  start_parameters = c(
    rep(0, nrow(tmb_data$XTp)), # beta
    rep(0, nrow(tmb_data$ATp)), # gamma
    theta_info$init[!duplicated(theta_info$map)]
    
  )
  
  if (for_dev)
    return(
      list(
        start_parameters = start_parameters,
        theta_info = theta_info,
        tmb_data = tmb_data,
        terms = terms,
        config = config,
        data = data
      )
    )
  
  
  
  if (is.null(tmb_parameters)) {
    tmb_parameters = list()
  }
  tmbParametersDefault = list(beta = 0,
                              gamma = 0,
                              theta = theta_info$init)
  tmb_parameters = c(tmb_parameters, tmbParametersDefault[setdiff(names(tmbParametersDefault), names(tmb_parameters))])
  
  tmb_parameters$beta  = rep_len(tmb_parameters$beta, ncol(X))
  tmb_parameters$gamma = rep_len(tmb_parameters$gamma, ncol(A))
  tmb_parameters$theta = rep_len(tmb_parameters$theta, length(theta_info$init))
  
  
  objectiveFunction = function(parameters, data, 
                               config) {
                                 data = formatHpolData(data)
                                 result = objectiveFunctionC(parameters, data, config)
                                 try(result$hessian <- do.call(Matrix::sparseMatrix, resut$hessian))
                                 result
                               }
parameters = c(tmb_parameters$beta, tmb_parameters$gamma, tmb_parameters$theta)
stuff = objectiveFunction(parameters, data=tmb_data, config)


  
  if (!'map' %in% names(tmb_parameters)) {
    map <- list(theta = factor(theta_info$map))
  } else {
    map = tmb_parameters$map
  }
  
  
  optim_inline_parameters = optim_parameters[intersect(names(optim_parameters), c('upper', 'lower', 'method'))]
  theMax = apply(tmb_data$X, 2, function(xx)
    quantile(abs(xx), 0.99))
  
  Nthetas = length(unique(na.omit(theta_info$map)))
  if (is.null(optim_inline_parameters$upper))
    optim_inline_parameters$upper = c(5 / theMax, rep(5, Nthetas))
  if (is.null(optim_inline_parameters$lower))
    optim_inline_parameters$lower =  c(-5 / theMax, rep(1e-3, Nthetas))
  if (!'parscale' %in% names(optim_parameters)) {
    optim_parameters$parscale = c(theMax, rep(1 / 100, Nthetas))
  }
  
  r <- if (length(tmb_parameters$gamma) > 0) {
    "gamma"
  } else{
    NULL
  }
  
  # Return the result
  return(
    list(
      NULL
    )
  )#, funNoRandom = funNoRandom))
}
