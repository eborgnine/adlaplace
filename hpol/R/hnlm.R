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
 dirichlet = FALSE,
 tmb_parameters = NULL,
 optim_parameters = list(eval.max = 2000, iter.max = 2000),
 optimizer = c('nlminb', 'optim'),
 config=list(),
 control=list(),
 control_inner = control,
 for_dev = FALSE,
 verbose = FALSE,
 ...) {

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
  

  data.table::setDT(data)
  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)

  strat_time_vars = strat_time_vars[
    order(sapply(data[, strat_time_vars, with = FALSE], function(xx) length(unique(xx))), decreasing=FALSE)
  ]


  data.table::setorderv(data, strat_time_vars)
  
  # setup the data for case-crossover
  if (verbose) {
    cat("setting strata\n")
  }
  
  cc_matrix <- setStrata(
    cc_design = cc_design, 
    data = data, 
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

  theta_info$model = c(theta_info$model, '')
  theta_info$psd_scale = exp(theta_info$psd_scale_log)

  if (dirichlet) {
    theta_info$var = c(theta_info$var, 'overdisp')
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
  
  config = c(
    config[
    setdiff(names(config), c('verbose', 'dirichlet'))
    ], 
    list(
      verbose = (verbose>1),  
      dirichlet = as.integer(dirichlet))
  )
  
  if(verbose) cat("formatting data..")
  tmb_data = formatHpolData(tmb_data)
  if(verbose) cat("done\n")

 
    configDefaults = list(
      verbose=FALSE, transform_theta=TRUE,
      num_threads = 1
    )

  configDefaults = configDefaults[setdiff(names(configDefaults), names(config))]
  config = c(config, configDefaults)
  controlInner = control_inner

  Sgamma = seq(nrow(tmb_data$XTp)+1, len=nrow(tmb_data$ATp))
  start_beta = rep(1e-2, nrow(tmb_data$XTp))
  start_gamma=    rep(1e-2, nrow(tmb_data$ATp)) 

  start_theta = theta_info$init[!duplicated(theta_info$map)]
  if(config$transform_theta)
    start_theta = log(start_theta)    

  parameters = c(start_beta, start_theta)
  start_parameters = c(start_beta, start_gamma, start_theta)


  if(verbose) cat("getting groups..")
  groups = try(sparsity_grouped(start_parameters, tmb_data, config, verbose))
  if(verbose) cat("done\n")
  if(any(class(groups) == 'try-error')) groups = list()
  config$sparsity = groups$sparsity
  config$groups = groups$groups
  config$group_sparsity = groups$group_sparsity

  config$beta = start_beta
  config$theta = start_theta

  getAdfunHere = function() {getAdfun(start_gamma, data=tmb_data, config=config)}
  env <- new.env(parent = asNamespace("hpolcc"))
  env$start_gamma <- start_gamma
  env$tmb_data    <- tmb_data
  env$config      <- config
  environment(getAdfunHere) <- env
  lockEnvironment(env, bindings = TRUE)
  
  cache = new.env(parent = emptyenv())
  assign("Nfun", 0, cache)
  assign("Ngr", 0, cache)
  assign("gamma_start", start_gamma, cache)


  if(for_dev)
    return(
      list( 
        start_gamma = start_gamma,
        parameters = parameters, 
        parameters_for_sparsity = start_parameters,
        theta_info = theta_info,
        tmb_data = tmb_data,
        terms = terms,
        config = config,
        control = control,
        control_inner = control_inner,
        data = data,
        groups = groups[setdiff(names(groups), "sparsity")],
        adfun = getAdfunHere,
        cache = cache
      )
    )
  






  if(verbose) cat("optimizing")

    mle = trustOptim::trust.optim(
    x = parameters,
    fn = wrappers_outer$fn,
    gr = wrappers_outer$gr,
    method = 'BFGS',
    control = control,
    data=tmb_data, config = config, cache =  cache, controlInner = control_inner
  )

  if(verbose) cat("done")

    mle$extra = try(loglik(
      mle$solution, 
      get("gamma_start", cache), 
      tmb_data, config, control_inner, check=TRUE))

  mle$gamma_hat = mle$extra$solution

  
  return(mle)
}
