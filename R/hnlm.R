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
#' @import TMB
#'
#' @examples
#' # See vignette for basic usage
#'
#' @useDynLib hpoltest
#' @export
hnlm <- function(formula,
                 data,
                 cc_design = ccDesign(),
                 weight_var,
                 dirichlet = FALSE,
                 tmb_parameters = NULL,
                 optim_parameters = list(eval.max = 2000, iter.max = 2000),
                 optimizer = c('nlminb', 'optim'),
                 for_dev = FALSE,
                 verbose = FALSE,
                 ...) {
  setDT(data)
  
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
  setorderv(data, strat_time_vars)
  
  #  strat_time_vars <- c(cc_design$strat_vars, cc_design$time_var)
  #  new_order <- eval(str2lang(paste0("order(", paste0("data[['", strat_time_vars, "']]", collapse = ", "), ")")))
  #  data <- data[new_order,]
  
  
  # setup the data for case-crossover
  if (verbose)
    cat("setting strata")
  cc_matrix <- setStrata(cc_design = cc_design, data = data)
  if (verbose)
    cat(".\n")
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
  if (verbose)
    cat('collecting terms ')
  while (k <= length(terms)) {
    if (verbose)
      cat(k, ' ')
    
    term <- hpoltest:::getExtra(terms[[k]], data = data, cc_matrix = cc_matrix)
    term$id <- k
    if (!is.factor(data[[term$var]][1]) &&
        !is.character(data[[term$var]][1]) &&
        is.null(term$range))
      range <- term$range <- range(data[[term$var]])
    
    if (term$run_as_is) {
      Xsub <- sparse.model.matrix(term$f, data)
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
    theta_info$psd_scale_log <- c(theta_info$psd_scale_log,
                                  c(theta_setup$psd_scale_log, 0)[1])
    
    # update term with new elements
    terms[[k]] <- term
    k <- k + 1
  } # done k loop
  # final element of theta is the dirichlet SD
  theta_info$var = c(theta_info$var, 'overdisp')
  theta_info$psd_scale = exp(theta_info$psd_scale_log)
  if (dirichlet) {
    theta_info$map = c(theta_info$map, max(theta_info$map) + 1)
    dirichletStart = 0.01
  } else {
    dirichletStart = 0
    theta_info$map = c(theta_info$map, NA)
  }
  theta_info$init <- c(theta_info$init, dirichletStart)
  
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
  
  y <- data[[all.vars(formula)[1]]]
  tmb_data <- list(
    X = X,
    A = A,
    y = y,
    # gamma_nreplicate = gamma_info$nreplicate, # **** when hiwp, reuse the Q matrix for all (split gamma in nreplicate equal parts). gamma_nreplicate=nlevel+1
    # Q = do.call(bdiag, Qs[!unlist(lapply(Qs, is.null))]), #Qs |> .bdiag(),
    Q =  bdiag(Qs[!sapply(Qs, is.null)]),
    gamma_split = gamma_info$split,
    psd_scale = theta_info$psd_scale,
    # theta_id = theta_info$id
    cc_matrix = cc_matrix,
    dirichlet = as.integer(dirichlet)
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
    optim_inline_parameters$lower =  c(-5 / theMax, rep(1e-5, Nthetas))
  if (!'parscale' %in% names(optim_parameters)) {
    optim_parameters$parscale = c(theMax, rep(1/100, Nthetas))
  }
  
  r <- if (length(tmb_parameters$gamma) > 0) {
    "gamma"
  } else{
    NULL
  }
  
  if (for_dev)
    return(
      list(
        X = X,
        A = A,
        gamma_split = gamma_info$split,
        Qs = Qs,
        theta_info = theta_info,
        #new_order = new_order,
        tmb_parameters = tmb_parameters,
        Alist = Alist,
        tmb_data = tmb_data,
        map = map,
        optim_parameters = optim_parameters,
        optim_inline_parameters = optim_inline_parameters,
        r = r,
        dirichlet = dirichlet,
        verbose = verbose
      )
    )
  
  
  # OPTIMIZATION ----
  # # preliminary run fixing the random effects for iwp, hiwp and od
  # # (but not the corresponding random slopes)
  # to_rm_ids <- which(theta_info$model %in% c("od", "iwp", "hiwp"))
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
  
  if (verbose)
    message("making AD function")
  obj <- MakeADFun(
    data = tmb_data,
    parameters = tmb_parameters[c('beta', 'gamma', 'theta')],
    random = r,
    map = map,
    intern = FALSE,
    type = 'ADFun',
    DLL = "hpoltest",
    ...
  )
  if (verbose)
    message("first evaluation")
  firstEval = try(obj$fn(obj$par))
  if (any(class(firstEval) == 'try-error')) {
    warning('errrors in first evaluation')
    return(obj)
  }
  
  
  if (verbose)
    message("beginning optimization")
  if (optimizer[1] == 'nlminb') {
    if (verbose)
      message("nlminb")
    opt <- stats::nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      upper = optim_inline_parameters$upper,
      lower = optim_inline_parameters$lower,
      control = optim_parameters[setdiff(names(optim_parameters),
                                         c('parscale', 'lower', 'upper', 'scale'))]
    )
  } else if (optimizer[1] == 'optim') {
    if (verbose)
      message("optim")
    optim_parameters = optim_parameters[setdiff(names(optim_parameters), c('lower', 'upper', 'method'))]
    if (!length(optim_inline_parameters$method))
      optim_inline_parameters$method = 'L-BFGS-B'
    names(optim_parameters) = gsub("^iter.max$", "maxit", names(optim_parameters))
    opt <- try(stats::optim(
      par = pmin(pmax(obj$par,optim_inline_parameters$lower), optim_inline_parameters$upper),
      fn = obj$fn,
      gr = obj$gr,
      method = optim_inline_parameters$method,
      upper = optim_inline_parameters$upper,
      lower = optim_inline_parameters$lower,
      control = optim_parameters[setdiff(names(optim_parameters), 'eval.max')],
      hessian = TRUE
    ))
    if (any(class(opt) == 'try-error')) {
      warning('errrors in first evaluation')
      return(obj)
    }
    if (FALSE) {
      obj$fn(obj$env$last.par[grep("gamma", names(obj$env$last.par), invert = TRUE)])
      obj$env$last.par[1:20]
      obj$env$last.par[grep("gamma", names(obj$env$last.par), invert = TRUE)]
      obj$env$f(obj$env$last.par)
      thebad = which.min(obj$report()$eta)
      plot(obj$report()$eta)
      obj$fn(c(0.29968561, -0.01796229,  1.28560509,  0.16854090,  0.00001000 ))
      obj$gr(obj$env$last.par[grep("gamma", names(obj$env$last.par), invert=TRUE)])
    }
    opt$control = c(optim_parameters, optim_inline_parameters)
  } else {
    warning("optimizer must be nlminb or optim")
    return(list(
      obj = obj,
      formula = formula,
      terms = terms,
      cc_design = cc_design,
    ))
  }
  if (verbose)
    message("done optimization")
  
  
  fitList = try(formatResult(obj))
  
  
  # Return the result
  return(
    list(
      obj = obj,
      formula = formula,
      terms = terms,
      cc_design = cc_design,
      beta_info = beta_info,
      gamma_info = gamma_info,
      theta_info = theta_info,
      est = fitList,
      #order = new_order,
      opt = opt
    )
  )#, funNoRandom = funNoRandom))
}
