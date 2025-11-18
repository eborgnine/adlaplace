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
hnlm <- function(
 formula,
 data,
 cc_design = ccDesign(),
 weight_var,
 tmb_parameters = NULL,
 optim_parameters = list(eval.max = 2000, iter.max = 2000),
 optimizer = c('nlminb', 'optim'),
 config=list(dirichlet = TRUE, boundary_is_random=FALSE),
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
  

  configDefaults = list(
    verbose=FALSE, transform_theta=TRUE,
    num_threads = 1, dirichlet=TRUE
  )

  configDefaults = configDefaults[setdiff(names(configDefaults), names(config))]
  config = c(config, configDefaults)

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
    order(sapply(data[, strat_time_vars, with = FALSE],
     function(xx) length(unique(xx))), decreasing=FALSE)
    ]


  data.table::setorderv(data, strat_time_vars)
  
  # setup the data for case-crossover
  if (config$verbose) {
    cat("setting strata\n")
  }
  
  cc_matrix <- hpolcc:::setStrata(
    cc_design = cc_design, 
    data = data, 
    outcome = all.vars(formula)[1])

  if (config$verbose) {
    cat("collecting terms\n")
  }
  # setup of the design matrices and other parameters
  # terms carries all the information throughout
  terms1 <- hpolcc:::collectTerms(formula)
  terms = lapply(terms1, hpolcc:::getExtra, data=data, cc_matrix=ccMatrix)
  for(D in 1:length(terms)) {
    terms[[D]]$id = D
  }

  isFpoly = which(unlist(lapply(terms, "[[", "model")) == "fpoly")
  isHrpoly = which(unlist(lapply(terms, "[[", "model")) == "hrpoly")
  isAsis = which(unlist(lapply(terms, "[[", "run_as_is")))
  isRandom = setdiff(1:length(terms), c(isFpoly, isHrpoly, isAsis))

  XasIs = lapply(terms[isAsis], function(xx, data) {
    res = Matrix::sparse.model.matrix(xx$f, data, drop.unused.levels=TRUE)
    if(is.factor(data[[xx$var]])) {
      res = res[,-1, drop=FALSE]
    }
    res
  }, data=data)
  names(XasIs) = unlist(lapply(terms[isAsis], '[[', 'var'))

  XfPoly = lapply(terms[isFpoly], function(xx, data) {
    Xsub = as(
      poly(
        data[[xx$var]] - xx$ref_value,
        degree = xx$p,
        raw = TRUE,
        simple = TRUE
      ),
      "TsparseMatrix"
    )
    colnames(Xsub) <- paste0(xx$var, "_fpoly_",seq(
      from = 1,
      len = ncol(Xsub)
    ))
    Xsub
  },
  data=data
)
  names(XfPoly) = unlist(lapply(terms[isFpoly], '[[', 'var'))


  Arandom <- parallel::mclapply(
    terms[c(isHrpoly, isRandom)],
    hpolcc:::getDesign,
    data=data,  mc.cores= config$num_threads)

  Qs = lapply(terms[c(isHrpoly, isRandom)], hpolcc:::getPrecision)
  names(Qs) = names(Arandom) = paste(
    unlist(lapply(terms[isRandom], '[[', 'var')),
    unlist(lapply(terms[isRandom], '[[', 'model')), sep='_')


  if(identical(config$boundary_is_random, TRUE) ) {
    if(is.null(config$prec_boundary)){
      config$prec_boundary = 0
    }
    QfPoly = lapply(XfPoly, function(xx, prec) Matrix::Diagonal(ncol(xx), x=prec), prec= config$prec_boundary)
    isBoundary = isFpoly
    Xlist = XasIs
    Alist = c(XfPoly, Arandom)
    Qs = c(QfPoly, Qs)

  } else {
    isBoundary = c()
    Xlist = c(XasIs, XfPoly)
    Alist = Arandom
    Qall = Qs
  }

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

  Q =  Matrix::bdiag(Qs[!sapply(Qs, is.null)])

  theta_setup = lapply(terms[c(isHrpoly, isRandom)], hpolcc:::getThetaSetup, theta_info = list())

  if (config$dirichlet) {
    if(is.null(config$dirichlet_init)) {
      config$dirichlet_init = 0.1
    }
    theta_setup = c(theta_setup, 
      list(data.frame(var='overdisp', 
        model='overdisp', global=NA, order=NA, 
        init=config$dirichlet_init, name='overdisp')))
  } 
  theta_info = do.call(rbind, theta_setup)

  gamma_setup <- lapply(terms[c(isBoundary, isHrpoly, isRandom)], hpolcc:::getGammaSetup)
  if(!all(unlist(lapply(gamma_setup, nrow)) == unlist(lapply(Alist, ncol)))) {
    warning("gamma and A dont match")
  }
  gamma_info = do.call(rbind, gamma_setup)

  if(length(setdiff(colnames(A), gamma_info$name)) | length(setdiff(gamma_info$name, colnames(A)))) {
    warning("some columns of design matrix not found in gamma")
  }

  gamma_info$global = gamma_info$group == 'GLOBAL'
  gamma_info[gamma_info$model != 'hiwp','global'] = NA

  gamma_theta = merge(gamma_info, theta_info, 
    by = c('var','model','global','order'), all.x=TRUE, all.y=TRUE,
    suffixes = c("_gamma","_theta"))

  gamma_theta$matchA = match(gamma_theta$name_gamma, colnames(A))
  gamma_theta = gamma_theta[order(gamma_theta$matchA),]
  gamma_theta$matchTheta = match(gamma_theta$name_theta, theta_info$name)
  gamma_theta[gamma_theta$var == 'overdisp','matchTheta'] = NA


  NAtheta = which(is.na(gamma_theta$matchTheta))
  if(!all(NAtheta == 
    c(seq(1,len=length(NAtheta) - config$dirichlet), nrow(gamma_theta)[config$dirichlet])
    )) {
    warning("fpoly thetas not at top")
  }

  if(config$transform_theta) {
    theta_info$init = pmax(-15, log(theta_info$init))
  }

  tmb_data <- list(
    X = X,
    A = A,
    y = data[[all.vars(formula)[1]]],
    Q =  Q,
    map = na.omit(gamma_theta$matchTheta)-1,
    cc_matrix = cc_matrix
  )
  tmb_data = formatHpolData(tmb_data)

verboseOrig = config$verbose
  config$verbose = config$verbose > 1

  Sgamma = seq(nrow(tmb_data$XTp)+1, len=nrow(tmb_data$ATp))

  start_beta = rep(0, nrow(tmb_data$XTp))
  start_gamma=    rep(0, nrow(tmb_data$ATp)) 
  start_theta = theta_info$init
  config$beta = start_beta
  config$theta = start_theta

  parameters = c(start_beta, start_theta)
  full_parameters = c(start_beta, start_gamma, start_theta)


  if(verboseOrig) {
    cat("getting groups..")
  }
  groups = sparsity_grouped(x=full_parameters, data=tmb_data, config, verbose=verboseOrig)
  if(verbose) {
    cat("done\n")
  }
  if(any(class(groups) == 'try-error')) {
    groups = list()
  }

  config$sparsity = groups$sparsity
  config$groups = groups$groups
  config$group_sparsity = groups$group_sparsity


  cache = new.env(parent = emptyenv())
  assign("Nfun", 0, cache)
  assign("Ngr", 0, cache)
  assign("gamma_start", start_gamma, cache)


  if(for_dev)
    return(
      list( 
        gamma_start = start_gamma,
        parameters = parameters, 
        parameters_and_gamma = start_parameters,
        theta_info = theta_info,
        tmb_data = tmb_data,
        formula = formula,
        terms = terms, 
        config = config,
        control = control,
        control_inner = control_inner,
        groups = groups[setdiff(names(groups), "sparsity")],
        cache = cache
      )
    )


  if(verboseOrig) cat("optimizing")

    mle = trustOptim::trust.optim(
      x = parameters,
      fn = wrappers_outer$fn,
      gr = wrappers_outer$gr,
      method = 'BFGS',
      control = control,
      data=tmb_data, config = config, cache =  cache, controlInner = control_inner
    )

  if(verboseOrig) cat("done")

    mle$extra = try(loglik(
      mle$solution, 
      get("gamma_start", cache), 
      tmb_data, config, control_inner, check=TRUE))

  mle$gamma_hat = mle$extra$solution
  mle$formula = formula

  return(mle)
}
