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
 optim_parameters = list(eval.max = 2000, iter.max = 2000),
 config=list(dirichlet = TRUE, boundary_is_random=FALSE, transform_theta=TRUE),
 control=list(maxit=2000, start.trust.radius = 0.1, report.level=4, report.freq=1, report.header.freq=10, report.precision=5),
 control_inner = list(report.level=0),
 for_dev = FALSE,
 ...) {

  # Check inputs
  if (!is(formula, "formula"))
    stop("formula must be a formula.")
  #  if (!is.data.frame(data)) stop("data must be a data.frame.")
  

  configDefaults = list(
    verbose=FALSE, transform_theta=TRUE,
    num_threads = 1, dirichlet=TRUE
  )

  configDefaults = configDefaults[setdiff(names(configDefaults), names(config))]
  config = c(config, configDefaults)

  # Order the rows of data appropriately.
  if (is.character(cc_design)) {
    cc_design = hpolcc:::ccDesign(strat_vars = cc_design)
  }
  if (is.null(cc_design$strat_vars) &
    is.null(cc_design$time_var)) {
    stop("Provide a valid stratification (or time) variable.")
  }

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
    res = Matrix::sparse.model.matrix(xx$f, data, drop.unused.levels=FALSE)
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
    # boundary X's go in A
    if(is.null(config$prec_boundary)){
      config$prec_boundary = 0
    }
    QfPoly = lapply(XfPoly, function(xx, prec) Matrix::Diagonal(ncol(xx), x=prec), prec= config$prec_boundary)
    isBoundary = isFpoly
    Xlist = XasIs
    XfPoly = do.call(cbind, XfPoly)
    if(!is.null(XfPoly)) {
      colnames(XfPoly) = gsub("_fpoly_", "_fpoly_GLOBAL_", colnames(XfPoly))
      Alist = c(list(XfPoly), Arandom)
    } else {
      Alist = Arandom      
    }
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

#  if (config$dirichlet) {
    if(is.null(config$dirichlet_init)) {
      config$dirichlet_init = 0.1
    }
    theta_setup = c(theta_setup, 
      list(data.frame(var='overdisp', 
        model='overdisp', global=NA, order=NA, 
        init=config$dirichlet_init, name='overdisp')))

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
  tmb_data = hpolcc:::formatHpolData(tmb_data)
  gamma_info$matchA = match(gamma_info$name, rownames(tmb_data$ATp))


  verboseOrig = config$verbose
  config$verbose = config$verbose > 1

  Sgamma = seq(nrow(tmb_data$XTp)+1, len=nrow(tmb_data$ATp))

  start_beta = rep(0, nrow(tmb_data$XTp))
  config$start_gamma=    rep(0, nrow(tmb_data$ATp)) 
  start_theta = theta_info$init
  config$beta = start_beta
  config$theta = start_theta

  parameters = c(start_beta, start_theta)
  full_parameters = c(start_beta, config$start_gamma, start_theta)

  if(verboseOrig) {
    cat("getting groups...")
  }

  sparsity_raw = hpolcc::sparsity(tmb_data, config)
  # last two are Q and extra
  sparsity_raw = sparsity_raw[seq(1, len=length(sparsity_raw)-2)]

  firstDerivList = lapply(sparsity_raw, '[[', 'grad')
  firstDeriv = Matrix::sparseMatrix(
    i = unlist(firstDerivList),
    j = rep(seq(0, len=length(firstDerivList)), unlist(lapply(firstDerivList, length))),
    x = rep(1, length(unlist(firstDerivList))),
    dims = c(length(full_parameters), length(firstDerivList)),
    index1=FALSE
  )

 if(verboseOrig) {
    cat(",")
  }
  config$groups = adlaplace::adFun_groups(ncol(firstDeriv), firstDeriv)
  sparsity_raw = hpolcc::sparsity(tmb_data, config)
  sparsity_list = adlaplace::group_sparsity(
    data=tmb_data,  
    config=config, 
    sparsity_list = sparsity_raw)

  config = modifyList(config, sparsity_list)

  cache = new.env()
  assign('start_gamma', config$start_gamma, cache)



if(FALSE) {
  theInner = hpolcc::inner_opt(
    config$start_gamma,
    tmb_data,
    control = control_inner, 
    config=config 
  )

  onel = adlaplace::logLik(
    parameters, tmb_data, config,
    config$start_gamma,
    control = control_inner,
    deriv=TRUE, 
    package='hpolcc'
  )
}

  if(for_dev)
    return(
      list( 
        tmb_data = tmb_data,
        config = config,
        formula = formula,
        terms = terms, 
        theta_info = theta_info,
        gamma_info = gamma_info,
        control = control,
        control_inner = control_inner,
        cache = cache
      )
    )


  if(verboseOrig) {
    cat("optimizing")
  }

adFunFull = hpolcc::getAdFun(tmb_data, config,  inner=FALSE)

mle = trustOptim::trust.optim(
  c(config$beta, config$theta),
  adlaplace::outer_fn, adlaplace::outer_gr,
  method='SR1',
  data=tmb_data, config=config, cache=cache, 
  adPack = adFunFull, package = 'hpolcc',
  control_inner = control_inner,
  control = list(report.level=4, report.freq=1)
  )


  result = list(opt = mle, 
    objects = list(
      tmb_data=tmb_data, config=config, formula=formula, terms = terms,
      theta_info = theta_info, gamma_info = gamma_info))


  if(verboseOrig) {
    cat("done")
  }
  result$extra = try(adlaplace::logLik(
    mle$solution, 
    start_gamma=get("start_gamma", cache), 
    data=tmb_data, config=config, control = control_inner, 
    package = 'hpolcc',
    adPack = adFunFull,
    deriv=0))

  if(FALSE) {
  result$hessian_parameters = try(
    numDeriv::jacobian(
      adlaplace::outer_gr,
      x= mle$solution,
      package='hpolcc',
      data = tmb_data, config=config, control_inner=control_inner, adPack=adFunFull, cache=cache
    )
  )
  }

  result$extra$parameters = hpolcc::formatParameters(result$extra$fullParameters, result$objects)

#  mle$parameters = hpolcc::formatParameters(mle$extra$fullParameters, listres, TRUE)

  Ngamma = nrow(result$extra$parameters$gamma)


  Nsim = c(config$Nsim, 500)[1]
  simInd = matrix(rnorm(Nsim * Nparams), Nparams, Nsim)

  simGamma = as.matrix(Matrix::crossprod(result$extra$inner$cholHessian$halfH, simInd))
  rownames(simGamma) = rownames(tmb_data$ATp)

  Sref = unlist(lapply(terms, '[[', "ref_value"))
  Svar = unlist(lapply(terms, '[[', "var"))
  Smodel = unlist(lapply(terms, '[[', "model"))

  isHiwp = which(Smodel %in% c('iwp', 'hiwp'))
  Sref = Sref[isHiwp]
  Svar = Svar[isHiwp]
  Srange = lapply(terms[isHiwp], '[[', 'range')
  predSeq = lapply(Srange, function(xx) seq(min(xx), max(xx), len=100))
  Sgroup = lapply(terms[isHiwp], '[[', 'group_var')

  names(predSeq) =names(Sgroup) = names(Sref) = Svar

  newConstr = Sref # replace by new constraints


  newXA = fixedPart = newConstrIndex= simGlobal = simF= list()

  for(D in names(predSeq)) {
    newDf = data.frame(x = predSeq[[D]], group = NA)
    names(newDf)[2] = Sgroup[D]
    colnames(newDf)[1] = D
    newXA[[D]] = hpolcc:::getNewXA(
      terms = terms,
      df= newDf,
      boundary_is_random= result$objects$config$boundary_is_random
    )

    newColNamesBeta = gsub("_fpoly_", "", names(result$extra$parameters$beta))
    namesBoth = intersect(newColNamesBeta, colnames(newXA[[D]]$X))    

    fixedPart[[D]] = as.vector(newXA[[D]]$X[,namesBoth, drop=FALSE] %*% 
      result$extra$parameters$beta[match(namesBoth, newColNamesBeta)])

    gamma_info[gamma_info$model == 'iwp','global'] = TRUE
    colsA = gamma_info[gamma_info$var == D & gamma_info$global,'name']

    testA = setdiff(colsA, colnames(newXA[[D]]$A))
    if(length(testA)) {
      warning("missing A ", paste(testA, collapse=','))
    }
    testA = setdiff(colsA, rownames(simGamma))
    if(length(testA)) {
      warning("missing A ", paste(testA, collapse=','))
    }

    simF[[D]] = as.matrix(newXA[[D]]$A[,colsA,drop=FALSE] %*% simGamma[colsA,,drop=FALSE])

    simGlobalHere= simF[[D]] + fixedPart[[D]]
    newConstrIndex[[D]] = which.min(abs(predSeq[[D]] - newConstr[D]))
    toSubtract = matrix(simGlobalHere[newConstrIndex[[D]],],
      nrow(simGlobalHere), ncol(simGlobalHere), byrow=TRUE)
    simGlobal[[D]] = simGlobalHere - toSubtract
  }

  Sregions1 = unique(gsub("_[[:digit:]]+$", "", rownames(simGamma)))
  Sregions = setdiff(unique(gsub(".*_", "", Sregions1)), "global")

  result$sample = list(gamma=simGamma, global=simGlobal, groups = Sregions, newXA = newXA, x = predSeq)

# matplot(result$sample$x[[1]], exp(result$sample$global[[1]])-1, type='l', lty=1, col='#00000020', ylim = c(-0.25, 1)) 

  return(result)
}
