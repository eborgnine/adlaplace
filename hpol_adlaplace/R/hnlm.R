#' Fit (some) hierarchical non-linear models (that we like)
#'
#' @description This function fits a hierarchical model using the specified formula and data, incorporating case-crossover designs, fixed and random effects, and precision matrices for random effects. It allows for flexible inclusion of stratification variables, time variables, and complex random effect structures.
#'
#' @param formula A formula object specifying the model to be fitted.
#' @param data A data frame containing the variables specified in the formula and any additional variables required for the model.
#' @param cc_design An object specifying the case-crossover design, including stratification and time variables. Defaults to the output of `ccDesign()`.
#' @param dirichlet if TRUE, fit dirichlet-multinomial model, with time-within-strata gamma random effect, otherwise multinomial.
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
    num_threads = 1, dirichlet=TRUE,
    dirichlet_init = 0.1, Ngroups = 1e4
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
  strat_time_vars <- unique(c(cc_design$strat_vars, cc_design$time_var))

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
  terms = lapply(terms1, hpolcc:::getExtra, data=data, cc_matrix=cc_matrix)
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


  fPolyRandom = unlist(lapply(terms[isFpoly], '[[', "boundary_is_random") )
  isBoundary = isFpoly[fPolyRandom]
  XfPolyRandom = XfPoly[fPolyRandom]
  XfPolyFixed = XfPoly[!fPolyRandom]

  if(any(fPolyRandom) & is.null(config$prec_boundary)){
    config$prec_boundary = 0
  }

    # random boundary X's go in A
  QfPoly = lapply(XfPolyRandom, function(xx, prec) Matrix::Diagonal(ncol(xx), x=prec), prec= config$prec_boundary)
  isBoundary = isFpoly
  Xlist = XasIs


  Alist = c(XfPolyRandom, Arandom)
  Qall = c(QfPoly, Qs)
  Q =  Matrix::bdiag(Qall[!sapply(Qall, is.null)])


  Xlist = c(XasIs, XfPolyFixed)


# A, Gamma
  if (length(Alist)) {
    A = do.call(cbind, Alist) |> as("TsparseMatrix")
  } else {
    A = matrix(nrow = 0, ncol = 0) |> as("TsparseMatrix")
  }

  gamma_setup <- lapply(terms[c(isBoundary, isHrpoly, isRandom)], hpolcc:::getGammaSetup)
  gamma_info = do.call(rbind, gamma_setup)
  gamma_info$global = as.logical(pmin(gamma_info$group == 'GLOBAL', TRUE, na.rm=TRUE))

  missingGammaInfo = setdiff(colnames(A), gamma_info$name)
  missingGammaA = setdiff(gamma_info$name, colnames(A))
  if(length(missingGammaInfo)) {
    warning(" missing gamma info ", paste(missingGammaInfo, collapse=','))
  }
  if(length(missingGammaA)) {
    warning(" missing gamma A ", paste(missingGammaA, collapse=','))
  }

  gamma_info = gamma_info[order(match(gamma_info$name, colnames(A))), ]
  gamma_info$gamma_id = seq(0L, len=nrow(gamma_info))
  if(any(is.na(gamma_info$var))) {
    warning("some columns of design matrix not found in gamma")
  }

  if (length(Xlist)) {
    X = do.call(cbind, Xlist)
    beta_info = data.frame(
      var = rep(names(Xlist), unlist(lapply(Xlist, ncol))),
      name = colnames(X)
    )
  } else {
    X = matrix(nrow = nrow(data), ncol = 0)
    beta_info = data.frame()
  }


# theta setup
  theta_setup = lapply(terms[c(isHrpoly, isRandom)], hpolcc:::getThetaSetup, theta_info = list())

  theta_setup = c(theta_setup, 
    list(data.frame(var='overdisp', 
      model='overdisp', global=TRUE, order=NA, 
      init=config$dirichlet_init, name='overdisp')))

  theta_info = do.call(rbind, theta_setup)
  if(config$transform_theta) {
    theta_info$init = pmax(-15, log(theta_info$init))
    theta_info$log = TRUE
  } else {
    theta_info$log = FALSE
  }
  theta_info$theta_id = seq(0L, len=nrow(theta_info))

  gamma_theta = merge(gamma_info, theta_info, 
    by = c('var','model','order','global'), all.x=TRUE, all.y=TRUE,
    suffixes = c("_gamma","_theta"))

  if(sum(is.na(gamma_theta$theta_id))>10){
    warning("problem matching gamma and theta")
  }


  anyNA = is.na(gamma_theta$theta_id) | is.na(gamma_theta$gamma_id)
#  gamma_theta[anyNA,]

  gamma_theta_both = gamma_theta[!anyNA,]
  #map matrix column theta, row gamma 
  gamma_theta_map = Matrix::sparseMatrix(
    i = gamma_theta_both$gamma_id,
    j = gamma_theta_both$theta_id,
    x = rep(1L, nrow(gamma_theta_both)),
    index1 = FALSE,
    dims = c(nrow(gamma_info), nrow(theta_info))
  )

  tmb_data <- list(
    X = X,
    A = A,
    y = data[[all.vars(formula)[1]]],
    Q =  Q,
    map = gamma_theta_map,
    elgm_matrix = cc_matrix
  )
  tmb_data = hpolcc:::formatHpolData(tmb_data)
  gamma_info$matchA = match(gamma_info$name, rownames(tmb_data$ATp))

  verboseOrig = config$verbose
  config$verbose = config$verbose > 1

  Sgamma = seq(nrow(tmb_data$XTp)+1, len=nrow(tmb_data$ATp))

  if(!length(config$beta)) {
    config$beta = rep(0, nrow(tmb_data$XTp))
  } else {
    config$beta = rep_len(config$beta, nrow(tmb_data$XTp))
  }
  if(!length(config$theta)) {
    config$theta = theta_info$init
  } else {
    config$theta = rep_len(config$theta, length(theta_info$init))
  }
  if(!length(config$gamma)) {
    config$gamma = rep(0, nrow(tmb_data$ATp))
  } else {
    config$gamma = rep_len(config$gamma, nrow(tmb_data$ATp))
  }

  parameters = c(config$beta, config$theta)
  full_parameters = c(config$beta, config$gamma, config$theta)

  parameters_info = list(
    beta = beta_info,
    gamma= gamma_info,
    theta = theta_info
  )

  if(verboseOrig) {
    cat("getting groups...")
  }

  forGroups = tmb_data$ATp %*% tmb_data$elgm_matrix  
  forGroups@x = rep(1, length(forGroups@x))
  config$groups = adlaplace::adFun_groups(
    ATp = forGroups, Ngroups=config$Ngroups)
  if(verboseOrig) {
    cat("done.")
  }

  config$package = "hpolcc"

  cache = new.env()
  assign('gamma', config$gamma, cache)
  # some checks
  if(!all(parameters_info$gamma$name == rownames(tmb_data$ATp))) {
    warning("names of gamma don't match up")
    setdiff(parameters_info$gamma$name, rownames(tmb_data$ATp))
  }

  if(verboseOrig) {
    cat(for_dev)
  }
  if(for_dev) {
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
        cache = cache,
        parameters_info = parameters_info
      )
    )
  }

  if(verboseOrig) {
    cat("getting AD fun...")
  }

  adFun = adlaplace::getAdFun(tmb_data, config, package = "hpolcc")

cache <- new.env(parent = emptyenv())
cache$gamma <- config$gamma

x0 <- c(config$beta, config$theta)

#adlaplace::outer_fn(x=x0, cache=cache, config=forres$config, adFun = adFun)
#adlaplace::outer_gr(x=x0, cache=cache, config=forres$config, adFun = adFun)


  if(verboseOrig) {
    cat("optimizing")
  }
mle = trustOptim::trust.optim(
  x = x0,
  fn = adlaplace::outer_fn,
  gr = adlaplace::outer_gr,
  method = "SR1",
  config = config,
  adFun = adFun,
  cache = cache,
  control = control,
  control_inner = control_inner)


  config$gamma = get("gamma", cache)

  result = list(opt = mle, 
    objects = list(
      tmb_data=tmb_data, config=config, formula=formula, terms = terms,
      parameters_info = parameters_info, gamma_info = gamma_info, 
      control_inner = control$inner)
  )
  if(verboseOrig) {
    cat("done")
  }


if(FALSE) {
   adFun = adlaplace::getAdFun(
    result$objects$tmb_data, 
    result$objects$config, 
    package = "hpolcc")
} 

  result$extra = try(adlaplace::logLikLaplace(
    result$opt$solution, 
    gamma=result$objects$config$gamma, 
    data=result$objects$tmb_data, 
    config=result$objects$config, 
    control = control_inner, 
    adFun = adFun,
    deriv=1))


  result$parameters = try(formatParameters(
    x=result$extra$fullParameters, 
    result$objects$parameters_info))

  if(FALSE) {
    result$hessian_parameters = try(
      numDeriv::jacobian(
        adlaplace::outer_gr,
        x= mle$solution,
        package='hpolcc',
        data = tmb_data, config=config, control_inner=control_inner, adFun=adFun, cache=cache
      )
    )
  }


  result$sample = try(condSimIwp(
    fit=result$extra, 
    terms = result$objects$terms, 
    parameters_info = result$objects$parameters_info,
    Nsim = c(result$objects$config$Nsim, 500)[1]
  ))

  return(result)
}
