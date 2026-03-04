
  getTermsPred = function(terms) {
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

  predDf = list()
  for(D in names(predSeq)) {
    newdf = data.frame(x=predSeq[[D]])
    colnames(newdf) = D
    if(!all(is.na(Sgroup[[D]]))) {
      newdf[[Sgroup[[D]]]]= rep(NA, nrow(newdf))
    }
    predDf[[D]] = newdf
  }


  list(predSeq = predSeq, predDf = predDf, Sgroup = Sgroup, Sref = Sref)
  }

condSimGamma = function(fit, Nsim) {

  halfH = fit$extra$halfHinv
  # note tcrossprod(halfH) = Hinv
  Ngamma = nrow(halfH)

  gammaHat = fit$opt$solution

  simInd = matrix(rnorm(Nsim * Ngamma), Ngamma, Nsim)

  simGamma1 = as.matrix(halfH %*% simInd)

  simGamma = simGamma1 + matrix(gammaHat, length(gammaHat), ncol(simGamma1))
  rownames(simGamma) = names(gammaHat)

  simGamma
}


#' @export
condSim = function(fit, term, newx, Nsim=500) {
  termsHere =   grep(term, unlist(lapply(fit$terms, function(xx)
    xx$var)))
  modelHere = unlist(lapply(fit$terms[termsHere], function(xx)
    xx$model))
  if(any(modelHere %in% c('hiwp','iwp'))) {
    result = condSimIwp(fit, term, newx, Nsim)
  }
  if(any(modelHere %in% c('iid'))) {
    result = condSimIid(fit, term, Nsim)
  }

  result  
  
}
#' @export
condSimIwp = function(fit, terms, 
  parameters_info, Nsim, newx, newConstr) {
 
  # fit needs full_parameters, inner (all of it), 

  beta = fit$fullParameters[seq(1, len=nrow(parameters_info$beta))]

  simGamma = condSimGamma(fit, Nsim)
  rownames(simGamma) = parameters_info$gamma$name

  termsPred = getTermsPred(terms)

  if(!missing(newx)) {
    termsPred$predDf = newx
  }

  if(missing(newConstr)) {
    newConstr = termsPred$Sref # replace by new constraints
  }

  newXA = mapply(
#    hpolcc:::getNewXA,
    getNewXA,
    df= termsPred$predDf,
    MoreArgs = list(
      terms = terms
    ),
    SIMPLIFY=FALSE
  )

  simF = mapply(
    function(XA, gamma, beta) {
      namesBothBeta = intersect(names(beta), colnames(XA$X))    
      namesBothGamma = intersect(rownames(gamma), colnames(XA$A))    
      fixedPart = 
        XA$X[,namesBothBeta, drop=FALSE] %*% 
        beta[namesBothBeta]
  
      randomPart = 
        XA$A[,namesBothGamma, drop=FALSE] %*% 
        gamma[namesBothGamma,,drop=FALSE]

      randomPart + fixedPart[,rep(1, ncol(randomPart)), drop=FALSE]
    },
    XA = newXA,
    MoreArgs = list(beta = beta, gamma = simGamma)
  )

  list(sim=simF, x = termsPred, gamma=simGamma, XA = newXA)

 }

condSimIid = function(fit, term, Nsim) {
  termsHere =   grep(term, unlist(lapply(fit$terms, function(xx)
    xx$var)))
  modelHere = unlist(lapply(fit$terms[termsHere], function(xx)
    xx$model))
  
  if (!(all(modelHere == 'iid'))) {
    warning("model should be iid to use condSimIid")
  }
  gammaHere = grep(term, fit$gamma_info$var)
  Sx1 = colnames(fit$obj$env$data$A)[gammaHere]
  Sx = gsub(paste0("(factor[(])?", term, "[)]?"), "", Sx1)

  gammaSim = getGammaSim(fit, term, Nsim)
  result = list(x = Sx,
                y = gammaSim,
                fixed_mean = 0)
  result
}



