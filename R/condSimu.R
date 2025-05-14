

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
condSimIwp = function(fit, term, newx, Nsim) {
  terms = fit$terms
  termsHere =   grep(term, unlist(lapply(terms, function(xx)
    xx$var)))
  modelHere = unlist(lapply(terms[termsHere], function(xx)
    xx$model))
  groupsHere = lapply(terms[termsHere], function(xx)
    xx$groups)
  groupsHere = groupsHere[unlist(lapply(groupsHere, length))>0]
  knotsHere = lapply(terms[termsHere], function(xx)
    xx$knots)
  knotsHere = knotsHere[unlist(lapply(knotsHere, length))>0][[1]]
  groupVarHere = unique(unlist(lapply(terms[termsHere], function(xx)
    xx$group_var)))
  
    
  if (!(any(modelHere %in% c('hiwp', 'iwp')))) {
    warning("model should be iwp to use condSimIwp")
  }
  
  ref_value = unique(unlist(lapply(terms[termsHere], 
                                   function(xx) xx$ref_value)))
  if (length(ref_value) > 1)
    warning("can't figure out ref_value for ", term, "\n")
  if (missing(newx)) {
    Nvalues = 101
    theRange = range(c(knotsHere, range(fit$obj$env$.data$X[, term]) + ref_value))
    Sx = seq(theRange[1], theRange[2], len = Nvalues)
    if (any(modelHere == 'hiwp')) {
      iwpTerm = fit$terms[[termsHere[modelHere == 'hiwp']]]
      newx = expand.grid(newx = Sx, groupvar = levels(groupsHere[[1]]))
      names(newx) = gsub("^groupvar$", groupVarHere, names(newx))
    } else {
      newx = data.frame(newx = Sx)
    }
  }
  names(newx) = gsub("newx", term, names(newx))
  
  newXA = getNewXA(terms, newx)

  
  if (!is.null(newXA$X)) {
    betaHat = fit$est$param$est[colnames(fit$obj$env$data$X)]
    parametricPart = newXA$X[, grep(term, colnames(newXA$X)), drop = FALSE] %*%
      betaHat[grep(term, names(betaHat))]
  } else {
    parametricPart = 0
  }
  
  # approx gamma | Y ~ N(gamma hat, hessian^(-1))
  
  gammaSim = getGammaSim(fit, term, Nsim)
  
  if (!all(colnames(newXA$A) %in% rownames(gammaSim)))
    warning("some predicted values not in the original data")
  # to do, unconditional simulation
  
  iwpConstr = as.matrix(newXA$A %*% gammaSim[colnames(newXA$A), ])
  
#  toPlot = grep("hrpoly_2_", colnames(newXA$A))
#  matplot(Sx, newXA$A[newx$ecoRegion == '2',toPlot] %*% gammaSim[toPlot, ], type='l')
  
  iwpSim = as.matrix(parametricPart + iwpConstr)
  
  if (length(groupsHere) >= 1) {
    iwpSim2 = cbind(newx, as.data.frame(iwpSim))
    iwpSim = reshape2::acast(
      reshape2::melt(
        iwpSim2,
        id.vars = colnames(newx),
        variable.name = 'sim'
      ),
      as.formula(paste(c(
        colnames(newx), 'sim'
      ), collapse = '~'))
    )
  }
  
  result = list(x = Sx,
                y = iwpSim,
                fixed_mean = drop(parametricPart[newx[[iwpTerm$group_var]] == newx[[iwpTerm$group_var]][1], ]))
  result
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




getGammaSim = function(fit, term, Nsim) {

  if(any(names(fit) == 'obj')) {
    obj = fit$obj
    gammaHat = fit$est$random$est
    precMat = fit$est$random$hessian
  } else {
    obj = fit
    gammaHat = obj$env$last.par.best[grep("^gamma$",
                                          names(obj$env$last.par.best))]
    # obj$env$last.par.best[grep("^theta$",names(obj$env$last.par.best))]
    precMat = obj$env$spHess(par = obj$env$last.par.best,random=TRUE)
  }
  Shere = grep(term, colnames(obj$env$data$A))
  randomMean = gammaHat[Shere]
  
  cholSchur = chol(rje::schur(
    as.matrix(Matrix::forceSymmetric(precMat)),
    Shere,
    Shere,
    setdiff(1:length(randomMean), Shere)
  ))
  cholSchurInv = solve(cholSchur)
  
  gammaSim = randomMean +
    cholSchurInv %*% matrix(rnorm(Nsim * nrow(cholSchur)),
                            ncol = Nsim,
                            dimnames = list(rownames(cholSchur), 1:Nsim))
  gammaSim
}
