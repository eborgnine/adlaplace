#' @export
condSimIwp = function(fit,
                   term,
                   x,
                   Nsim = 500,
                   Nvalues = 1001) {
  termsHere =   grep(term, unlist(lapply(fit$terms, function(xx)
    xx$var)))
  modelHere = unlist(lapply(fit$terms[termsHere], function(xx)
    xx$model))
  
  if(!(all(modelHere == 'iwp'))) {
    warning("model should be iwp to use condSimIwp")
  }
  
    ref_value = unique(unlist(lapply(fit$terms[termsHere], function(xx)
      xx$ref_value)))
    if (length(ref_value) > 1)
      warning("can't figure out ref_value for ", term, "\n")
    
    theRange = range(fit$obj$env$.data$X[, term]) + ref_value
    Sx = seq(theRange[1], theRange[2], len = Nvalues)
    newDf = data.frame(newx = Sx)
    names(newDf) = term

  newXA = hpoltest:::getNewXA(fit$terms, newDf)
  
  
  betaHat = fit$est$param$est[grep("theta", names(fit$est$param$est), invert =
                                     TRUE)]
  names(betaHat) = colnames(fit$obj$env$.data$X)
  
  if (!is.null(newXA$X)) {
    parametricPart = newXA$X[, grep(term, colnames(newXA$X)), drop = FALSE] %*%
      betaHat[grep(term, names(betaHat))]
  } else {
    parametricPart = 0
  }
  
  # approx gamma | Y ~ N(gamma hat, hessian^(-1))
  randomMean = fit$est$random$est
  precMat = Matrix::forceSymmetric(fit$est$random$hessian)
  
  Shere = grep(term, names(randomMean))
  cholSchur = chol(rje::schur(as.matrix(precMat), Shere, Shere, setdiff(1:length(randomMean), Shere)))
  cholSchurInv = solve(cholSchur)
  
  
  gammaSim = cholSchurInv %*% matrix(rnorm(Nsim * nrow(cholSchur)),
                                     ncol = Nsim,
                                     dimnames = list(rownames(cholSchur), 1:Nsim))
  if (!all(colnames(newXA$A) %in% rownames(gammaSim)))
    warning("some predicted values not in the original data")
  # to do, unconditional simulation
  
  iwpConstr = as.matrix(newXA$A %*% gammaSim[colnames(newXA$A), ])
  
  iwpSim = parametricPart + iwpConstr
  
  result = list(x = Sx, y = iwpSim, mean = parametricPart)
}

condSimHiwp = function(fit,
                      term,
                      Nsim = 500,
                      Nvalues = 1001) {
 
  simIwp = condSimIwp(fit, term, Nsim, Nvalues)
   
}
  
condSimIid = function(fit,
                      term,
                      Nsim = 500) {
  
  termsHere =   grep(term, unlist(lapply(fit$terms, function(xx)
    xx$var)))
  modelHere = unlist(lapply(fit$terms[termsHere], function(xx)
    xx$model))
  
  if(!(all(modelHere == 'iid'))) {
    warning("model should be iid to use condSimIid")
  }
    gammaHere = grep(term, fit$gamma_info$var)
    Sx1 = colnames(fit$obj$env$data$A)[gammaHere]
    Sx = gsub(paste0("(factor[(])?", term, "[)]?"), "", Sx1)
    term2 = gsub(Sx[1], "", Sx1[1])
    newDf = data.frame(newx = factor(Sx))
    names(newDf) = term
}

# function not done yet
}
                      