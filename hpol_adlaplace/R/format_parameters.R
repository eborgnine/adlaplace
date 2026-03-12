#' parameters_info is a list with elements beta, gamma, theta
#' each needs columns var, name, theta needs log

#' @export
format_parameters = function(x, parameters_info) {

  Ntheta = nrow(parameters_info$theta)
  Nbeta = nrow(parameters_info$beta)
  Ngamma = nrow(parameters_info$gamma)

  betaNames =gsub("_fpoly_", "", parameters_info$beta$name)
  thetaNames = parameters_info$theta$name
  gammaNames = parameters_info$gamma$name

  result = list(
    theta = x[seq(to=length(x), len=Ntheta)],
    beta = x[seq(1, len=Nbeta)]
    )

  names(result$beta) = betaNames
  names(result$theta) = thetaNames

  isLogged = parameters_info$theta$log
  if(any(isLogged)) {
    result$theta[isLogged] = exp(result$theta[isLogged])
  }


  theGamma = data.frame(
    name = gammaNames, 
    value = x[seq(Nbeta+1, len=Ngamma)], 
    term = gsub("(_|[[:digit:]]|global)+$", "", gammaNames),
    ext = gsub("^_", "", regmatches(gammaNames, regexpr("(_|[[:digit:]]|global)+$", gammaNames)))
    )

    gammaCat = strsplit(theGamma$ext, '_')
  gammaCatN <- lapply(gammaCat, function(x) {
    c(NA,x)[seq(to=length(x)+1, by=1, len=2)]
  })
  gammaCat = as.data.frame(do.call(rbind, gammaCatN))
  names(gammaCat) = c('group','index')
  theGamma = cbind(theGamma, gammaCat)


  rownames(theGamma) = theGamma$name
  theGamma$index = as.numeric(theGamma$index)
  theGamma[is.na(theGamma$group), 'group'] = 'global'
  theGamma$group = factor(theGamma$group,
      levels= c('global', sort(setdiff(unique(theGamma$group), 'global'))))

  theGamma = theGamma[order(theGamma$term, theGamma$index, theGamma$group),]

  result$gamma = theGamma[,c('term','group','index','value')]

  result$info = parameters_info
  result

}