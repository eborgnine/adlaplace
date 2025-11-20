#' @export
formatParameters = function(x, obj, logscale_theta = obj$config$transform_theta) {
  Ntheta = length(obj$theta_info$name)
  Nbeta = nrow(obj$tmb_data$XTp)
  Ngamma = nrow(obj$tmb_data$ATp)

  result = list(
    theta = x[seq(to=length(x), len=Ntheta)],
    beta = x[seq(1, Nbeta)]
    )
  names(result$beta) = rownames(obj$tmb_data$XTp)
  names(result$theta) = obj$theta_info$name
  if(logscale_theta ) {
    result$theta = exp(result$theta)
  }

  gammaNames = rownames(obj$tmb_data$ATp)
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

  result

}