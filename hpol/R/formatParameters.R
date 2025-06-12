
formatParameters = function(x, obj, logscale_theta = FALSE) {
  Ntheta = length(obj$theta_info$name)
  Nbeta = nrow(obj$tmb_data$XTp)
  Ngamma = nrow(obj$tmb_data$ATp)

  result = list(
    theta = x[seq(to=length(x), len=Ntheta)],
    beta = x[seq(1, Nbeta)]
    )
  names(result$beta) = rownames(obj$tmb_data$XTp)
  names(result$theta) = obj$theta_info$names
  if(logscale_theta | any(result$theta < 0)) {
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
    length(x) <- 2  # this will pad with NAs
    x
  })
  gammaCat = as.data.frame(do.call(rbind, gammaCatN))
  names(gammaCat) = c('group','order')
  theGamma = cbind(theGamma, gammaCat)

  rownames(theGamma) = theGamma$name
  theGamma$order = as.numeric(theGamma$order)
  theGamma$group = factor(theGamma$group,
      levels= c('global', sort(as.numeric(setdiff(unique(theGamma$group), 'global')))))

  theGamma = theGamma[order(theGamma$term, theGamma$order, theGamma$group),]

  result$gamma = theGamma[,c('term','group','order','value')]

  result

}