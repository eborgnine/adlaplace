# cache must have gamma_start, Nfun, Nge
#' @export
outer_fn = function(x, data, config, adFunFull, control_inner, cache, parscale = rep(1, length(x))) {
    assign("Nfun", get("Nfun", cache)+1, cache)
    assign("last.par", x, envir=cache)
    xScale = x*parscale
    result=try(loglik(xScale,
      gamma_start = get("gamma_start", envir=cache), 
      data=data, config=config, 
        adFunFull = adFunFull,
      control=control_inner, 
      deriv=0))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA)}
        cat(c(get("Nfun", cache), result$minusLogLik, NA, formatC(x, format = "f", digits = 5), '\n'), file = get('file', cache), append=TRUE)
      }
      if(length(result$solution)) assign("gamma_start", result$solution, envir=cache)
      result$minusLogLik
    }
#' @export
outer_gr = function(x, data, config, adFunFull, control_inner, cache, parscale = rep(1, length(x))) {
    assign("Ngr", get("Ngr", cache)+1, cache)
    xScale = x*parscale
    result= try(loglik(xScale,
        gamma_start = get("gamma_start", envir=cache), 
        adFunFull = adFunFull,
        data=data, config=config, control=control_inner))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA, deriv = list(dL = NA))}
        cat(c(get("Ngr", cache), result$minusLogLik, sqrt(sum(result$deriv$dL^2)), formatC(x, format = "f", digits = 5), '\n'), file = get('file', cache), append=TRUE)
      }
  if(length(result$solution)) assign("gamma_start", result$solution, envir=cache)
  result$deriv$dL * parscale
}
 


