# cache must have gamma_start, Nfun, Nge
#' @export
wrappers_outer = list( 
  fn = function(x, data, config, adFunFull, control_inner, cache) {
    assign("Nfun", get("Nfun", cache)+1, cache)
    assign("last.par", x, envir=cache)
    result=try(loglik(x,
      gamma_start = get("gamma_start", envir=cache), 
      data=data, config=config, 
        adFunFull = adFunFull,
      control=control_inner, 
      deriv=0))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA)}
        cat(c(0, get("Nfun", cache), result$minusLogLik, NA, x, '\n'), file = get('file', cache), append=TRUE)
      }
      assign("gamma_start", result$solution, envir=cache)
      result$minusLogLik
    },
  gr = function(x, data, config, adFunFull, control_inner, cache) {
    assign("Ngr", get("Ngr", cache)+1, cache)
    result= try(loglik(x,
        gamma_start = get("gamma_start", envir=cache), 
        adFunFull = adFunFull,
        data=data, config=config, control=control_inner))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA, deriv = list(dL = NA))}
        cat(c(1, get("Ngr", cache), result$minusLogLik, sqrt(sum(result$deriv$dL^2)), x, '\n'), file = get('file', cache), append=TRUE)
      }
  assign("gamma_start", result$solution, envir=cache)
  result$deriv$dL
}
 
 )


