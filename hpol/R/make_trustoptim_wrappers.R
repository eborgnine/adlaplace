# cache must have gamma_start, Nfun, Nge
#' @export
wrappers_outer = list( 
  fn = function(x, data, config, controlInner, cache) {
    assign("Nfun", get("Nfun", cache)+1, cache)
    assign("last.par", x, envir=cache)
    result=try(loglik(x,
      gamma_start = get("gamma_start", envir=cache), 
      data=data, config=config, control=controlInner, 
      deriv=0))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA)}
        cat(c(0, get("Nfun", cache), result$minusLogLik, NA, x, '\n'), file = get('file', cache), append=TRUE)
      }
      assign("gamma_start", result$solution, envir=cache)
      result$minusLogLik
    },
  gr = function(x, data, config, controlInner, cache) {
    assign("Ngr", get("Ngr", cache)+1, cache)
    result= try(loglik(x,
        gamma_start = get("gamma_start", envir=cache), 
        data=data, config=config, control=controlInner))
      if("file" %in% ls(cache)) {
        if('try-error' %in% class(result) ) {result = list(minusLogLik = NA, deriv = list(dL = NA))}
        cat(c(1, get("Ngr", cache), result$minusLogLik, sqrt(sum(result$deriv$dL^2)), x, '\n'), file = get('file', cache), append=TRUE)
      }
  assign("gamma_start", result$solution, envir=cache)
  result$deriv$dL
}
 
 )

hessianGamma = function(x, data, config) {
  res = objectiveFunctionHessian(x, data, config)

  hesDfPairs = data.frame(
    j = config$sparsity$random$jLong,
    i = config$sparsity$random$i,
    x = res$off_diag)
  hesDfDiag = data.frame(
    j = c(config$sparsity$random$SjNotPairOrDiag, config$sparsity$random$onesJ),
    x = c(res$colSum[1+2*config$sparsity$random$SjNotPairOrDiag], res$diag)
  )
  hesDfDiag$i = hesDfDiag$j

  hesDfAll1 = rbind(hesDfPairs, hesDfDiag[,names(hesDfPairs)])

  Ngamma = length(x)
  hessianR1 = Matrix::sparseMatrix(
    i = pmax(hesDfAll1$i, hesDfAll1$j),
    j = pmin(hesDfAll1$i, hesDfAll1$j),
    x = hesDfAll1$x,
    index1=FALSE, symmetric=TRUE,
    dims = c(Ngamma, Ngamma)
  )

# elements with one necessary off-diagonal, the diagonal j was computed with forward 2
  hesDfOnes = data.frame(
    j = config$sparsity$random$onesJ,
    i = config$sparsity$random$onesI,   
    x= res$colSum[1+2*config$sparsity$random$onesJ]- apply(hessianR1[, 1+config$sparsity$random$onesJ], 2, sum)
  )

  hesDfAll2 = rbind(hesDfAll1, hesDfOnes[,names(hesDfAll1)])
  hessianR2 = Matrix::sparseMatrix(
    i = pmax(hesDfAll2$i, hesDfAll2$j),
    j = pmin(hesDfAll2$i, hesDfAll2$j),
    x = hesDfAll2$x,
    index1=FALSE, symmetric=TRUE
  )

# needed one more off-diagonal but diagonal wasn't computed.  
  hesDfInotJ = data.frame(
    i = config$sparsity$random$SiNotJ, j=config$sparsity$random$SiNotJ,
    x = res$colSum[1+2*config$sparsity$random$SiNotJ] - apply(hessianR2[, 1+config$sparsity$random$SiNotJ], 2, sum)
  )

  hesDfAll = rbind(hesDfAll2, hesDfInotJ[,names(hesDfAll2)])
  hessianR = Matrix::sparseMatrix(
    i = pmax(hesDfAll$i, hesDfAll$j), 
    j = pmin(hesDfAll$i, hesDfAll$j),
    x = hesDfAll$x,
    index1=FALSE, symmetric=TRUE
  )
  as(as(hessianR, 'CsparseMatrix'),'generalMatrix')
}

#' @export
wrappers_gamma = list( 
  fn = objectiveFunctionNoDiff,
  gr = objectiveFunctionGrad,
  hs = hessianGamma
)
