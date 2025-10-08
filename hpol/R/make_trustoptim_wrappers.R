
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

#' @export
make_trustoptim_wrappers_old <- function(data, 
  config = list(dirichelet=TRUE), 
  obj_fn = objectiveFunctionC,
  obj_fn_noad = objectiveFunctionNoDiff,
  debug=FALSE) {

  # Create environment to store the last evaluated x and result
  cache_env <- new.env()
  cache_env$data = data
  cache_env$config1 = config[setdiff(names(config), 'maxDeriv')]
  if(is.null(config$debugfile)) cache_env$debugfile = 'hpoldebug.rds'

  get_result <- function(x, maxDeriv) {
    result <- obj_fn(x, cache_env$data, 
      c(cache_env$config1, list(maxDeriv = unname(maxDeriv))))
    return(result)
  }
  
  fn_wrapper <- function(x, ...) {
#   get_result(x, maxDeriv=0)$value
    objectiveFunctionNoDiff(x, cache_env$data, cache_env$config1)
  }
  gr_wrapper <- function(x, ...) {
    get_result(x,  maxDeriv=1)$grad
  }
  hs_wrapper <- function(x, ...) {
    result = get_result(x,  maxDeriv=2)
    if(identical(cache_env$config1$debug, TRUE)) { 
      saveRDS(c(result, list(x=x)), file=cache_env$config1$debugfile)  
    }
    result = result$hessian

    if(is.null(cache_env$config1$sparsity)) {
      # save the sparsity pattern
      sparsity = list(
        second=list(
          list(
            i=result@i, j=result@j,
            p=as(result, "CsparseMatrix")@p)
        )
      )
      names(sparsity$second) = c("full", "random")[1+(
        length(x) == nrow(cache_env$data$ATp)
      )]
      cache_env$config1$sparsity = sparsity
    }
  # trustoptim needs general matrix
    as(as(result, 'CsparseMatrix'),'generalMatrix')

  }
  list(fn = fn_wrapper, gr = gr_wrapper, hs = hs_wrapper)
}

