#' @export
make_trustoptim_wrappers <- function(data, 
    config = list(dirichelet=TRUE), 
    obj_fn = objectiveFunctionC,
    debug=FALSE) {
  # Create environment to store the last evaluated x and result
  cache_env <- new.env()
  cache_env$data = data
  cache_env$config1 = config[setdiff(names(config), 'maxDeriv')]
  if(is.null(config$debugfile)) cache_env$debugfile = 'hpoldebug.rds'

  get_result <- function(x, maxDeriv) {
      result <- obj_fn(x, cache_env$data, 
        c(cache_env$config1, list(maxDeriv = unname(maxDeriv))))

      if(debug) saveRDS(c(result, list(x=x)), file=cache_env$debugfile)  

      return(result)
  }
  
  fn_wrapper <- function(x, ...) {
    get_result(x, maxDeriv=0)$value
  }
  gr_wrapper <- function(x, ...) {
    get_result(x,  maxDeriv=1)$grad
  }
  hs_wrapper <- function(x, ...) {
    result = get_result(x,  maxDeriv=2)$hessian
 
    if(is.null(cache_env$hessian)) {
      hessT= as(as(result, 'TsparseMatrix'),'generalMatrix')
      # fancy stuff to force matrix to be integers
      cache_env$hessian = cbind(
        i = as.integer(hessT@i), 
        j = as.integer(hessT@j), 
        idx = as.integer(hessT@i + hessT@Dim[1] * hessT@j))

    } 
    as(as(Matrix::sparseMatrix(
        i = cache_env$hessian[,'i'],
        j = cache_env$hessian[,'j'],
        x = result[cache_env$hessian[,'idx']+1],
        index1 = FALSE), 
      'CsparseMatrix'), 'generalMatrix')
  }
  list(fn = fn_wrapper, gr = gr_wrapper, hs = hs_wrapper)
}