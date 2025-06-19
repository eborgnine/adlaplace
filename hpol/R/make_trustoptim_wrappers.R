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
    if(debug) saveRDS(c(result, list(x=x)), file=cache_env$debugfile)  

    if(is.null(cache_env$config1$sparsity)) {
      # save the sparsity pattern
      hess2 = as(as(result, 'TsparseMatrix'), 'generalMatrix')
      sparsity = data.frame(
        i=result@i, j=result@j
      )
      reordered = order( sparsity$j,  sparsity$i)
      config1$sparsity = sparsity[reordered,]
    }

    as(as(result, 'CsparseMatrix'),'generalMatrix')

  }
  list(fn = fn_wrapper, gr = gr_wrapper, hs = hs_wrapper)
}