#' @export
make_trustoptim_wrappers <- function(data, 
    config = list(dirichelet=TRUE), 
    obj_fn = objectiveFunctionC,
    debug=FALSE) {
  # Create environment to store the last evaluated x and result
  cache_env <- new.env()
  cache_env$last_x <- NULL
  cache_env$last_result <- NULL
  cache_env$data = data
  cache_env$niter = c(f = 0, g=0, h=0)
  cache_env$config1 = config[setdiff(names(config), 'hessMax')]
  cache_env$config2 = list()
  if(!is.null(config$hessMax)) cache_env$config2$hessMax = config$hessMax
  if(is.null(config$debugfile)) cache_env$debugfile = 'hpoldebug.rds'

  get_result <- function(x) {
    if (!is.null(cache_env$last_x) && all(x == cache_env$last_x)) {
      return(cache_env$last_result)
    } else {
      result <- obj_fn(x, cache_env$data, 
        c(cache_env$config1, cache_env$config2))
      cache_env$last_x <- x
      result$hessian <-as(
        as(result$hessian, 'CsparseMatrix'), 
        'generalMatrix')
      cache_env$last_result <- result
      if(debug) saveRDS(c(result, list(x=x)), file=cache_env$debugfile)  
      return(result)
    }
  }
  
  fn_wrapper <- function(x, ...) {
    cache_env$niter['f'] = cache_env$niter['f'] + 1
    get_result(x)$value
  }
  gr_wrapper <- function(x, ...) {
    cache_env$niter['g'] = cache_env$niter['g'] + 1
    get_result(x)$grad
  }
  hs_wrapper <- function(x, ...) {
    cache_env$niter['h'] = cache_env$niter['h'] + 1
    get_result(x)$hessian
  }
  list(fn = fn_wrapper, gr = gr_wrapper, hs = hs_wrapper)
}