#' @export
make_trustoptim_wrappers <- function(data, 
    config = list(dirichelet=TRUE), 
    obj_fn = objectiveFunctionC) {
  # Create environment to store the last evaluated x and result
  cache_env <- new.env()
  cache_env$last_x <- NULL
  cache_env$last_result <- NULL
  cache_env$data = data
  cache_env$config1 = config[setdiff(names(config), 'hessMax')]
  cache_env$config2 = list()
  if(!is.null(config$hessMax)) cache_env$config2$hessMax = config$hessMax
  
  get_result <- function(x) {
    if (!is.null(cache_env$last_x) && all(x == cache_env$last_x)) {
      return(cache_env$last_result)
    } else {
      result <- obj_fn(x, cache_env$data, c(cache_env$config1, cache_env$config2))
      cache_env$last_x <- x
      cache_env$config2$hessMax = ceiling(1.1*length(result$hessian$x))
      result$hessian <- as(
      	result$hessian,
      	'dgCMatrix')
      cache_env$last_result <- result
      return(result)
    }
  }
  
  fn_wrapper <- function(x, ...) get_result(x)$value
  gr_wrapper <- function(x, ...) get_result(x)$grad
  hs_wrapper <- function(x, ...) get_result(x)$hessian
  
  list(fn = fn_wrapper, gr = gr_wrapper, hs = hs_wrapper)
}