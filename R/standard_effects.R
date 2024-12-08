#' Title of the Function
#'
#' @description A brief description of what the function does.
#'
#' @param x A description of the `x` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param y A description of the `y` parameter. Mention its type and purpose (e.g., a numeric vector).
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A description of the return value, including its type (e.g., a numeric vector, a data frame, etc.).
#' @export
#'
#' @details Provide any additional details about the function, such as edge cases, assumptions, or implementation notes.
#'
#' @examples
#' # Basic usage
#' my_function(1:10, 2:11)
#'
bs <- function(x, o = 2, ref_value, range = NULL, nknots = NULL, run_as_is = T, ...) {
  
  stop("bs() is not implemented.")

  l <- c(list(var = deparse(substitute(x)), 
         type = "bs", 
         o = o, 
         ref_value = ref_value,
         nknots = nknots,
         range = range,
         run_as_is = run_as_is,
         prefix = "splines::"),
    list(...))
  
  if(!is.null(range)){
    if(ref_value < range[1] | ref_value > range[2])
      stop("ref_value not in the range of the data for", var, " (bs).")
  }
  
  return(l)
}
