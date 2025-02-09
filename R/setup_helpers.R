#' Initialize Theta Parameters for Hierarchical Models
#'
#' @description This function initializes the \eqn{\theta} parameters for hierarchical models by processing terms in a formula object. It supports multiple term types such as fixed polynomials (\code{fpoly}), integrated Wiener processes (\code{iwp}), hierarchical versions of these models (\code{hiwp}, \code{hrpoly}), and independent random effects (\code{iid}). It also handles reference values and degrees for polynomial terms.
#'
#' @param formula A formula object specifying the hierarchical model structure.
#'
#' @return A data frame containing information about the \eqn{\theta} parameters, including:
#' \item{variable}{The variable name associated with the term.}
#' \item{level}{The level of the term (e.g., \code{GLOBAL} or \code{LOCAL} for hierarchical terms).}
#' \item{model}{The type of model for the term (e.g., \code{iwp}, \code{poly}).}
#' \item{degree}{The degree for polynomial terms (if applicable) or p of IWP(p).}
#' \item{ref_value}{The reference value used for the term.}
#' \item{initial_value}{The initial value assigned to the \eqn{\theta} parameter.}
#'
#' @details The function iteratively processes the terms in the formula. For each term, it expands hierarchical and random polynomial terms into their detailed components, assigning appropriate reference values and initial values for \eqn{\theta} parameters. Unsupported term types are skipped.
#'
#' @examples
#' # Example usage:
#' formula <- ~ fpoly(x, degree = 2) + iwp(z) + hrpoly(y, degree = 3)
#' theta_init <- thetaInitInfo(formula)
#' print(theta_init)
#'
#' @export
thetaInitInfo <- function(formula) {
  
  # Check inputs
  if (!is(formula, "formula")) stop("formula must be a formula.")
  terms <- collectTerms(formula)
  
  k <- 1
  while(k <= length(terms)){
    term <- terms[[k]]
    if(term$run_as_is){k <- k+1; next}
    if(term$model %in% "fpoly"){k <- k+1; next}
    terms <- c(terms[1:k], addFPoly(term), addRPoly(term), terms[-(1:k)])
    k <- k+1
  }
  
  lapply(terms, \(term){
    
    list2env(term, envir = environment())
    
    df <- if(term$model == "iwp"){
        data.frame(variable = var, level = NA, model = model, degree = p, ref_value = ref_value)
      }else if(term$model == "hiwp"){
        data.frame(variable = var, 
                   level = if(include_global){c("GLOBAL", "LOCAL")}else{"LOCAL"}, 
                   model = "iwp", 
                   degree = p, 
                   ref_value = ref_value)
      }else if(term$model == "iid"){
        data.frame(variable = var, level = NA, model = model, degree = NA, ref_value = NA)
      }else if(term$model == "rpoly"){
        data.frame(variable = var, level = NA, model = "poly", degree = 1:p, ref_value = ref_value)
      }else if(term$model == "hrpoly"){
        data.frame(variable = var, 
                   level = rep(if(include_global){c("GLOBAL", "LOCAL")}else{"LOCAL"}, each = p),
                   model = "poly", 
                   degree = rep(1:p, 1+include_global), 
                   ref_value = term$ref_value)
      }else{
        return(NULL)
      }
    
    df$initial_value <- .my_theta_init
    return(df)
  }) |> do.call(what = "rbind")
}

