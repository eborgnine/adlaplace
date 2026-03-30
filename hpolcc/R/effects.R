#' Functions to Specify Terms of the Formula
#'
#' @description This suite of functions implements various models (fixed and random effects), 
#' polynomials, integrated Wiener processes, and hierarchical extensions. The primary 
#' interface for these models is the `f()` function, which allows specification of 
#' model terms via the `model` argument.
#'
#' @param x A numeric vector representing the data to which the model is applied.
#' @param model A character string specifying the model of model term. Options include:
#'   - `"fpoly"`: Fixed polynomial effect terms.
#'   - `"rpoly"`: Random polynomial effect terms.
#'   - `"hrpoly"`: Random hierarchical polynomial effect terms.
#'   - `"iwp"`: Integrated Wiener process terms (random).
#'   - `"hiwp"`: Hierarchical integrated Wiener process terms (random).
#' @param p A numeric value specifying the polynomial degree (default: varies by function).
#' @param ref_value A numeric value specifying the reference point for the basis functions.
#' @param knots A numeric vector specifying the knot locations for spline-based models.
#' @param range A numeric vector of length 2 specifying the range of the data (optional).
#' @param group_var A factor or grouping variable used to define hierarchical structures (optional).
#' @param include_global A logical value indicating whether to include a global component in hierarchical models (default: `TRUE`).
#' @param rpoly_p A numeric value specifying the degree for (random effects) polynomial components to add to the model as a separate term.
#' @param fpoly_p A numeric value specifying the degree for (fixed effects) polynomial components to add to the model as a separate term.
#' @param hrpoly_p A numeric value specifying the degree for hierarchical (random effects) polynomial components to add to the model as a separate term.
#' @param ... Additional arguments passed to other methods or functions.
#'
#' @return A list or matrix, depending on the function, representing the design matrix, precision matrix, or parameter information for the specified model term.
#'
#' @details
#' The `f()` function serves as the primary entry point for specifying model terms. Internally, it dispatches to specific functions based on the `model` argument:
#' - `"fpoly"`: Calls the `fpoly()` function to specify fixed polynomial effect terms.
#' - `"rpoly"`: Calls the `rpoly()` function to specify random polynomial effect terms.
#' - `"hrpoly"`: Calls the `hrpoly()` function to specify hierarchical random polynomial effect terms.
#' - `"iwp"`: Calls the `iwp()` function to construct integrated Wiener process terms.
#' - `"hiwp"`: Calls the `hiwp()` function for hierarchical integrated Wiener process terms.
#'
#' @examples
#' # Example usage with f()
#' term <- f(x = 0:10, model = "iwp", ref_value = 5, knots = seq(0,10,2))
#'
#' @export
f <- function(x, model="iid",  ...) {
  model_fun = get(model)
  x <- deparse(substitute(x))
  model_fun(x, ...)
}

.my_beta_init <- 0
.my_beta_lower <- -Inf
.my_beta_upper <- Inf
.my_beta_parscale <- 1

.my_theta_init <- 0.02
.my_theta_lower <- 1e-9
.my_theta_upper <- Inf
.my_theta_parscale <- 1

ref_align <- function(ref_value, knots) {
  knots[which.min(abs(knots - ref_value))]
}

# Collect terms from formula
collect_terms <- function(formula) {
  term_labels <- attr(stats::terms(formula), "term.labels")

  terms_1 <- lapply(term_labels, function(lab) {
    if (grepl("f\\(", lab)) {
      tryCatch({
        term_obj <- eval(parse(text = lab))
        if (is.null(term_obj)) {
          stop(paste("Failed to parse term:", lab))
        }
        return(term_obj)
      }, error = function(e) {
        warning(paste("Error parsing term", lab, ":", e$message))
        return(list(linear(lab)))
      })
    } else {
      result = list(linear(lab))
      names(result) = result[[1]]@term
      return(result)
    }
  })
  terms_all = do.call(c, terms_1)

  terms_all
}


get_by_levels = function(term, data) {
  if(length(term@by)) {
  if(!length(term@by_levels)) {
    unique_values <- unique(data[[term@by]])
    if(is.numeric(unique_values)) {
      unique_values_string <- 
          formatC(unique_values, width = ceiling(log10(max(unique_values))), flag = "0")
    } else {
      unique_values_string = as.character(unique_values)
    }
    term@by_levels = unique_values
    term@by_labels = unique_values_string
  }
  }
  term
}
