#' @include 000.R
#' @include rpoly.R
#' Hierarchical Random Polynomial Model Term
#'
#' @description Creates and manages hierarchical random polynomial model terms.
#' @name hrpoly-class
#' @docType class
#' @exportClass hrpoly
#'
#' @section Methods:
#' The following methods are available for `hrpoly` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for hrpoly term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for hrpoly term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
#' @examples
#' # Create an hrpoly term for variable 'age' with grouping by 'site'
#' term <- hrpoly(x = "age", p = 2, ref_value = 0, by = "site")
#' str(term)
#' # Create sample data
#' set.seed(42)
#' dat <- data.frame(
#'   age = rep(seq(0, 80, by = 10), each = 3),
#'   site = rep(c("A", "B", "C"), times = 9)
#' )
#' by_group(term@by@term, dat)
#' # Compute the design matrix
#' design_mat <- design(term, dat)
#' print(design_mat)
NULL

setClass("hrpoly",
         representation = representation(
           by = "by_group"
         ),
         contains = "model",
         prototype = list(
           by = methods::new("by_group"),
           knots = numeric(0),
           type = factor("random", levels = .type_factor_levels)
         )
)

#' Hierarchical Random Polynomial Term Constructor
#'
#' @description Creates a hierarchical random polynomial model term.
#' @rdname hrpoly-class
#' @param x Variable name.
#' @param p Polynomial degree (default: 1).
#' @param ref_value Reference value for the polynomial.
#' @param by Grouping variable for hierarchical structure.
#' @param init Initial value for theta parameter.
#' @param lower Lower bound for theta parameter.
#' @param upper Upper bound for theta parameter.
#' @param parscale Parameter scale for optimization.
#' @return An `hrpoly` term object.
#' @export
hrpoly <- function(
  x,
  p = 1,
  ref_value = 0,
  by,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale
) {
  # Check all arguments are of length 1
  if (length(x) != 1) stop("x must be a single variable name")
  if (length(p) != 1) stop("p must be a single value")
  if (length(ref_value) != 1) stop("ref_value must be a single value")
  if (length(init) != 1) stop("init must be a single value")
  if (length(lower) != 1) stop("lower must be a single value")
  if (length(upper) != 1) stop("upper must be a single value")
  if (length(parscale) != 1) stop("parscale must be a single value")
  
  # Check that numeric arguments are valid
  if (p <= 0) stop("p must be positive")
  if (any(lower >= upper)) stop("lower bounds must be less than upper bounds")
  if (any(parscale <= 0)) stop("parscale must be positive")

  # Handle by: can be character or by_group
  if (is.character(by)) {
    by <- by_group(term = by)
  }
  # If by is already a by_group, use as is

  methods::new("hrpoly",
    term = x,
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    p.order = as.integer(p),
    ref_value = ref_value,
    by = by,
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale
    # type is already set in prototype
  )
}


#' @rdname hrpoly-class
#' @param term An `hrpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A design matrix for the hrpoly term.
#' @export
setMethod("design", "hrpoly", function(term, data) {

  # Create an rpoly version of the term (without by slot)
  rpoly_term <- methods::new("rpoly",
    term = term@term,
    formula = term@formula,
    p.order = term@p.order,
    ref_value = term@ref_value,
    sd = rep_len(Inf, term@p.order)
  )
  a_base <- design(rpoly_term, data)[,term@p.order]

  a_split = mapply(
    function(x, a_base, id) {
      the_i = which(x==id)
      data.frame(
        i = the_i,
        j_orig = id,
        x = a_base[the_i]
      )
    },
    id = term@by@levels,
    MoreArgs = list(x=data[[term@by@term]], a_base = a_base),
    SIMPLIFY=FALSE
  )
  a_df = do.call(rbind, a_split)
  a_df$j = match(a_df$j_orig, term@by@levels)

  result = Matrix::sparseMatrix(i=a_df$i, j=a_df$j, 
    x = a_df$x, dims = c(length(a_base), length(term@by@levels)),
    dimnames = list(NULL, paste0(
    term@term,
    "_hrpoly_",
    term@p.order,
    "_g",
    term@by@labels
  )))

  result
})

#' @rdname hrpoly-class
#' @param term An `hrpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A precision matrix for the hrpoly term.
#' @export
setMethod("precision", "hrpoly", function(term, data) {
  if (term@p.order == 0) {
    return(NULL)
  }

  result = Matrix::Diagonal(length(term@by@levels), 1)
  dimnames(result) = list(
    paste0(term@term, "_hrpoly_", term@p.order, "_g", term@by@labels)
  )[c(1,1)]
  result

})

#' @rdname hrpoly-class
#' @param term An `hrpoly` term object.
#' @return A data frame containing theta parameter information for the hrpoly term.
#' @export
setMethod("theta_info", "hrpoly", function(term) {
  result <- data.frame(
    term = term@term, model = "hrpoly", 
    label = paste(c(term@term, "hrpoly", term@p.order), collapse = "_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

#' @rdname hrpoly-class
#' @param term An `hrpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return NULL (hrpoly terms don't have beta parameters).
#' @export
setMethod("beta_info", "hrpoly", function(term, data) {
  # HRpoly terms don't have beta parameters (random effects only)
  return(NULL)
})

#' @rdname hrpoly-class
#' @param term An `hrpoly` term object.
#' @param data A data frame containing the variables used in the term.
#' @return A data frame containing random effects information for the hrpoly term.
#' @export
setMethod("random_info", "hrpoly", function(term, data) {
  basis <- NA

  result <- expand.grid(
    term = term@term,
    model = "hrpoly",
    label = paste(c(term@term, "hrpoly", term@p.order), collapse = "_"),
    by = term@by@levels,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )
  result$by_labels <- term@by@labels[match(result$by, term@by@levels)]
  result$gamma_label <- paste0(result$label,  "_g", result$by_labels)
  result
})


