#' IID Random Effects Term
#'
#' @description Creates and manages IID (independent and identically distributed) random effects terms.
#' @name iid-class
#' @aliases iid
#' @exportClass iid
#' @export
setClass("iid",
  slots = list(),
  contains = "model",
  prototype = prototype(
    ref_value = numeric(0),
    p.order = as.integer(0),
    knots = numeric(0),
    type = factor("random", levels = .type_factor_levels),
    ad_fun = "random_diagonal",
    ad_kind = "random"
  )
)

#' @param x Variable name (typically a factor).
#' @param init Initial value for theta parameter.
#' @param lower Lower bound for theta parameter.
#' @param upper Upper bound for theta parameter.
#' @param parscale Parameter scale for optimization.
#' @param log Whether theta is optimized on the log scale.
#' @return An `iid` term object (in a named list).
#' @rdname iid-class
#' @export
iid <- function(x,
                init = .my_theta_init,
                lower = .my_theta_lower,
                upper = .my_theta_upper,
                parscale = .my_theta_parscale,
                log = TRUE) {
  x <- strip_term_name(as.character(x))
  result <- list(methods::new("iid",
    term = x,
    label = paste(x, "iid", sep = "_"),
    formula = stats::as.formula(paste0("~ 0 + ", x), env = new.env()),
    init = init,
    lower = lower,
    upper = upper,
    parscale = parscale,
    log = log
  ))
  names(result) <- result[[1]]@term
  result
}

#' @describeIn iid-class Creates design matrix for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "iid", function(term, data) {
  if (is.numeric(data[[term@term]])) {
    data[[term@term]] <- factor(data[[term@term]])
  }
  result <- methods::as(Matrix::sparse.model.matrix(term@formula, data), "TsparseMatrix")
  colnames(result) <- gsub(paste0("^", term@term), paste0(term@term, "_iid_"), colnames(result))
  result
})

#' @describeIn iid-class Creates precision matrix for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "iid", function(term, data) {
  # Identity matrix for iid terms
  term_here <- data[[term@term]]
  if (is.factor(term_here)) {
    n <- nlevels(term_here)
    labels_here <- levels(term_here)
  } else {
    labels_here <- as.character(unique(term_here))
    n <- length(labels_here)
  }
  result <- Matrix::Diagonal(n, 1)
  result@x <- rep(1, n)
  dimnames(result) <- list(paste0(term@term, "_iid_", labels_here))[c(1, 1)]
  result
})

#' @describeIn iid-class Extracts random effects information for iid term
#' @param term An iid term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "iid", function(term, data) {
  term_values <- data[[term@term]]
  basis_labels <- if (is.factor(term_values)) {
    levels(term_values)
  } else {
    as.character(sort(unique(term_values)))
  }

  result <- data.frame(
    term = term@term,
    model = "iid",
    label = term@label,
    by = NA,
    by_labels = NA,
    basis = basis_labels,
    order = NA
  )
  result$gamma_label <- paste(result$label, result$basis, sep = "_")

  result
})

#' @describeIn iid-class Extracts theta parameter information for iid term
#' @param term An iid term object
#' @export
setMethod("theta_info", "iid", function(term) {
  result <- data.frame(
    term = term@term, model = "iid",
    label = term@label,
    init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type,
    log = term@log
  )
  return(result)
})

