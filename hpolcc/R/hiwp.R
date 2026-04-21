#' Hierarchical Integrated Wiener Process Model Term
#'
#' @description Creates and manages hierarchical integrated Wiener process (HIWP) model terms.
#' @name hiwp-class
#' @aliases hiwp
#' @docType class
#' @title Hierarchical Integrated Wiener Process Model Term
#' @importFrom adlaplace iwp
#' @importFrom adlaplace .type_factor_levels
#' @importFrom adlaplace rpoly
#' @importFrom adlaplace fpoly
#'
#' @section Methods:
#' The following methods are available for `hiwp` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for HIWP term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for HIWP term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

setClass("hiwp",
  representation = representation(
    by = "character",
    by_levels = "integer",
    by_labels = "character",
    init = "numeric"
  ),
  contains = "model",
  prototype = prototype(
    by = character(0),
    by_levels = integer(0),
    by_labels = character(0),
    init = numeric(0),
    type = factor("random", levels = adlaplace::.type_factor_levels)
  )
)

#' @export
#' @rdname hiwp-class
#' @param x Variable name
#' @param p Order of the integrated Wiener process (default: 2)
#' @param ref_value Reference value for the basis
#' @param knots Vector of knot locations
#' @param range Range of the data (optional)
#' @param by Grouping variable for hierarchical structure
#' @param init Initial values for theta parameters
#' @param lower Lower bounds for theta parameters
#' @param upper Upper bounds for theta parameters
#' @param parscale Parameter scales for optimization
#' @param boundary_is_random Whether boundary should be treated as random
#' @param include_poly Whether to include polynomial terms
#' @param include_global Whether to include global component
#' @return A list containing hiwp term object and optionally related terms
#' @examples
#' # Example usage:
#' # knots <- seq(0, 1, length.out = 5)
#' # hiwp_term <- hiwp(x = "age", knots = knots, by = "group")
hiwp <- function(
  x, p = 2, ref_value = 0, knots, range = NULL,
  by,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE,
  include_global = TRUE
) {
  init <- rep_len(init, 2 * p + 4)
  lower <- rep_len(lower, 2 * p + 4)
  upper <- rep_len(upper, 2 * p + 4)
  parscale <- rep_len(parscale, 2 * p + 4)

  ref_value <- adlaplace::ref_align(ref_value, knots)

  the_f <- stats::as.formula(paste0("~ 0 + ", x), env=new.env())
  result <- list()
  hiwp_name <- paste(c(x, "hiwp"), collapse = "_")
  result[[hiwp_name]] <- new("hiwp",
    term = x,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    knots = knots,
    by = by,
    by_levels = integer(0),
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
  )

  if (include_global) {
    iwp_name <- paste(c(x, "iwp"), collapse = "_")
    result[[iwp_name]] <- new("iwp",
      term = x,
      formula = the_f,
      p.order = as.integer(p),
      ref_value = ref_value,
      knots = knots,
      init = init[2],
      lower = lower[2],
      upper = upper[2],
      parscale = parscale[2]
    )
  }

  if (include_poly) {
    for (D_poly in seq(1, len = p - 1)) {
      hrpoly_name <- paste(c(x, "hrpoly", D_poly), collapse = "_")
      result[[hrpoly_name]] <- hrpoly(
        x = x, p = D_poly, ref_value = ref_value,
        by = result[[1]]@by,
        init = init[2 + D_poly],
        lower = lower[2 + D_poly],
        upper = upper[2 + D_poly],
        parscale = parscale[2 + D_poly]
      )
    }
  }
  if (include_global & include_poly) {
    if (boundary_is_random) {
      rpoly_name <- paste(c(x, "rpoly"), collapse = "_")

      result[[rpoly_name]] <- adlaplace::rpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    } else {
      fpoly_name <- paste(c(x, "fpoly"), collapse = "_")
      result[[fpoly_name]] <- adlaplace::fpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    }
  }
  result
}

setMethod("design", "hiwp", function(term, data) {

  by_stuff = get_by_levels(term, data)

  term_iwp <- methods::as(term, "iwp")
  A0 <- adlaplace::design(term_iwp, data)
  if(!all(unique(data[[term@by]] %in% term@by_levels))) {
    warning("by levels in data not in the model")
  }
  id_split <- split(1:nrow(data),
    factor(data[[term@by]], levels = term@by_levels),
    drop = FALSE
  )

  A0split <- mapply(
    function(AA, xx) {
      res <- methods::as(AA[xx, , drop = FALSE], "TsparseMatrix")
      cbind(i = xx[1 + res@i], j = res@j + 1, x = res@x)
    },
    xx = id_split, MoreArgs = list(AA = A0), SIMPLIFY = FALSE
  )

  A0combine <- cbind(
    as.data.frame(do.call(rbind, A0split)),
    split = rep(seq(0, len = length(A0split)), unlist(lapply(A0split, nrow)))
  )

  A0combine[, "j2"] <- A0combine[, "j"] + ncol(A0) * (A0combine[, "split"])

  Afinal <- Matrix::sparseMatrix(
    i = A0combine$i, j = A0combine$j2, x = A0combine$x,
    dims = c(nrow(data), ncol(A0) * length(id_split)),
    dimnames = list(
      rownames(data),
      paste0(term@term, "_hiwp_k",
        rep(
          formatC(1:ncol(A0), width = ceiling(log10(ncol(A0))), flag = "0"),
          length(id_split)
        ),
      "_g",
                rep(term@by_labels, each = ncol(A0))
      )
    )
  )

  Afinal
})

setMethod("precision", "hiwp", function(term, data) {
  term = get_by_levels(term, data)

  iwp_precision <- precision(methods::as(term, "iwp"), data)
  result <- Matrix::.bdiag(replicate(length(term@by_levels), iwp_precision))
  dimnames(result) <- list(
  paste0(
    rep(gsub("_iwp_k", "_hiwp_k", colnames(iwp_precision)), length(term@by_levels)),
    "_g",
    rep(term@by_labels, each = ncol(iwp_precision))
  )
)[c(1, 1)]
  result
})

setMethod("theta_info", "hiwp", function(term) {
  # Global and local parameters
  result <- data.frame(
    term = term@term,
    model = "hiwp",
    label = paste(c(term@term, "hiwp"), collapse = "_"),
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )

  result
})

setMethod("beta_info", "hiwp", function(term) {
  # HIWP terms don't have beta parameters
  return(NULL)
})

setMethod("random_info", "hiwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "hiwp",
    label = paste(c(term@term, "hiwp"), collapse = "_"),
    by = term@by_levels,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )

  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$by_labels <- term@by_labels[match(result$by, term@by_levels)]
  result$gamma_label <- paste0(result$label, "_k", bnumPad, "_g", result$by_label)

  result
})

# Register the coercion method properly
methods::setAs(
  "hiwp", "iwp",
  function(from) {
    # Create a new iwp object with the same basic properties
    methods::new("iwp",
      term = from@term,
      formula = from@formula,
      knots = from@knots,
      ref_value = from@ref_value,
      p.order = from@p.order,
      by = character(0), # hiwp has hierarchical structure, iwp doesn't
      init = from@init,
      lower = from@lower,
      upper = from@upper,
      parscale = from@parscale
    )
  }
)

