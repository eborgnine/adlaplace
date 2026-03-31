#' Hierarchical Integrated Wiener Process Model Term
#'
#' @description Creates a hierarchical integrated Wiener process (HIWP) model term.
#'
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
#' @param prefix Optional prefix for term names
#'
#' @return A list containing hiwp term object and optionally related terms

# HIWP class definition
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
           type = factor("random", levels = .type_factor_levels)
         )
)


hiwp <- function(
  x, p = 2, ref_value = 0, knots, range = NULL,
  by,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE,
  include_global = TRUE,
  prefix = NULL
) {
  init <- rep_len(init, 2 * p + 4)
  lower <- rep_len(lower, 2 * p + 4)
  upper <- rep_len(upper, 2 * p + 4)
  parscale <- rep_len(parscale, 2 * p + 4)

  ref_value <- ref_align(ref_value, knots)

  the_f <- formula(paste0("~ 0 + ", prefix, x))
  result <- list()
  hiwp_name <- paste(c(prefix, x, "hiwp"), collapse = "_")
  result[[hiwp_name]] <- new("hiwp",
    term = x,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    knots = knots,
    by = by, 
    by_levels = integer(0),  # Will be set later when data is available
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
    # type is already set in prototype
  )

  if (include_global) {
    iwp_name <- paste(c(prefix, x, "iwp"), collapse = "_")
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
      # type is already set in iwp prototype
    )
  }

  if (include_poly) {
    for (D_poly in seq(1, len = p - 1)) {
      hrpoly_name <- paste(c(prefix, x, "hrpoly", D_poly), collapse = "_")
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
      rpoly_name <- paste(c(prefix, x, "rpoly"), collapse = "_")

      result[[rpoly_name]] <- rpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    } else {
      fpoly_name <- paste(c(prefix, x, "fpoly"), collapse = "_")
      result[[fpoly_name]] <- fpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    }
  }
  result
}

# Design matrix for hiwp terms
setMethod("design", "hiwp", function(term, data) {

  by_stuff = get_by_levels(term, data)

  term_iwp <- as(term, "iwp")
  A0 <- design(term_iwp, data)
  if(!all(unique(data[[term@by]] %in% term@by_levels))) {
    warning("by levels in data not in the model")
  }
  id_split <- split(1:nrow(data),
    factor(data[[term@by]], levels = term@by_levels),
    drop = FALSE
  )

  A0split <- mapply(
    function(AA, xx) {
      res <- as(AA[xx, , drop = FALSE], "TsparseMatrix")
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

# Precision matrix for hiwp terms
setMethod("precision", "hiwp", function(term, data) {
  term = get_by_levels(term, data)

  iwp_precision <- precision(as(term, "iwp"), data)
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

# Theta info for hiwp terms
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


  return(result)
})

# Beta info for hiwp terms
setMethod("beta_info", "hiwp", function(term) {
  # HIWP terms don't have beta parameters
  return(NULL)
})

# Gamma info for hiwp terms
setMethod("random_info", "hiwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "hiwp",
    label = paste(c(term@term,"hiwp"), collapse = "_"),
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


