#' Functions to Specify Terms of the Formula (and Associated Utility Functions)
#' @rdname effects_and_utilities
#' @description This suite of functions implements various models (fixed and random effects), polynomials, integrated Wiener processes, and hierarchical extensions. 
#' The primary interface for these models is the `f()` function, which allows specification of model terms via the `model` argument. 
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
#' Utility functions:
#' - **Design Functions (`*Design`)**: Construct design matrices tailored to specific modeling terms.
#' - **Precision Matrix Functions (`*Precision`)**: Construct precision matrices corresponding to specific modeling terms.
#' - **Theta Functions (`*Theta`)**: Handle the variance parameters associated with random effects.
#'
#' These utility functions are specifically used within the `hnlm` framework and are not exported.
#'
#' @examples
#' These were generated with chatGPT and have not been reviewed. See vignette for a human generated example.
#' 
#' # Example usage with f()
#' term <- f(x = 0:10, model = "iwp", ref_value = 5, knots = seq(0,10,2))
#' term <- f(x = 0:10, model = "hiwp", ref_value = 5, knots = seq(0,10,2), group_var = "city")
#'
#' @export
f <- function(x, 
model = c("iwp", "hiwp", "fpoly", "rpoly", "hrpoly", "iid"), ...
) {
  model <- match.arg(model)
  x <- deparse(substitute(x))

  switch(model,
    iwp = iwp(x, ...),
    hiwp = hiwp(x, ...),
    fpoly = fpoly(x, ...),
    rpoly = rpoly(x, ...),
    hrpoly = hrpoly(x, ...),
    iid = iid(x, ...),
    stop("Unknown model")
  )
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

#' @rdname effects_and_utilities
#' @export
linear <- function(x, prefix = NULL) {
  new("linear",
    var = x,
    name = paste(c(prefix, x, "linear"), collapse = "_"),
    f = formula(paste0("~ 0 + ", x))
  )
}

#' @rdname effects_and_utilities
#' @export
iid <- function(x,
                init = .my_theta_init,
                lower = .my_theta_lower,
                upper = .my_theta_upper,
                parscale = .my_theta_parscale,
                prefix = NULL) {
  new("iid",
    var = x,
    name = paste(c(prefix, x, "iid"), collapse = "_"),
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1],
    f = formula(paste0("~ 0 + ", prefix, x))
  )
}

#' @rdname effects_and_utilities
#' @export
fpoly <- function(x, p = 2, ref_value = 0,
                  init = .my_beta_init,
                  lower = .my_beta_lower,
                  upper = .my_beta_upper,
                  parscale = .my_beta_parscale,
                  prefix = NULL) {
  new("fpoly",
    var = x,
    name = paste(c(prefix, x, "fpoly"), collapse="_"),
    f = formula(paste0("~ 0 + ", prefix, x)),
    p = p,
    ref_value = ref_value,
    init = rep_len(init, p),
    lower = rep_len(lower, p),
    upper = rep_len(upper, p),
    parscale = rep_len(parscale, p)
  )
}

#' @rdname effects_and_utilities
#' @export
rpoly <- function(x, p = 2, ref_value, sd = Inf, prefix = NULL) {
  new("rpoly",
    var = x,
    name = paste(c(prefix, x, "rpoly"), collapse = "_"),
    f = formula(paste0("~ 0 + ", prefix, x)),
    p = p,
    ref_value = ref_value,
    sd = rep_len(sd, p)
  )
}

#' @rdname effects_and_utilities
#' @export
hrpoly <- function(
  x,
  p = 1,
  ref_value,
  group_var,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  prefix = NULL
) {
  new("hrpoly",
    var = x,
    name = paste(c(prefix, x, "hrpoly", p), collapse = "_"),
    f = formula(paste0("~ 0 + ", prefix, x)),
    p = p,
    ref_value = ref_value,
    group_var = group_var,
    groups = integer(0),  # Will be set later when data is available
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
  )
}

#' @rdname effects_and_utilities
#' @export
iwp <- function(
  x, p = 2,
  ref_value, knots, range = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE,
  prefix = NULL
) {
  the_f <- formula(paste0("~ 0 + ", prefix, x))
  ref_value <- ref_align(ref_value, knots)
  result <- list()
  iwp_name <- paste(c(prefix, x, "iwp"), collapse = "_")
  result[[iwp_name]] <- new("iwp",
    var = x,
    name = iwp_name,
    f = the_f,
    p = p,
    ref_value = ref_value,
    knots = knots,
    range = range,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
  )
  if (include_poly) {
    poly_name <- paste(c(prefix, x, "poly"), collapse = "_")
    if (boundary_is_random) {
      result[[poly_name]] <- rpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    } else {
      result[[poly_name]] <- fpoly(
        x = x,
        p = p - 1,
        ref_value = ref_value
      )
    }
  }
  result
}


#' @rdname effects_and_utilities
#' @export
hiwp <- function(
  x, p = 2, ref_value, knots, range = NULL,
  group_var,
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
    var = x,
    name = hiwp_name,
    f = the_f,
    p = p,
    ref_value = ref_value,
    knots = knots,
    range = range,
    group_var = group_var, 
    groups = integer(0),  # Will be set later when data is available
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
  )

  if (include_global) {
    iwp_name <- paste(c(prefix, x, "iwp"), collapse = "_")
    result[[iwp_name]] <- new("iwp",
      name = iwp_name,
      var = x,
      f = the_f,
      p = p,
      ref_value = ref_value,
      knots = knots,
      range = range,
      init = init[2],
      lower = lower[2],
      upper = upper[2],
      parscale = parscale[2]
    )
  }

  if (include_poly) {
    for (D_poly in seq(1, len = p - 1)) {
      hrpoly_name <- paste(c(prefix, x, "hrpoly", D_poly), collapse = "_")
      result[[hrpoly_name]] <- hrpoly(
        x = x, p = D_poly, ref_value = ref_value,
        group_var = result[[1]]@group_var,
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

#' @rdname effects_and_utilities 
#' Generic design function
#'
#' @param term A model term object
#' @param data A data frame containing the data
#' @return A design matrix
#' @rdname design
setGeneric("design", function(term, data) {
  standardGeneric("design")
})

#' @rdname design
#' @export
setMethod("design", "iwp", function(term, data) {
  res <- local_poly(
    knots = term@knots - term@ref_value,
    refined_x = data[[term@var]] - term@ref_value, p = term@p
  )
  res <- methods::as(res, "TsparseMatrix")

  bnumPad <- formatC(1:ncol(res),
    width = max(ceiling(c(1, log10(ncol(res)))), na.rm = TRUE),
    flag = "0"
  )
  dimnames(res) <- list(
    rownames(data),
    paste(term@name, bnumPad, sep = "_")
  )
  res
})
#' @rdname design
#' @export
setMethod("design", "rpoly", function(term, data) {
  if (term@p == 0) {
    return(NULL)
  }

  D <- poly(data[[term@var]]-term@ref_value, raw=TRUE, degree = term@p)
  D <- D[,1:ncol(D),drop=FALSE]
  colnames(D) <- paste(term@name, 1:ncol(D), sep="_")
  D
})
#' @rdname design
#' @export
setMethod("design", "fpoly", function(term, data) {
  if (term@fpoly_p == 0) {
    return(NULL)
  }
  D <- poly(data[[term@var]]-term@ref_value, degree = term@p)
  D <- D[,1:ncol(D),drop=F]
  colnames(D) <- paste0(term@var, c('', seq(from=1, by=1, len=ncol(D)-1)))
  D
})

#' @rdname design
#' @export
setMethod("design", "iid", function(term, data){
  ff <- paste0("~ 0 + factor(", term@var, ")") |> formula()
  Matrix::sparse.model.matrix(ff, data) |> as("TsparseMatrix")
})

#' @rdname design
#' @export
setMethod("design", "linear", function(term, data){
    res <- Matrix::sparse.model.matrix(term@f, data, drop.unused.levels = FALSE)
    if (is.factor(data[[term@var]])) {
      res <- res[, -1, drop = FALSE]
    }
    res
})


#' @rdname design
#' @export
setMethod("design", "hrpoly", function(term, data) {

  id_split <- split(1:nrow(data), 
    factor(data[[term@group_var]], levels = term@groups), 
  drop = FALSE)
  mm <- c(0, sapply(id_split, length)) |> cumsum()

  A0 <- drop(poly(
    data[[term@var]] - term@ref_value, 
    raw = TRUE, degree = term@p
  )[,term@p])
    
  Afinal <- Matrix::Matrix(0, nrow = nrow(data), ncol = length(term@groups)) |> as("TsparseMatrix")

  for (k in seq_along(id_split)) {
    Afinal[id_split[[k]], k] <- A0[id_split[[k]]]
  }
  colnames(Afinal) <- paste(
    term@name,
    term@groups_string,
    sep = "_"
  )

  Afinal
})

#' @rdname design
#' @export
setMethod("design", "hiwp", function(term, data) {
  term_iwp <- as(term, "iwp")
  A0 <- design(term_iwp, data)
  if(!all(unique(data[[term@group_var]] %in% term@groups))) {
    warning("groups in data not in the model")
  }
  id_split <- split(1:nrow(data),
    factor(data[[term@group_var]], levels = term@groups),
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
      paste(term@name,
        rep(term@groups_string, each = ncol(A0)),
        rep(
          formatC(1:ncol(A0), width = ceiling(log10(ncol(A0))), flag = "0"),
          length(id_split)
        ),
        sep = "_"
      )
    )
  )

  Afinal
})

#' Generic precision function
#'
#' @param term A model term object
#' @param data A data frame containing the data
#' @return A precision matrix
#' @rdname precision
setGeneric("precision", function(term, data) {
  standardGeneric("precision")
})

#' @rdname precision
#' @export
setMethod("precision", "iwp", function(term, data) {
  result = as(compute_weights_precision(knots = term@knots), "TsparseMatrix")
  dimnames(result) = list(paste(
    term@name,
    1:ncol(result),
    sep="_"
  ))[c(1,1)]
  result
})

#' @rdname precision
#' @export
setMethod("precision", "hiwp", function(term, data) {
  unique_values <- term@groups
  if (!length(unique_values)) {
    unique_values <- unique(data[[term@var]])
  }
  iwp_precision <- precision(as(term, "iwp"), data)
  result <- Matrix::.bdiag(replicate(length(unique_values), iwp_precision))
dimnames(result) <- list(
  paste(
    rep(colnames(iwp_precision), length(unique_values)),
    rep(unique_values, each = ncol(iwp_precision)),
    sep = "_"
  )
)[c(1, 1)]
  result
})

#' @rdname precision
#' @export
setMethod("precision", "rpoly", function(term, data) {
  if (term@p == 0) {
    return(NULL)
  }
  Matrix::Diagonal(term@p, 1)
})

#' @rdname precision
#' @export
setMethod("precision", "hrpoly", function(term, data) {
  if (term@p == 0) {
    return(NULL)
  }
  unique_values <- term@groups
  if (!length(unique_values)) {
    unique_values <- unique(data[[term@var]])
  }
  Matrix::Diagonal(length(unique_values), 1)
})

#' @rdname precision
#' @export
setMethod("precision", "fpoly", function(term, data) {
  # Fixed effects don't have precision matrices
  NULL
})

#' @rdname precision
#' @export
setMethod("precision", "iid", function(term, data) {
  # Identity matrix for iid terms
  n <- length(unique(data[[term@var]]))
  Matrix::Diagonal(n, 1)
})

#' @rdname precision
#' @export
setMethod("precision", "linear", function(term, data) {
  # Linear terms don't have precision matrices
  NULL
})

#' Generic theta_info function
#'
#' @param theta_info Additional theta information (not currently used but kept for compatibility)
#' @param term A model term object
#' @return A data frame with parameter information for the term
#' @rdname theta_info
setGeneric("theta_info", function(term) {
  standardGeneric("theta_info")
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "iwp", function(term) {
  result <- data.frame(
    var = term@var, model = "iwp", name = term@name,
    global = TRUE, order = NA, init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale
  )
  return(result)
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "hiwp", function(term) {
  # Global and local parameters
  result <- data.frame(
    var = term@var,
    model = "hiwp",
    name = term@name,
    global = TRUE,
    order = NA,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )


  return(result)
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "iid", function(term) {
  result <- data.frame(
    var = term@var, model = "iid", name = term@name,
    global = TRUE, order = NA, init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale
  )
  return(result)
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "fpoly", function(term) {
  # not thetas for fpoly
  NULL
  })

#' @rdname theta_info
#' @export
setMethod("theta_info", "rpoly", function(term) {
  NULL
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "hrpoly", function(term) {
  result <- data.frame(
    var = term@var, model = "hrpoly", name = term@name,
    global = FALSE, order = term@p,
    init = term@init,
    lower = term@lower,
    upper = term@upper,
    parscale = term@parscale
  )
  return(result)
})

#' @rdname theta_info
#' @export
setMethod("theta_info", "linear", function(term) {
  # Linear terms don't have random effects parameters
  return(NULL)
})

#' Generic gamma_info function
#'
#' @param term A model term object
#' @param data A data frame containing the data
#' @return A data frame with gamma information for the term
#' @rdname gamma_info
setGeneric("gamma_info", function(term, data) {
  standardGeneric("gamma_info")
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "iwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)
  order <- NA

  result <- expand.grid(
    var = term@var,
    model = "iwp",
    name = term@name,
    groups = NA,
    groups_string = NA,
    basis = basis,
    order = term@p
  )
  bnumPad <- formatC(basis,
    width = max(ceiling(c(1, log10(max(basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$gamma_name <- paste(result$name, bnumPad, sep = "_")

  result
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "hiwp", function(term, data) {
  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    var = term@var,
    model = "hiwp",
    name = term@name,
    groups = term@groups,
    groups_string = term@groups_string,
    basis = basis,
    order = term@p
  )
  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$gamma_name <- paste(result$name, result$groups_string, bnumPad, sep = "_")

  result
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "fpoly", function(term, data) {
  
# no gammas for fpoly
NULL
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "hrpoly", function(term, data) {
  basis <- NA

  result <- expand.grid(
    var = term@var,
    model = "hrpoly",
    name = term@name,
    groups = term@groups,
    groups_string = term@groups_string,
    basis = basis,
    order = term@p
  )
  result$gamma_name <- paste(result$name, result$groups_string, sep = "_")
  result
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "rpoly", function(term, data) {
  order <- seq_len(term@p)
  basis <- NA

  result <- expand.grid(
    var = term@var,
    model = "rpoly",
    name = term@name,
    groups = NA,
    groups_string = NA,
    basis = basis,
    order = order
  )
  result$gamma_name <- paste(result$name, result$order, sep = "_")

  result
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "iid", function(term, data) {
    
  result <- expand.grid(
    var = term@var,
    model = "iid",
    name = term@name,
    groups = NA,
    groups_string = NA,
    basis = sort(unique(data[[term@var]])),
    order = NA
  )
  result$gamma_name = paste(result$name, result$basis, sep="_")

  result
})

#' @rdname gamma_info
#' @export
setMethod("gamma_info", "linear", function(term, data) {
  
  NULL})
