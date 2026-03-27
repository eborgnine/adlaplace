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
    name = paste(c(prefix, x, "linear"), sep = "_"),
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
    group_var = deparse(substitute(group_var)), ,
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
      result[[poly_name]] <- new("rpoly",
        var = x,
        name = poly_name,
        p = p - 1,
        f = the_f,
        ref_value = ref_value
      )
    } else {
      result[[poly_name]] <- new("fpoly",
        var = x,
        name = poly_name,
        p = p - 1,
        f = the_f,
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
    group_var = deparse(substitute(group_var)),
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
      result[[hrpoly_name]] <- new("hrpoly",
        var = x, p = D_poly, ref_value = ref_value,
        name = hrpoly_name,
        f = the_f,
        group_var = deparse(substitute(group_var)),
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

      result[[rpoly_name]] <- new("rpoly",
        var = x,
        name = rpoly_name,
        p = p - 1,
        ref_value = ref_value,
        f = the_f
      )
    } else {
      fpoly_name <- paste(c(prefix, x, "fpoly"), collapse = "_")
      result[[fpoly_name]] <- new("fpoly",
        var = x,
        name = fpoly_name,
        p = p - 1,
        ref_value = ref_value,
        f = the_f
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
  D[,1:ncol(D),drop=FALSE]

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
  the_levels <- unique(data[[term@group_var]])
  id_split <- split(1:nrow(data), 
    factor(data[[term@group_var]], levels = the_levels), 
  drop = FALSE)
  mm <- c(0, sapply(id_split, length)) |> cumsum()

  A0 <- drop(poly(
    data[[term@var]] - term@ref_value, 
    raw = TRUE, degree = term@p
  )[,term@p])
    
  Afinal <- Matrix::Matrix(0, nrow = nrow(data), ncol = length(the_levels)) |> as("TsparseMatrix")

  for (k in seq_along(id_split)) {
    Afinal[id_split[[k]], k] <- A0[id_split[[k]]]
  }
  colnames(Afinal) <- paste(
    term@name,
    the_levels,
    sep = "_"
  )

  Afinal
})

#' @rdname design
#' @export
setMethod("design", "hiwp", function(term, data) {
  term_iwp <- as(term, "iwp")
  A0 <- design(term_iwp, data)
  the_levels <- unique(data[[term@group_var]])
  id_split <- split(1:nrow(data),
    factor(data[[term@group_var]], levels = the_levels),
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
        rep(names(id_split), each = ncol(A0)),
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

#' @rdname effects_and_utilities
#' @export
 iwpPrecision <- function(term) {
   as(compute_weights_precision(knots = term@knots), "TsparseMatrix")
 }

#' @rdname effects_and_utilities
#' @export
 iwpTheta <- function(theta_info, term) {
   result <- data.frame(
     var = term$var, model = term$model,
     global = TRUE, order = NA, init = term$init,
     lower = term$lower, upper=term$upper,
     parscale = ifelse(is.null(term$parscale), 1, term$parscale)
   )
   result$name <- apply(result[, c("var", "model")], 1, paste, collapse = "_")
   return(result)
 }

# 
# 


#' @rdname effects_and_utilities 
hiwpPrecision <- function(term){
  list2env(term, envir = environment())
  Matrix::.bdiag(replicate(include_global+ngroups, iwpPrecision(term)))
}

#' @rdname effects_and_utilities 
hiwpTheta <- function(theta_info, term){

  # Global and local parameters
  result_global = data.frame(
    var=term$var, model=term$model, 
    global=TRUE, 
    order=NA, 
    init= term$init[1], 
    lower = term$lower[1], 
    upper = term$upper[1],
    parscale = ifelse(is.null(term$parscale), 1, term$parscale[1])
  )
  result_global$name = paste(result_global$var, result_global$model, "GLOBAL", sep='_')

  # Local (hrpoly) parameters
  if(term$hrpoly_p > 0) {
    result_local = data.frame(
      var=term$var, model=term$model, 
      global=FALSE, 
      order=seq(1, term$hrpoly_p), 
      init= term$init_hrpoly,
      lower = term$lower_hrpoly, 
      upper = term$upper_hrpoly,
      parscale = ifelse(is.null(term$parscale_hrpoly), rep(1, term$hrpoly_p), term$parscale_hrpoly)
    )
    result_local$name = paste(result_local$var, result_local$model, "LOCAL", result_local$order, sep='_')
    result = rbind(result_global, result_local)
  } else {
    result = result_global
  }

return(result)
}





#' @rdname effects_and_utilities 
fpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, degree = p)
  D <- D[,1:ncol(D),drop=F]
  colnames(D) <- paste0(term$var, c('', seq(from=1, by=1, len=ncol(D)-1)))
  
  D
}




#

#' @rdname effects_and_utilities 
rpolyDesign <- function(term, data){
  list2env(term, envir = environment())
  D <- poly(data[[var]]-ref_value, raw=T, degree = p)
  D[,1:ncol(D),drop=F]
}

#' @rdname effects_and_utilities 
rpolyPrecision <- function(term){
  as(matrix(1), "TsparseMatrix")
}

#' @rdname effects_and_utilities 
rpolyTheta <- function(theta_info, term){

  result = data.frame(
    var=term$var, model=term$model, 
    global=NA, 
    order=seq(1, term$p),
  init = term$init,
  parscale = ifelse(is.null(term$parscale), rep(1, term$p), rep_len(term$parscale, term$p)))

    result$name = apply(result[,c('var','model','order')],1,paste, collapse='_')

return(result)
}



#' @rdname effects_and_utilities 
hrpolyPrecision <- function(term){
  list2env(term, envir = environment())
  pp <- (include_global+ngroups)*p
  Matrix::sparseMatrix(i=1:pp, j=1:pp, rep = "T") |> as("TsparseMatrix")
}

#' @rdname effects_and_utilities 
hrpolyTheta <- function(theta_info, term){

  result = data.frame(
    var=term@var, model=term$model, 
    global=FALSE, 
    order=seq(1, term$p),
  init = term$init,
  lower = term$lower, upper=term$upper,
  parscale = ifelse(is.null(term$parscale), rep(1, term$p), rep_len(term$parscale, term$p)))
  result$name = apply(result[,c('var','model','order')],1,paste, collapse='_')
  return(result)

}






#' @rdname effects_and_utilities 
iidDesign <- function(term, data){
  list2env(term, envir = environment())
  ff <- paste0("~ 0 + factor(", var, ")") |> formula()
  Matrix::sparse.model.matrix(ff, data) |> as("TsparseMatrix")
}

#' @rdname effects_and_utilities 
iidPrecision <- function(term){
  list2env(term, envir = environment())
  Matrix::sparseMatrix(x=1, i=1:n, j=1:n, rep = "T") |> as("TsparseMatrix")
}



iidTheta <- function(theta_info, term){
  result = data.frame(
    var = term@var, model=term$model,global=NA, order=NA, init=term$init,
    parscale = ifelse(is.null(term$parscale), .my_theta_parscale, term$parscale)
  )
  result$name = paste0(result$var, '_', result$model)
  return(result)
}

#' @rdname effects_and_utilities 
fpolyTheta <- function(theta_info, term){
  result = data.frame(
    var = term@var, model=term$model, global=NA, 
    order=seq(1, term$p), init=term$init, lower=term$lower, upper=term$upper,
    parscale = ifelse(is.null(term$parscale), rep(.my_beta_parscale, term$p), rep_len(term$parscale, term$p))
  )
  result$name = apply(result[,c('var','model','order')],1,paste, collapse='_')
  return(result)
}
