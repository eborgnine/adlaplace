# IWP class definition
setClass("iwp",
  representation = representation(
    # Inherits all slots from model class
  ),
  contains = "model",
  prototype = prototype(
    # Default values for iwp-specific behavior
    by = character(0),
    type = factor("random", levels = .type_factor_levels)
  )
)

#' @title Create an Integrated Wiener Process Term
#'
#' @description
#' Creates an integrated Wiener process (IWP) model term.
#'
#' @param x Variable name
#' @param p Order of the integrated Wiener process (default: 2)
#' @param ref_value Reference value for the basis
#' @param knots Vector of knot locations
#' @param range Range of the data (optional)
#' @param init Initial values for theta parameters
#' @param lower Lower bounds for theta parameters
#' @param upper Upper bounds for theta parameters
#' @param parscale Parameter scales for optimization
#' @param boundary_is_random Whether boundary should be treated as random
#' @param include_poly Whether to include polynomial terms
#'
#' @return A list containing iwp term object and optionally polynomial terms

iwp <- function(
  x, p = 2,
  ref_value = 0, 
  knots, 
  range = NULL,
  init = .my_theta_init,
  lower = .my_theta_lower,
  upper = .my_theta_upper,
  parscale = .my_theta_parscale,
  boundary_is_random = TRUE,
  include_poly = TRUE
) {
  # Check all arguments other than knots are length 1. knots must be length > 1
  if (length(x) != 1) stop("x must be a single variable name")
  if (length(p) != 1) stop("p must be a single value")
  if (length(ref_value) != 1) stop("ref_value must be a single value")
  if (!is.null(range) && length(range) != 2) stop("range must be a vector of length 2")
  if (length(init) != 1) stop("init must be a single value")
  if (length(lower) != 1) stop("lower must be a single value")
  if (length(upper) != 1) stop("upper must be a single value")
  if (length(parscale) != 1) stop("parscale must be a single value")
  if (length(knots) < 2) stop("knots must have length >= 2")
  if (length(boundary_is_random) != 1) stop("boundary_is_random must be a single value")
  if (length(include_poly) != 1) stop("include_poly must be a single value")

  the_f <- formula(paste0("~ 0 + ", x))
  result <- list()
  iwp_name <- paste("model", x, sep = "_")

  result[[iwp_name]] <- new("iwp",
    term = x,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    knots = knots,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
    # type is already set in prototype, no need to repeat
  )
  if (include_poly) {
    poly_name <- paste(c(x, "poly"), collapse = "_")
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

# Design matrix for iwp terms
setMethod("design", "iwp", function(term, data) {
  refined_x <- data[[term@term]] - term@ref_value
  basis <- local_poly(term@knots, refined_x, term@p.order)
  result <- basis[,1:ncol(basis),drop=FALSE]

  knots_string = formatC(seq.int(ncol(result)), 
    width = ceiling(log10(ncol(result))), flag = "0")

  colnames(result) <- paste0(term@term, "_iwp_k", knots_string)
  result
})

# Precision matrix for iwp terms
setMethod("precision", "iwp", function(term, data) {
  result = Matrix::Matrix(compute_weights_precision(term@knots))
  
  knots_string = formatC(seq.int(nrow(result)), 
    width = ceiling(log10(nrow(result))), flag = "0")

  dimnames(result) = list(
    paste0(term@term, "_iwp_k", knots_string)
  )[c(1,1)]
  result
})

# Theta info for iwp terms
setMethod("theta_info", "iwp", function(term) {
  result <- data.frame(
    term = term@term, model = "iwp", 
    label = paste(c(term@term, "iwp"), collapse = "_"),
    init = term@init,
    lower = term@lower, upper = term@upper,
    parscale = term@parscale,
    type = term@type
  )
  return(result)
})

# Beta info for iwp terms
setMethod("beta_info", "iwp", function(term) {
  # IWP terms don't have beta parameters
  return(NULL)
})

# Gamma info for iwp terms
setMethod("random_info", "iwp", function(term, data) {

  basis <- seq(1, len = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "iwp",
    label = paste(c(term@term, "iwp"), collapse = "_"),
    by = NA,  # iwp doesn't have hierarchical structure
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )

  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$by_labels <- NA  # iwp doesn't have by_labels
  result$gamma_label <- paste0(result$label, "_k", bnumPad)

  result
})

# Adaptation of the local_poly function from the OSplines packages --------
local_poly <- function(knots, refined_x, p) {
  if (min(knots) >= 0) {
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >=
          knots[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x[j] -
            knots[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] -
            knots[i + 1])^(p - k)) / (factorial(k) * factorial(p -
            k)))
        }
      }
    }
  } else if (max(knots) <= 0) {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] &
          refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] -
            knots_neg[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] -
            knots_neg[i + 1])^(p - k)) / (factorial(k) *
            factorial(p - k)))
        }
      }
    }
  } else {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] &
          refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] -
            knots_neg[i])^p
        } else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] -
            knots_neg[i + 1])^(p - k)) / (factorial(k) *
            factorial(p - k)))
        }
      }
    }
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        } else if (refined_x_pos[j] <= knots_pos[i + 1] &
          refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1 / factorial(p)) * (refined_x_pos[j] -
            knots_pos[i])^p
        } else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] -
            knots_pos[i + 1])^(p - k)) / (factorial(k) *
            factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D
}


compute_weights_precision <- function(knots) {
  if (min(knots) >= 0) {
    as(diag(diff(knots)), "matrix")
  } else if (max(knots) < 0) {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    as(diag(diff(knots_neg)), "matrix")
  } else {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    d1 <- diff(knots_neg)
    d2 <- diff(knots_pos)
    Precweights1 <- diag(d1)
    Precweights2 <- diag(d2)
    methods::as(Matrix::bdiag(Precweights1, Precweights2), "matrix")
  }
}
