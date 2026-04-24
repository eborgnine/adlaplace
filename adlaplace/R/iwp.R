#' @include 000.R
#' @include rpoly.R
NULL

setClass("iwp",
  slots = list(),
  contains = "model",
  prototype = prototype(
    # Default values for iwp-specific behavior
    by = character(0),
    type = factor("random", levels = .type_factor_levels)
  )
)

#' Integrated Wiener Process Term
#'
#' @description Creates and manages integrated Wiener process (IWP) model terms.
#' @name iwp-class
#' @docType class
#' @exportClass iwp
#'
#' @section Methods:
#' The following methods are available for `iwp` objects:
#' \describe{
#'   \item{\code{design(term, data)}}{Creates design matrix for IWP term}
#'   \item{\code{precision(term, data)}}{Creates precision matrix for IWP term}
#'   \item{\code{theta_info(term)}}{Extracts theta parameter information}
#'   \item{\code{beta_info(term, data)}}{Extracts beta parameter information}
#'   \item{\code{random_info(term, data)}}{Extracts random effects information}
#' }
NULL

#' Integrated Wiener Process Term Constructor
#'
#' @description Creates an integrated Wiener process (IWP) model term.
#' @name iwp
#' @param x Variable name.
#' @param p Order of the integrated Wiener process (default: 2).
#' @param ref_value Reference value for the basis.
#' @param knots Vector of knot locations.
#' @param range Range of the data (optional).
#' @param init Initial values for theta parameters.
#' @param lower Lower bounds for theta parameters.
#' @param upper Upper bounds for theta parameters.
#' @param parscale Parameter scales for optimization.
#' @param boundary_is_random Whether boundary should be treated as random.
#' @param include_poly Whether to include polynomial terms.
#' @return A list containing the `iwp` term object and optionally polynomial terms.
#' @examples
#' # Example usage:
#' # knots <- seq(0, 1, length.out = 5)
#' # iwp_term <- iwp(x = "age", knots = knots)
#' @export
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

  the_f <- stats::as.formula(paste0("~ 0 + ", x), env = new.env())
  result <- list()
  iwp_name <- paste("iwp", x, sep = "_")

  ref_value = ref_align(ref_value, knots)
  knots_ref = knots - ref_value

  result[[iwp_name]] <- methods::new("iwp",
    term = x,
    formula = the_f,
    p.order = as.integer(p),
    ref_value = ref_value,
    knots = knots_ref,
    init = init[1],
    lower = lower[1],
    upper = upper[1],
    parscale = parscale[1]
    # type is already set in prototype, no need to repeat
  )
  if(p==1) {
    include_poly = FALSE
  }
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
    methods::as(diag(diff(knots)), "matrix")
  } else if (max(knots) < 0) {
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    methods::as(diag(diff(knots_neg)), "matrix")
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

#' Align Reference Value to Nearest Knot
#'
#' @description Aligns a reference value to the nearest knot in a set of knots.
#' This is useful for ensuring that reference values used in polynomial terms
#' are aligned with the actual knot positions.
#'
#' @param ref_value A numeric value representing the desired reference point.
#' @param knots A numeric vector of knot positions.
#'
#' @return The knot value that is closest to the specified reference value.
#'
#' @examples
#' # Align reference value to nearest knot
#' knots <- c(10, 20, 30, 40, 50)
#' ref_align(23, knots)  # Returns 20
#'
#' @export
ref_align <- function(ref_value, knots) {
  knots[which.min(abs(knots - ref_value))]
}

# Method definitions for iwp class
#' @describeIn iwp-class Creates design matrix for IWP term
#' @param term An iwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("design", "iwp", function(term, data) {
  refined_x <- data[[term@term]] - term@ref_value

  if(any(refined_x < min(term@knots)) || any(refined_x > max(term@knots))) {
    warning("knots don't span the range of x ", term@term)
  }

  basis <- local_poly(term@knots, refined_x, term@p.order)
  result <- basis[, 1:ncol(basis), drop = FALSE]

  knots_string <- formatC(seq.int(ncol(result)),
    width = ceiling(log10(ncol(result))), flag = "0"
  )

  colnames(result) <- paste0(term@term, "_iwp_k", knots_string)
  result
})

#' @describeIn iwp-class Creates precision matrix for IWP term
#' @param term An iwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("precision", "iwp", function(term, data) {
  result <- Matrix::Matrix(compute_weights_precision(term@knots))

  knots_string <- formatC(seq.int(nrow(result)),
    width = ceiling(log10(nrow(result))), flag = "0"
  )

  dimnames(result) <- list(
    paste0(term@term, "_iwp_k", knots_string)
  )[c(1, 1)]
  result
})

#' @describeIn iwp-class Extracts theta parameter information for IWP term
#' @param term An iwp term object
#' @export
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

#' @describeIn iwp-class Extracts beta parameter information for IWP term
#' @param term An iwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("beta_info", "iwp", function(term, data) {
  # IWP terms don't have beta parameters
  return(NULL)
})

#' @describeIn iwp-class Extracts random effects information for IWP term
#' @param term An iwp term object
#' @param data A data frame containing the term variable
#' @export
setMethod("random_info", "iwp", function(term, data) {
  basis <- seq(1, length.out = length(term@knots) - 1)

  result <- expand.grid(
    term = term@term,
    model = "iwp",
    label = paste(c(term@term, "iwp"), collapse = "_"),
    by = NA, # iwp doesn't have hierarchical structure
    by_labels = NA,
    basis = basis,
    order = term@p.order,
    stringsAsFactors = FALSE
  )

  bnumPad <- formatC(result$basis,
    width = max(ceiling(c(1, log10(max(result$basis)))), na.rm = TRUE),
    flag = "0"
  )
  result$gamma_label <- paste0(result$label, "_k", bnumPad)

  result
})
