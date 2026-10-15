# Verbatim copies of Ziang Zhang / BayesGP near-monotone reference helpers
# for regression tests only. Source:
#   https://github.com/Bayes-GP/BayesGP (branch bayesgp-near-monotone-usability)
#   R/05_near_monotone_fem.R, R/04_near_monotone_state_space.R, R/07_near_monotone_helpers.R
# Not part of the adlaplaceFem runtime API. FEM tests skip without Suggests: fda.

## let a(.) be a given function
Compute_Prec <- function(a = NULL, alpha = NULL, c, k, region, accuracy = 0.01, boundary = TRUE){
  a <- resolve_curvature_parameter(a = a, alpha = alpha)
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  x <- fem_integration_grid(region = region, accuracy = accuracy, c = c)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
  B1matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE))
  B2matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE))
  a_func <- function(x) {-1/(a*(x+c))}
  a_matrix <- a_func(x)
  B1a <- Matrix::Matrix(apply(B1matrix, 2, function(x) x * a_matrix), sparse = TRUE)
  B2matrix <- Matrix::Matrix(B2matrix, sparse = TRUE)

  ### Compute I, L, T:
  Numerical_I <- Matrix::Diagonal(x = c(diff(c(0, x))))

  G <- Matrix::t(B1a) %*% Numerical_I %*% B1a
  C <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(Matrix::t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}

fem_integration_grid <- function(region, accuracy, c){
  region <- range(region)
  lower <- min(region)
  upper <- max(region)

  if(isTRUE(all.equal(c, 0)) && isTRUE(all.equal(lower, 0))){
    x <- seq(lower + accuracy, upper, by = accuracy)
  } else {
    x <- seq(lower, upper, by = accuracy)
  }
  if(length(x) == 0 || tail(x, 1) < upper){
    x <- c(x, upper)
  }
  x
}

Compute_Design <- function(x, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
}

Compute_Prec_rev <- function(a = NULL, alpha = NULL, c, k, region, accuracy = 0.01, boundary = TRUE){
  a <- resolve_curvature_parameter(a = a, alpha = alpha)
  ss <- function(M) {Matrix::forceSymmetric(M + Matrix::t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE))
  B1matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE))
  B2matrix <- as.matrix(fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE))
  a_func <- function(x) {1/(a*(c-x))}
  a_matrix <- a_func(x)
  B1a <- Matrix::Matrix(apply(B1matrix, 2, function(x) x * a_matrix), sparse = TRUE)
  B2matrix <- Matrix::Matrix(B2matrix, sparse = TRUE)

  ### Compute I, L, T:
  Numerical_I <- Matrix::Diagonal(x = c(diff(c(0, x))))

  G <- Matrix::t(B1a) %*% Numerical_I %*% B1a
  C <- Matrix::t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- ss(Matrix::t(B1a) %*% Numerical_I %*% B2matrix)
  Q <- G - M + C
  Matrix::forceSymmetric(Q)
}

SS_cov_c1 <- function(s,t,c = 1, alpha = 2){
  -1/3*(3*(2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*c^2*t + 3*(2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*c*t^2 + (2*alpha^4 + 5*alpha^3 + 2*alpha^2)*(c + s)^(2/alpha)*t^3 + (9*alpha^3*c^3 - 6*(alpha^4 - 2*alpha^3 + alpha^2)*c^2*s - 6*(alpha^4 - 2*alpha^3 + alpha^2)*c*s^2 - 2*(alpha^4 - 2*alpha^3 + alpha^2)*s^3)*(c + s)^(2/alpha) + 3*((2*alpha^4 + alpha^3)*c^3 + 2*(2*alpha^4 + alpha^3)*c^2*s + (2*alpha^4 + alpha^3)*c*s^2 + ((2*alpha^4 + alpha^3)*c^2 + 2*(2*alpha^4 + alpha^3)*c*s + (2*alpha^4 + alpha^3)*s^2)*t)*(c + t)^(2/alpha) - 6*(((alpha^4 + 2*alpha^3)*c + (alpha^4 + 2*alpha^3)*s)*(c + s)^(1/alpha)*t^2 + 2*((alpha^4 + 2*alpha^3)*c^2 + (alpha^4 + 2*alpha^3)*c*s)*(c + s)^(1/alpha)*t + ((alpha^4 + 2*alpha^3)*c^3 + (alpha^4 + 2*alpha^3)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))/((2*alpha^4 + alpha^3 - 6*alpha^2 + alpha + 2)*(c + s)^(2/alpha))
}

SS_cov_c2 <- function(s,t,c = 1){
  -1/3*c^3*log(c + s)^2 - 1/27*(9*log(c + s)^2 + 6*log(c + s) + 2)*t^3 - 2/9*c^3*log(c + s) + 2/9*c^2*s + 2/9*c*s^2 + 2/27*s^3 - 1/9*(9*c*log(c + s)^2 + 6*c*log(c + s) + 2*c)*t^2 - 1/3*(c^3 + 3*c^2*t + 3*c*t^2 + t^3)*log(c + t)^2 - 1/9*(9*c^2*log(c + s)^2 + 6*c^2*log(c + s) + 2*c^2)*t + 2/9*(t^3*(3*log(c + s) + 1) + 3*c^3*log(c + s) + c^3 + 3*(3*c*log(c + s) + c)*t^2 + 3*(3*c^2*log(c + s) + c^2)*t)*log(c + t)

}

SS_cov_c3 <- function(s,t,c = 1, alpha = -3){
  alpha <- abs(alpha)
  -1/3*(3*((2*alpha^4 - alpha^3)*c^2 + 2*(2*alpha^4 - alpha^3)*c*s + (2*alpha^4 - alpha^3)*s^2)*(c + s)^(2/alpha)*t + 3*((2*alpha^4 - alpha^3)*c^3 + 2*(2*alpha^4 - alpha^3)*c^2*s + (2*alpha^4 - alpha^3)*c*s^2)*(c + s)^(2/alpha) - (9*alpha^3*c^3 + 6*(alpha^4 + 2*alpha^3 + alpha^2)*c^2*s + 6*(alpha^4 + 2*alpha^3 + alpha^2)*c*s^2 + 2*(alpha^4 + 2*alpha^3 + alpha^2)*s^3 - 3*(2*alpha^4 - 5*alpha^3 + 2*alpha^2)*c^2*t - 3*(2*alpha^4 - 5*alpha^3 + 2*alpha^2)*c*t^2 - (2*alpha^4 - 5*alpha^3 + 2*alpha^2)*t^3)*(c + t)^(2/alpha) - 6*(((alpha^4 - 2*alpha^3)*c + (alpha^4 - 2*alpha^3)*s)*(c + s)^(1/alpha)*t^2 + 2*((alpha^4 - 2*alpha^3)*c^2 + (alpha^4 - 2*alpha^3)*c*s)*(c + s)^(1/alpha)*t + ((alpha^4 - 2*alpha^3)*c^3 + (alpha^4 - 2*alpha^3)*c^2*s)*(c + s)^(1/alpha))*(c + t)^(1/alpha))/((2*alpha^4 - alpha^3 - 6*alpha^2 - alpha + 2)*(c + t)^(2/alpha))
}

SS_cov_c4 <- function(s,t,c = 1){
  1/12*(4*c*s^3 + 3*s^4 - 4*c*t^3 - t^4 + 6*(2*c*s + s^2)*t^2 - 4*(3*c*s^2 + 2*s^3)*t)/(c + t)
}

SS_cov_c5 <- function(s,t,c = 1){
  -16/27*c^3 - 4/3*c^2*s - 4/3*c*s^2 - 4/9*s^3 - 4/9*c^2*t - 4/9*c*t^2 - 4/27*t^3 + 16/27*(c^2 + c*s + (c + s)*t)*sqrt(c + s)*sqrt(c + t) + 4/9*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + s) - 4/9*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + t)
}

SS_cov_c6 <- function(s,t,c = 1){
  1/27*(6*c^5*s + 15*c^4*s^2 + 20*c^3*s^3 + 15*c^2*s^4 + 6*c*s^5 + s^6 - 15*c^2*t^4 - 6*c*t^5 - t^6 - 2*(10*c^3 + 3*(c^3 + 3*c^2*s + 3*c*s^2 + s^3)*log(c + s))*t^3 - 3*(5*c^4 + 6*(c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*log(c + s))*t^2 - 6*(c^5 + 3*(c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*log(c + s))*t - 6*(c^6 + 3*c^5*s + 3*c^4*s^2 + c^3*s^3)*log(c + s) + 6*(c^6 + 3*c^5*s + 3*c^4*s^2 + c^3*s^3 + (c^3 + 3*c^2*s + 3*c*s^2 + s^3)*t^3 + 3*(c^4 + 3*c^3*s + 3*c^2*s^2 + c*s^3)*t^2 + 3*(c^5 + 3*c^4*s + 3*c^3*s^2 + c^2*s^3)*t)*log(c + t))/(c^3 + 3*c^2*t + 3*c*t^2 + t^3)
}

resolve_curvature_parameter <- function(a = NULL, alpha = NULL){
  if(!is.null(a) && !is.null(alpha) && !isTRUE(all.equal(a, alpha))){
    stop("Please supply only one curvature parameter or keep `a` and `alpha` identical.")
  }

  curvature <- if(!is.null(a)) a else alpha

  if(is.null(curvature)){
    stop("A curvature parameter must be supplied via `a` or `alpha`.")
  }

  curvature
}

SS_cov <- function(s,t,c = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(alpha != 1 & alpha > 0){
    SS_cov_c1(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == 1){
    SS_cov_c2(s = s, t = t, c = c)
  }
  else if(alpha < 0 & alpha != -1 & alpha != -2 & alpha != -(1/2)){
    SS_cov_c3(s = s, t = t, c = c, alpha = alpha)
  }
  else if(alpha == -1){
    SS_cov_c4(s = s, t = t, c = c)
  }
  else if(alpha == -2){
    SS_cov_c5(s = s, t = t, c = c)
  }
  else if(alpha == -1/2){
    SS_cov_c6(s = s, t = t, c = c)
  }
  else{
    message("Error: alpha value is not supported.")
  }
}

PSD_compute_mgp <- function(x, h, c = 1, sd = 1, a = NULL, alpha = NULL){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  abs(sd) * sqrt(SS_cov(s = (x + h), t = (x), c = c, alpha = alpha))
}

prior_conversion_mgp <- function(prior, d, h = NULL, c = 0.1, a = NULL, alpha = NULL, x = 0){
  alpha <- resolve_curvature_parameter(a = a, alpha = alpha)
  if(!is.null(h)){
    d <- h
  }
  correction_factor <- PSD_compute_mgp(x = x, h = d, c = c, a = alpha, sd = 1)
  prior_q <- list(alpha = prior$alpha, u = (prior$u/correction_factor))
  prior_q
}

resolve_reference_location <- function(x, initial_location){
  domain_range <- range(x)

  if(is.null(initial_location)){
    return(domain_range[1])
  }

  if(is.character(initial_location)){
    if(length(initial_location) != 1){
      stop("`initial_location` must be a scalar.")
    }
    initial_location <- tolower(initial_location)
    if(initial_location == "left"){
      return(domain_range[1])
    }
    if(initial_location == "middle"){
      return(mean(domain_range))
    }
    if(initial_location == "right"){
      return(domain_range[2])
    }
    stop("`initial_location` must be numeric or one of 'left', 'middle', or 'right'.")
  }

  if(!is.numeric(initial_location) || length(initial_location) != 1 || !is.finite(initial_location)){
    stop("`initial_location` must be numeric or one of 'left', 'middle', or 'right'.")
  }

  tolerance <- sqrt(.Machine$double.eps)
  if(initial_location < domain_range[1] - tolerance || initial_location > domain_range[2] + tolerance){
    stop("`initial_location` must lie inside the near-monotone domain.")
  }

  as.numeric(initial_location)
}

make_boxcox_transform <- function(a = NULL, alpha = NULL, c = 1,
                                  normalized = TRUE, ref_location = 0){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  lambda <- (curvature - 1) / curvature
  reference <- c + ref_location

  if(normalized){
    if(reference <= 0){
      stop("Normalized Box-Cox transforms require `c + ref_location > 0`.")
    }
    if(isTRUE(all.equal(lambda, 0))){
      return(function(x){
        reference * (log(x + c) - log(reference))
      })
    }

    return(function(x){
      ((x + c)^lambda - reference^lambda) / (lambda * reference^(lambda - 1))
    })
  }

  if(isTRUE(all.equal(lambda, 0))){
    return(function(x){
      log(x + c)
    })
  }

  function(x){
    (x + c)^lambda
  }
}

boundary_basis_matrix <- function(x, a = NULL, alpha = NULL, c = 1,
                                  initial_location = min(x),
                                  normalized = TRUE){
  curvature <- resolve_curvature_parameter(a = a, alpha = alpha)
  shifted_x <- x - initial_location
  transformed_arg <- shifted_x + c

  if(any(transformed_arg <= 0)){
    stop("All evaluation points must satisfy `x + c > 0` on the original scale.")
  }

  transform <- make_boxcox_transform(
    a = curvature,
    c = c,
    normalized = normalized,
    ref_location = 0
  )

  matrix(transform(shifted_x), ncol = 1)
}

