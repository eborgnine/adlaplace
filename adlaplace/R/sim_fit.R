#' @noRd
chol_inner_from_fit <- function(fit) {
  chol_inner <- fit$extra$hessian$chol_inner
  if (is.null(chol_inner) && !is.null(fit$hessian)) {
    chol_inner <- fit$hessian$chol_inner
  }
  if (is.null(chol_inner)) {
    stop(
      "fit must contain inner Cholesky factor in ",
      "$extra$hessian$chol_inner (from log_lik_laplace()) ",
      "or $hessian$chol_inner (from inner_opt())",
      call. = FALSE
    )
  }
  chol_inner
}

#' Simulate linear predictors on a new covariate grid
#'
#' Draws random effects from the Laplace approximation and evaluates the
#' linear predictor on \code{x} by summing fixed (\eqn{\beta}) and random
#' (\eqn{\gamma}) contributions from each model term.
#'
#' @param x Data frame of prediction covariates (same variables as the fitted model).
#' @param data Object returned by \code{\link{model_data}}.
#' @param fit Result from \code{\link{log_lik_laplace}} at the fitted outer parameters.
#' @param n Number of draws for random-effect simulation.
#'
#' @return A list with:
#' \describe{
#'   \item{eta}{Combined linear predictor, \code{nrow(x)} by \code{n} matrix
#'     (or vector if only fixed effects contribute).}
#'   \item{eta_fixed}{Fixed-effect part (\code{nrow(x)} by \code{1}).}
#'   \item{eta_random}{Random-effect part (\code{nrow(x)} by \code{n}).}
#'   \item{gamma_sims}{Draws from \code{\link{rmvnldl}} (\code{n} by \code{length(gamma)}).}
#'   \item{x}{The prediction grid \code{x} (echoed).}
#' }
#'
#' @export
sim_fit <- function(x, data, fit, n = 500L) {
  gamma_sims <- rmvnldl(
    n,
    mean = fit$opt$solution,
    chol_prec = chol_inner_from_fit(fit)
  )
  n_draws <- nrow(gamma_sims)
  n_grid <- nrow(x)

  if (is_model_data_bundle(data)) {
    info <- data$data$info
    term_list <- data$data$terms
  } else if (is(data, "ad_data") && is.list(data@precision) && !is.null(data@precision$info)) {
    info <- data@precision$info
    term_list <- list()
  } else {
    stop("`data` must be output from model_data()", call. = FALSE)
  }

  design_list <- lapply(
    term_list,
    function(term) {
      if (!any(term@term %in% names(x))) {
        return(NULL)
      }
      beta_here <- beta_info(term, data = x)
      if (!is.null(beta_here)) {
        if (any(beta_here$term %in% names(x))) {
          design_here <- design(term, data = x)
          par_id <- which(info$parameters$label == beta_here$label)
          return(
            as.matrix(
              design_here[, beta_here$beta_label, drop = FALSE] %*%
                fit$parameters[par_id]
            )
          )
        }
        return(NULL)
      }
      gamma_here <- random_info(term, data = x)
      if (!is.null(gamma_here)) {
        if (any(gamma_here$term %in% names(x))) {
          gamma_here_id <- match(
            gamma_here$gamma_label,
            info$gamma$gamma_label
          )
          if (any(is.na(gamma_here_id))) {
            warning("no match for labels: ", gamma_here$label[1])
            return(NULL)
          }
          design_here <- design(term, data = x)
          return(
            tcrossprod(
              design_here[, gamma_here$gamma_label, drop = FALSE],
              gamma_sims[, gamma_here_id, drop = FALSE]
            )
          )
        }
        return(NULL)
      }
      NULL
    }
  )
  design_list <- design_list[!vapply(design_list, is.null, logical(1))]

  eta_fixed <- matrix(0, nrow = n_grid, ncol = 1L)
  eta_random <- matrix(0, nrow = n_grid, ncol = n_draws)

  if (length(design_list)) {
    ncols <- vapply(design_list, ncol, integer(1))
    if (!all(ncols %in% c(1L, n_draws))) {
      warning("some design matrices have the wrong number of columns")
    }
    is_fixed <- ncols == 1L
    is_random <- ncols == n_draws
    if (any(is_fixed)) {
      eta_fixed <- Reduce(`+`, design_list[is_fixed])
      if (!is.matrix(eta_fixed)) {
        eta_fixed <- matrix(eta_fixed, ncol = 1L)
      }
    }
    if (any(is_random)) {
      eta_random <- Reduce(`+`, design_list[is_random])
    }
  }

  eta <- as.vector(eta_fixed) + eta_random

  eta
}
