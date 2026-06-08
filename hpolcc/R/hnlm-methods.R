#' @keywords internal
hnlm_check_fitted <- function(object) {
  if ("hnlm_dev" %in% class(object)) {
    stop("requires a fitted hnlm object", call. = FALSE)
  }
  invisible(object)
}

#' @keywords internal
hpolcc_pinv <- function(M, tol = NULL) {
  s <- svd(M)
  if (is.null(tol)) {
    tol <- max(dim(M)) * max(s$d) * .Machine$double.eps
  }
  d <- s$d
  d[d <= tol] <- 0
  positive <- d > 0
  d[positive] <- 1 / d[positive]
  s$v %*% diag(d, nrow = length(d), ncol = length(d)) %*% t(s$u)
}

#' @keywords internal
hnlm_outer_vcov <- function(object) {
  hnlm_check_fitted(object)
  H <- object$hessian$outer
  if (is.null(H) || inherits(H, "try-error")) {
    stop("outer Hessian not available for standard errors", call. = FALSE)
  }
  H <- as.matrix(H)
  par <- object$optim$par
  if (length(par) != nrow(H) || length(par) != ncol(H)) {
    stop(
      "outer Hessian dimension (", nrow(H), " x ", ncol(H),
      ") does not match length(optim$par) (", length(par), ")",
      call. = FALSE
    )
  }
  V <- tryCatch(
    solve(H),
    error = function(e) {
      warning(
        "outer Hessian is singular; using Moore-Penrose inverse for standard errors",
        call. = FALSE
      )
      hpolcc_pinv(H)
    }
  )
  dimnames(V) <- list(names(par), names(par))
  V
}

#' Variance-covariance matrix of outer parameters
#'
#' Inverse of the outer Hessian from profile likelihood optimization
#' (\code{hessian$outer}). Falls back to a Moore-Penrose inverse with a
#' warning if the Hessian is singular.
#'
#' @param object A fitted \code{hnlm} object.
#' @param ... Not used.
#' @return Variance-covariance matrix on the optimization scale of
#'   \code{optim$par}.
#' @export
vcov.hnlm <- function(object, ...) {
  hnlm_outer_vcov(object)
}

#' Summary method for \code{hnlm} fits
#'
#' Wald-style standard errors from \code{solve(hessian$outer)} (or a
#' Moore-Penrose inverse if singular). Estimates are on the optimization
#' scale (\code{optim$par}); transformed theta parameters are on the log scale.
#'
#' @param object A fitted \code{hnlm} object.
#' @param ... Not used.
#' @return An object of class \code{summary.hnlm} with coefficient table,
#'   variance matrix, log-likelihood, and convergence flag.
#' @export
summary.hnlm <- function(object, ...) {
  hnlm_check_fitted(object)
  V <- hnlm_outer_vcov(object)
  est <- object$optim$par
  se <- sqrt(diag(V))
  z <- est / se
  pval <- 2 * stats::pnorm(-abs(z))
  labels <- object$model_data$data$info$parameters$label
  coef_table <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = pval
  )
  if (length(labels) == nrow(coef_table)) {
    rownames(coef_table) <- labels
  }
  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      vcov = V,
      log_lik = object$log_lik,
      df = length(est),
      nobs = nrow(object$model_data$data$data),
      converged = object$converged,
      scale = "outer"
    ),
    class = "summary.hnlm"
  )
}

#' @rdname summary.hnlm
#' @param x A \code{summary.hnlm} object.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @export
print.summary.hnlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat("Log-likelihood:", format(x$log_lik, digits = digits), "\n")
  cat("Observations:", x$nobs, "\n")
  cat("Converged:", if (isTRUE(x$converged)) "yes" else "no", "\n")
  cat("\nOuter parameters (optimization scale):\n")
  print(
    signif(x$coefficients, digits = digits),
    print.gap = 2L
  )
  cat("\nStandard errors from inverse outer Hessian.\n")
  cat("Transformed theta parameters are on the log scale.\n")
  invisible(x)
}

#' @export
print.hnlm <- function(x, ...) {
  is_dev <- "hnlm_dev" %in% class(x)
  if (is_dev) {
    cat("hnlm development bundle (not fitted)\n")
  } else {
    cat("hnlm fit\n")
  }
  if (!is.null(x$formula)) {
    cat("  formula:", paste(deparse(x$formula), collapse = " "), "\n")
  }
  if (!is.null(x$model_data$data$data)) {
    cat("  n =", nrow(x$model_data$data$data), "\n")
  }
  if (!is_dev && !is.null(x$log_lik)) {
    cat("  log-lik:", format(x$log_lik, digits = 5), "\n")
  }
  if (!is_dev && !is.null(x$optim)) {
    conv <- if (isTRUE(x$converged)) "yes" else "no"
    cat("  converged:", conv, "\n")
  }
  if (!is_dev && !is.null(x$coefficients$beta) && nrow(x$coefficients$beta) > 0L) {
    cat("  beta:\n")
    print(x$coefficients$beta[, c("label", "mle"), drop = FALSE], row.names = FALSE)
  }
  if (!is_dev && !is.null(x$coefficients$theta) && nrow(x$coefficients$theta) > 0L) {
    cat("  theta:\n")
    print(x$coefficients$theta[, c("label", "mle"), drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

#' @export
coef.hnlm <- function(object, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("coef() requires a fitted hnlm object", call. = FALSE)
  }
  beta <- object$coefficients$beta
  theta <- object$coefficients$theta
  out <- numeric(0)
  if (nrow(beta) > 0L) {
    out <- stats::setNames(beta$mle, beta$beta_label)
  }
  if (nrow(theta) > 0L) {
    out <- c(out, stats::setNames(theta$mle, theta$label))
  }
  out
}

#' @export
logLik.hnlm <- function(object, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("logLik() requires a fitted hnlm object", call. = FALSE)
  }
  nobs <- nrow(object$model_data$data$data)
  structure(
    object$log_lik,
    df = length(coef(object)),
    nobs = nobs,
    class = "logLik"
  )
}

#' @export
simulate.hnlm <- function(object, nsim = 500, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("simulate() requires a fitted hnlm object", call. = FALSE)
  }
  adlaplaceHgp::cond_sim_iwp(
    laplace = object$laplace,
    model_data = object$model_data,
    n = nsim,
    ...
  )
}
