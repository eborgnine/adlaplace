#' @keywords internal
hnlm_check_fitted <- function(object) {
  if ("hnlm_dev" %in% class(object)) {
    stop("requires a fitted hnlm object", call. = FALSE)
  }
  invisible(object)
}

#' @keywords internal
hnlm_nobs <- function(object) {
  shards <- object$call$config$shards
  if (!is.null(shards)) {
    return(as.integer(nrow(shards)))
  }
  NA_integer_
}

#' @keywords internal
hnlm_model_data_stub <- function(object) {
  list(
    terms = object$call$terms,
    data = list(info = object$info)
  )
}

#' @keywords internal
hnlm_attach_parameter_se <- function(coefficients, hessian_outer, vcov = NULL) {
  if (inherits(coefficients, "try-error") ||
      is.null(coefficients$parameters) ||
      inherits(hessian_outer, "try-error") ||
      is.null(hessian_outer)) {
    return(coefficients)
  }
  H <- as.matrix(hessian_outer)
  n_par <- nrow(coefficients$parameters)
  if (nrow(H) != n_par || ncol(H) != n_par) {
    return(coefficients)
  }
  if (!is.null(vcov) && all(dim(vcov) == n_par)) {
    # vcov already computed (e.g. by adlaplace()); skip re-inverting.
    V <- as.matrix(vcov)
  } else {
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
  }
  se_opt <- sqrt(pmax(0, diag(V)))
  params <- coefficients$parameters
  se <- se_opt
  if ("transform" %in% names(params)) {
    idx <- which(params$transform %in% TRUE)
    if (length(idx) > 0L) {
      se[idx] <- params$mle[idx] * se_opt[idx]
    }
  }
  coefficients$parameters$se <- se
  coefficients
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
#' Returns \code{coefficients$parameters} from the fit, including \code{mle} and
#' \code{se} (when the outer Hessian was available at fit time).
#'
#' @param object A fitted \code{hnlm} object.
#' @param ... Not used.
#' @return An object of class \code{summary.hnlm} with the parameters table,
#'   variance matrix, log-likelihood, and convergence flag.
#' @export
summary.hnlm <- function(object, ...) {
  hnlm_check_fitted(object)
  params <- object$coefficients$parameters
  if (is.null(params) || nrow(params) == 0L) {
    stop("coefficients$parameters not available", call. = FALSE)
  }
  V <- tryCatch(
    hnlm_outer_vcov(object),
    error = function(e) NULL
  )
  structure(
    list(
      call = object$call$call,
      coefficients = params,
      vcov = V,
      log_lik = object$log_lik,
      df = length(object$optim$par),
      nobs = hnlm_nobs(object),
      converged = object$converged
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
  cat("\nParameters:\n")
  cols <- intersect(c("label", "term", "model", "mle", "se"), names(x$coefficients))
  tab <- x$coefficients[, cols, drop = FALSE]
  num_cols <- vapply(tab, is.numeric, logical(1))
  tab[num_cols] <- signif(tab[num_cols], digits = digits)
  print(tab, row.names = FALSE)
  if (!"se" %in% names(x$coefficients)) {
    cat("\nStandard errors not available (outer Hessian missing or invalid).\n")
  }
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
  if ("hnlm_dev" %in% class(x)) {
    if (!is.null(x$formula)) {
      cat("  formula:", paste(deparse(x$formula), collapse = " "), "\n")
    }
    if (!is.null(x$model_data$data$data)) {
      cat("  n =", nrow(x$model_data$data$data), "\n")
    }
  } else if (!is.null(x$call$call$formula)) {
    cat("  formula:", paste(deparse(x$call$call$formula), collapse = " "), "\n")
    nobs <- hnlm_nobs(x)
    if (!is.na(nobs)) {
      cat("  n =", nobs, "\n")
    }
  }
  if (!is_dev && !is.null(x$log_lik)) {
    cat("  log-lik:", format(x$log_lik, digits = 5), "\n")
  }
  if (!is_dev && !is.null(x$optim)) {
    conv <- if (isTRUE(x$converged)) "yes" else "no"
    cat("  converged:", conv, "\n")
  }
  if (!is_dev && length(x$coefficients) > 0L &&
      !is.null(x$coefficients$parameters) &&
      nrow(x$coefficients$parameters) > 0L) {
    params <- x$coefficients$parameters
    n_beta <- nrow(x$info$beta)
    n_theta <- nrow(x$info$theta)
    if (n_beta > 0L) {
      cat("  beta:\n")
      print(
        params[seq_len(n_beta), c("label", "mle"), drop = FALSE],
        row.names = FALSE
      )
    }
    if (n_theta > 0L) {
      cat("  theta:\n")
      print(
        params[seq(n_beta + 1L, length.out = n_theta), c("label", "mle"), drop = FALSE],
        row.names = FALSE
      )
    }
  }
  invisible(x)
}

#' @export
coef.hnlm <- function(object, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("coef() requires a fitted hnlm object", call. = FALSE)
  }
  params <- object$coefficients$parameters
  if (is.null(params) || nrow(params) == 0L) {
    return(numeric(0))
  }
  n_beta <- nrow(object$info$beta)
  n_theta <- nrow(object$info$theta)
  out <- numeric(0)
  if (n_beta > 0L) {
    out <- stats::setNames(
      params$mle[seq_len(n_beta)],
      object$info$beta$beta_label
    )
  }
  if (n_theta > 0L) {
    out <- c(
      out,
      stats::setNames(
        params$mle[seq(n_beta + 1L, length.out = n_theta)],
        object$info$theta$label
      )
    )
  }
  out
}

#' @export
logLik.hnlm <- function(object, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("logLik() requires a fitted hnlm object", call. = FALSE)
  }
  structure(
    object$log_lik,
    df = length(stats::coef(object)),
    nobs = hnlm_nobs(object),
    class = "logLik"
  )
}

#' @export
simulate.hnlm <- function(object, nsim = 500, ...) {
  if ("hnlm_dev" %in% class(object)) {
    stop("simulate() requires a fitted hnlm object", call. = FALSE)
  }
  adlaplaceHgp::cond_sim_iwp(
    fit = object$extra,
    model_data = hnlm_model_data_stub(object),
    n = nsim,
    ...
  )
}
