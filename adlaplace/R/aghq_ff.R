#' Build an \pkg{aghq} objective list from an adlaplace fit
#'
#' Returns a list suitable for \code{\link[aghq]{aghq}}. \code{fn}, \code{gr},
#' and \code{he} each call \code{\link{log_lik_laplace}} at the supplied outer
#' parameters, starting the inner solve at \code{fit$gamma} on every call.
#' They return the negated Laplace-approximate log posterior (same convention
#' as \code{\link{outer_fn}}), so call
#' \code{aghq::default_control(negate = TRUE)}.
#'
#' @param fit An \code{"adlaplace_fit"} from \code{\link{adlaplace}}.
#' @param control_inner Inner trust-region control list; defaults to
#'   \code{fit$control_inner} when present, otherwise a quiet short run.
#'
#' @return A list with
#' \describe{
#'   \item{fn}{Scalar negative log-likelihood from
#'     \code{log_lik_laplace(..., deriv = FALSE)}.}
#'   \item{gr}{Profile gradient from
#'     \code{log_lik_laplace(..., deriv = TRUE, return_hessians = FALSE)}.}
#'   \item{he}{Hessian of the profiled joint negative log-density in the outer
#'     parameters: the Schur complement of \code{hessian$outer} with respect
#'     to \code{gamma}, using \code{hessian$H_inv}. From
#'     \code{log_lik_laplace(..., deriv = TRUE, return_hessians = TRUE)}.}
#'   \item{mode_theta}{\code{fit$par_info$mle_internal}. \pkg{aghq} only calls
#'     \code{fn}, \code{gr}, and \code{he}, and keeps this element on
#'     \code{quad$optresults$ff}.}
#' }
#'
#' @seealso \code{\link{log_lik_laplace}}, \code{\link{adlaplace}}
#' @export
aghq_ff <- function(fit, control_inner = NULL) {
  if (!inherits(fit, "adlaplace_fit")) {
    stop("`fit` must be an adlaplace_fit object", call. = FALSE)
  }
  if (is.null(fit$ad_pack) || is.null(fit$config)) {
    stop("`fit` must contain ad_pack and config", call. = FALSE)
  }
  if (is.null(fit$gamma)) {
    stop("`fit` must contain gamma", call. = FALSE)
  }
  if (is.null(fit$par_info$mle_internal)) {
    stop("`fit$par_info$mle_internal` is missing", call. = FALSE)
  }
  if (is.null(control_inner)) {
    control_inner <- fit$control_inner
  }
  if (is.null(control_inner)) {
    control_inner <- list(maxit = 100L, report.level = 0, report.freq = 0)
  }

  ad_pack <- fit$ad_pack
  config <- fit$config
  gamma <- as.numeric(fit$gamma)
  sz <- ad_pack@sizes
  n_beta <- as.integer(sz["beta"])
  n_gamma <- as.integer(sz["gamma"])
  n_theta <- as.integer(sz["theta"])

  call_laplace <- function(theta, deriv, return_hessians) {
    log_lik_laplace(
      x = as.numeric(theta),
      gamma = gamma,
      ad_pack = ad_pack,
      config = config,
      control = control_inner,
      deriv = deriv,
      return_hessians = return_hessians
    )
  }

  list(
    fn = function(theta, ...) {
      as.numeric(call_laplace(theta, deriv = FALSE, return_hessians = TRUE)$neg_log_lik)
    },
    gr = function(theta, ...) {
      as.numeric(
        call_laplace(theta, deriv = TRUE, return_hessians = FALSE)$deriv$d_neg_log_lik
      )
    },
    he = function(theta, ...) {
      outer_profile_hessian(
        call_laplace(theta, deriv = TRUE, return_hessians = TRUE),
        n_beta = n_beta,
        n_gamma = n_gamma,
        n_theta = n_theta
      )
    },
    mode_theta = as.numeric(fit$par_info$mle_internal)
  )
}

#' Schur complement of the joint Hessian for the outer parameters
#'
#' \eqn{H_{\psi\psi} + H_{\psi\gamma} dU} with \eqn{dU = -H_{\gamma\gamma}^{-1} H_{\gamma\psi}},
#' the Hessian of the profiled joint negative log-density.
#'
#' @noRd
outer_profile_hessian <- function(res, n_beta, n_gamma, n_theta) {
  n_psi <- n_beta + n_theta
  if (n_psi < 1L) {
    return(matrix(numeric(), 0L, 0L))
  }
  dU <- res$gradient$outer$dU
  if (is.null(dU)) {
    stop(
      "log_lik_laplace did not return gradient$outer$dU; ",
      "call with deriv = TRUE and return_hessians = TRUE",
      call. = FALSE
    )
  }
  dU <- as.matrix(dU)
  if (nrow(dU) != n_gamma || ncol(dU) != n_psi) {
    stop(
      "gradient$outer$dU must be Ngamma by Nbeta + Ntheta",
      call. = FALSE
    )
  }
  H <- res$hessian$outer
  if (is.null(H)) {
    stop("log_lik_laplace did not return hessian$outer", call. = FALSE)
  }
  gamma_idx <- seq.int(n_beta + 1L, length.out = n_gamma)
  psi_idx <- c(
    seq_len(n_beta),
    seq.int(n_beta + n_gamma + 1L, length.out = n_theta)
  )
  Hpp <- as.matrix(H[psi_idx, psi_idx, drop = FALSE])
  Hpg <- as.matrix(H[psi_idx, gamma_idx, drop = FALSE])
  hess <- Hpp + Hpg %*% dU
  (hess + t(hess)) / 2
}
