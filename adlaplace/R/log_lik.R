#' Outer objective and gradient wrappers
#'
#' Convenience wrappers around \code{\link{log_lik_laplace}} for use with outer optimizers
#' that expect \code{fn} / \code{gr} callbacks. Both functions solve (or warm-start)
#' the inner problem over \code{gamma} via \code{\link{log_lik_laplace}} and update a mutable
#' \code{cache} environment with the latest inner solution.
#'
#' \describe{
#' \item{\code{outer_fn()}}{Returns the scalar objective negative log likelihood.}
#' \item{\code{outer_gr()}}{Returns the gradient.}
#' }
#'
#' @param x Numeric outer parameter vector \code{c(beta, theta)}.
#' @param config Configuration list passed to \code{\link{log_lik_laplace}}.
#' @param ad_pack \code{ad_pack} object from \code{\link{ad_pack}}.
#' @param ... Additional arguments forwarded to \code{\link{log_lik_laplace}}
#'   (for example \code{data} or \code{package}).
#' @param control_inner A list of control options forwarded to the \code{control}
#'   argument of \code{\link{log_lik_laplace}} for the inner optimization.
#' @param cache An \code{\link[base]{environment}} used to warm-start and store the inner
#'   \code{gamma} solution. If \code{cache$gamma} is missing or the wrong length,
#'   it is initialized from \code{config$gamma} (if present) or zeros of length
#'   \code{ad_pack@sizes["gamma"]}. Both functions update \code{cache$gamma} after
#'   each evaluation.
#'
#' @return
#' \itemize{
#' \item \code{outer_fn}: a numeric scalar \code{neg_log_lik}.
#' \item \code{outer_gr}: gradient of \code{neg_log_lik} (same sign as \code{outer_fn}).
#' }
#'
#' @seealso \code{\link{log_lik_laplace}}
#'
#' @examples
#' \dontrun{
#' cache <- new.env(parent = emptyenv())
#' cache$gamma <- rep(0, nrow(data$ATp))
#' ad_pack <- ad_pack(data, config)
#'
#' val <- outer_fn(x = x0, data = data, config = config, cache = cache, ad_pack = ad_pack)
#' gr <- outer_gr(x = x0, data = data, config = config, cache = cache, ad_pack = ad_pack)
#' }
#'
#' @name outer_optim_wrappers
#' @rdname outer_optim_wrappers
#' @export
outer_fn <- function(
  x, config, cache, ad_pack, control_inner = list(), ...
) {
  assign("last_par_fn", x, cache)
  num_gamma <- as.integer(ad_pack@sizes["gamma"])
  cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

  if (isTRUE(config$verbose)) {
    message("outer_fn: calling inner_opt(deriv=FALSE)...")
    utils::flush.console()
  }
  result <- adlaplace::inner_opt(
    parameters = x,
    gamma = cache$gamma,
    ad_pack = ad_pack,
    control = control_inner,
    deriv = FALSE,
    verbose = isTRUE(config$verbose)
  )
  if (isTRUE(config$verbose)) {
    message("outer_fn: inner_opt done neg_log_lik=", format(result$neg_log_lik))
    utils::flush.console()
  }

  assign("gamma", result$inner_opt$solution, cache)
  result$neg_log_lik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr <- function(
  x, config, cache, ad_pack, control_inner = list(), ...
) {
  assign("last_par_gr", x, cache)
  num_gamma <- as.integer(ad_pack@sizes["gamma"])
  cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

  if (isTRUE(config$verbose)) {
    message("outer_gr: calling log_lik_laplace(deriv=TRUE)...")
    utils::flush.console()
  }
  result <- adlaplace::log_lik_laplace(
    x = x, config = config,
    gamma = cache$gamma,
    control = control_inner,
    ad_pack = ad_pack,
    deriv = TRUE, ...
  )
  if (isTRUE(config$verbose)) {
    message(
      "outer_gr: done length(d_neg_log_lik)=",
      length(result$deriv$d_neg_log_lik)
    )
    utils::flush.console()
  }
  assign("gamma", result$inner_opt$solution, cache)
  result$deriv$d_neg_log_lik
}

#' Log-likelihood with inner Laplace optimization
#'
#' Evaluates the (profiled) log-likelihood for a hierarchical model by solving an
#' inner optimization problem for the latent vector \code{gamma} (e.g. random
#' effects) given outer parameters \code{x} (typically \code{beta} and \code{theta}).
#' The function delegates the inner optimization to \code{inner_opt()} from the
#' selected backend package (by default \pkg{adlaplace}).
#'
#' @param x Numeric vector of outer parameters (\code{beta}, then \code{theta}).
#'   Length must equal \code{num_beta + num_theta} from \code{ad_pack@sizes}.
#' @param config A list of backend options (\code{verbose}, \code{package}, etc.).
#'   When \code{ad_pack} is missing, must also include \code{beta}, \code{gamma},
#'   and \code{theta} so \code{ad_pack()} can be built from \code{data}.
#' @param gamma Optional numeric vector of starting values for the inner
#'   parameter \code{gamma}. If missing, uses \code{config$gamma} when present,
#'   otherwise a zero vector of length \code{ad_pack@sizes["gamma"]}.
#' @param control List of control parameters passed *as-is* to the backend inner
#'   optimizer (e.g., \code{report.level}, \code{report.freq}). See backend
#'   documentation (e.g., \pkg{trustOptim}) for supported options.
#' @param ad_pack Optional \code{ad_pack} object from \code{\link{ad_pack}}.
#'   This is a single backend handle (no separate inner/outer handles). If
#'   missing, it will be constructed automatically using \code{data}.
#' @param data Optional data list used to build \code{ad_pack} when \code{ad_pack}
#'   is not supplied.
#' @param package Character scalar naming the backend package to use for
#'   \code{ad_pack()} and \code{inner_opt()}. Defaults to \code{"adlaplace"}.
#'   Other backends must export a compatible \code{ad_pack()} builder and \code{inner_opt()}.
#' @param deriv Logical scalar. If \code{TRUE}, include derivative quantities in
#'   the output (gradient, intermediate derivatives).
#'
#' @details
#' The default \pkg{adlaplace} backend uses a single AD handle. This function
#' passes that handle to \code{inner_opt(..., ad_pack = ad_pack)}. Parameter block
#' sizes are read from \code{ad_pack@sizes}; \code{config} is not required to
#' contain \code{beta}, \code{gamma}, or \code{theta} when \code{ad_pack} is
#' supplied.
#'
#' When \code{deriv=FALSE}, the return value mirrors \code{\link{inner_opt}()}
#' with \code{inner_opt}, top-level \code{gradient}, and \code{hessian}.
#' When \code{deriv=TRUE}, the same list is augmented with \code{deriv} and
#' profile-derivative pieces in \code{gradient$outer}.
#'
#' @return A list. With \code{deriv=FALSE}, same structure as
#' \code{\link{inner_opt}} (without \code{parameters}). With \code{deriv=TRUE},
#' additionally:
#' \describe{
#'   \item{deriv}{Data frame of profile derivative pieces.}
#'   \item{gradient$outer$dU, gradient$outer$trace3}{Intermediate profile objects.}
#' }
#'
#' @note The Laplace approximation assumes the inner Hessian is positive definite
#'   at the optimum. If \code{\link{inner_opt}} fails to converge or returns a
#'   non-invertible Hessian, \code{log_lik_laplace()} issues a warning or error.
#'
#' @seealso
#' \code{\link{ad_pack}}, \code{\link{inner_opt}}
#'
#' @export
log_lik_laplace <- function(
  x,
  config = list(),
  gamma,
  control = list(report.level = 4, report.freq = 1),
  ad_pack,
  data,
  package = c(config[["package"]], "adlaplace")[1L],
  deriv = FALSE
) {
  if (missing(ad_pack)) {
    if (missing(data)) {
      stop("at least one of data and ad_pack must be supplied", call. = FALSE)
    }
    ad_pack <- adlaplace::ad_pack(data, config)
  }

  if (!is(ad_pack, "ad_pack")) {
    stop("ad_pack must be an ad_pack object", call. = FALSE)
  }
  sz <- ad_pack@sizes
  if (length(sz) < 3L || any(is.na(sz[c("beta", "gamma", "theta")]))) {
    stop(
      "ad_pack@sizes must contain finite beta, gamma, and theta",
      call. = FALSE
    )
  }
  num_beta <- as.integer(sz["beta"])
  num_theta <- as.integer(sz["theta"])
  num_gamma <- as.integer(sz["gamma"])

  n_outer <- length(x)
  n_outer_expected <- num_beta + num_theta
  if (n_outer != n_outer_expected) {
    stop(
      "length(x) (", n_outer, ") must equal num_beta + num_theta (",
      n_outer_expected, ")",
      call. = FALSE
    )
  }

  gamma_use <- if (!missing(gamma)) {
    gamma
  } else if (!is.null(config[["gamma"]]) && length(config[["gamma"]]) > 0L) {
    config[["gamma"]]
  } else {
    rep(0, num_gamma)
  }
  if (length(gamma_use) != num_gamma) {
    stop(
      "length(gamma) (", length(gamma_use), ") must equal num_gamma (",
      num_gamma, ")",
      call. = FALSE
    )
  }

  result_inner <- inner_opt(
    parameters = x,
    gamma = gamma_use,
    control = control,
    ad_pack = ad_pack,
    deriv = deriv,
    verbose = isTRUE(config[["verbose"]])
  )
  result <- restructure_laplace_result(result_inner, control_inner = control)

  if (!deriv) {
    return(result)
  }

  the_deriv <- log_lik_deriv(
    full_parameters = result$full_parameters,
    hessian_pack = result$hessian,
    grad = result$gradient$outer$grad,
    ad_pack = ad_pack,
    verbose = isTRUE(config[["verbose"]])
  )

  if (!is.null(the_deriv$extra$half_H_inv)) {
    result$hessian$half_H_inv <- the_deriv$extra$half_H_inv
  }
  if (!is.null(the_deriv$extra$H_inv)) {
    result$hessian$H_inv <- the_deriv$extra$H_inv
  }
  result$gradient$outer$trace3 <- the_deriv$extra$trace3
  result$gradient$outer$dU <- the_deriv$extra$dU
  if (!is.null(result$hessian$trace3)) {
    result$hessian$trace3 <- NULL
  }

  result$deriv <- the_deriv$deriv
  result
}

#' @noRd
restructure_laplace_result <- function(result_inner, control_inner = NULL) {
  keep <- setdiff(
    names(result_inner),
    c("opt", "inner_opt", "gradient", "hessian")
  )
  out <- result_inner[keep]
  inner <- result_inner$inner_opt
  if (is.null(inner) && !is.null(result_inner$opt)) {
    inner <- result_inner$opt
  }
  out$inner_opt <- inner
  if (!is.null(control_inner)) {
    out$inner_opt$control_inner <- control_inner
  }
  out$gradient <- list(
    inner = result_inner$gradient$inner,
    outer = list(grad = result_inner$gradient$outer)
  )
  out$hessian <- result_inner$hessian
  out
}

log_lik_deriv <- function(
  full_parameters,
  hessian_pack,
  grad,
  ad_pack,
  verbose = FALSE
) {
  if (!is(ad_pack, "ad_pack")) {
    stop("ad_pack must be an ad_pack object", call. = FALSE)
  }
  sz <- ad_pack@sizes
  if (length(sz) < 3L || any(is.na(sz[c("beta", "gamma", "theta")]))) {
    stop(
      "ad_pack@sizes must contain finite beta, gamma, and theta",
      call. = FALSE
    )
  }
  num_beta <- as.integer(sz["beta"])
  num_theta <- as.integer(sz["theta"])
  num_gamma <- as.integer(sz["gamma"])
  n_params <- num_beta + num_gamma + num_theta
  if (length(grad) != n_params) {
    stop(
      "length(grad) (", length(grad), ") must equal num_beta + num_gamma + ",
      "num_theta (", n_params, ")",
      call. = FALSE
    )
  }

  ldl <- as_ldl_list(hessian_pack$chol_inner)
  if (is.null(ldl$Linv)) {
    stop(
      "chol_inner$Linv required for log_lik_deriv (inner_opt with deriv = TRUE)",
      call. = FALSE
    )
  }
  half_H_inv <- hessian_pack$half_H_inv
  if (is.null(half_H_inv)) {
    half_H_inv <- half_H_inv_from_ldl(ldl)
  }
  H_inv <- hessian_pack$H_inv
  if (is.null(H_inv)) {
    H_inv <- Matrix::tcrossprod(half_H_inv)
  }
  hessian_inv <- list(
    half_H_inv = half_H_inv,
    H_inv = H_inv
  )

  seq_gamma1 <- seq.int(num_beta + 1L, length.out = num_gamma)

  the_trace <- hessian_pack$trace3
  if (is.null(the_trace)) {
    stop(
      "hessian$trace3 is missing; run inner_opt(..., deriv = TRUE) with ",
      "current adlaplace (trace is computed inside inner_opt)",
      call. = FALSE
    )
  }
  if (length(the_trace) != n_params) {
    stop(
      "length(hessian$trace3) (", length(the_trace), ") must equal num_beta + ",
      "num_gamma + num_theta (", n_params, ")",
      call. = FALSE
    )
  }

  dU <- -hessian_inv$H_inv %*% hessian_pack$outer[seq_gamma1, -seq_gamma1]

  result <- list(
    extra = c(
      list(dU = dU, trace3 = the_trace),
      hessian_inv
    )
  )

  result$deriv <- data.frame(
    d_det_upart = -as.vector(the_trace[seq_gamma1] %*% dU),
    d_det_tpart = -the_trace[-seq_gamma1]
  )
  result$deriv$grad_theta <- -grad[-seq_gamma1]
  result$deriv$grad_u <- as.vector(-grad[seq_gamma1] %*% dU)
  result$deriv$d_det <- result$deriv$d_det_upart + result$deriv$d_det_tpart
  result$deriv$d_neg_log_lik <-
    -result$deriv$grad_theta +
    result$deriv$d_det - result$deriv$grad_u
  result$deriv$d_log_lik <- -result$deriv$d_neg_log_lik

  result
}

#' Extract inner inverse Cholesky factor from Laplace output
#'
#' Returns \eqn{H^{-1/2}} for the inner (random-effect) Hessian from the output
#' of \code{\link{log_lik_laplace}}. Accepts the flat \code{half_H_inv} matrix
#' when \code{deriv = TRUE}, or reconstructs it from \code{hessian$chol_inner}.
#'
#' @param laplace Output of \code{log_lik_laplace(..., deriv = TRUE)} (or
#'   \code{inner_opt(..., deriv = TRUE)} with the same \code{hessian} layout).
#'
#' @return A sparse or dense matrix \eqn{H^{-1/2}} with \code{nrow = Ngamma}.
#'
#' @seealso \code{\link{log_lik_laplace}}, \code{\link{rmvnldl}}
#' @export
laplace_half_H_inv <- function(laplace) {
  if (!is.list(laplace)) {
    stop("laplace must be a list (output of log_lik_laplace)", call. = FALSE)
  }

  half_h <- laplace$hessian$half_H_inv
  if (is.null(half_h) && is.list(laplace$extra)) {
    half_h <- laplace$extra$half_H_inv
    if (is.null(half_h)) {
      half_h <- laplace$extra$extra$half_H_inv
    }
  }
  if (!is.null(half_h)) {
    return(half_h)
  }

  chol_inner <- NULL
  if (is.list(laplace$hessian)) {
    chol_inner <- laplace$hessian$chol_inner
  }
  if (is.null(chol_inner) && is.list(laplace$extra) && is.list(laplace$extra$hessian)) {
    chol_inner <- laplace$extra$hessian$chol_inner
  }
  if (is.null(chol_inner)) {
    stop(
      "laplace must contain hessian$half_H_inv (from log_lik_laplace(..., deriv = TRUE)) ",
      "or hessian$chol_inner",
      call. = FALSE
    )
  }

  ldl <- as_ldl_list(chol_inner)
  half_H_inv_from_ldl(ldl)
}
