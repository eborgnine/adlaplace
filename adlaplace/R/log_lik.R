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
#' @param ad_fun \code{ad_fun} object from \code{\link{ad_fun}}.
#' @param ... Additional arguments forwarded to \code{\link{log_lik_laplace}}
#'   (for example \code{data} or \code{package}).
#' @param control_inner A list of control options forwarded to the \code{control}
#'   argument of \code{\link{log_lik_laplace}} for the inner optimization.
#' @param cache An \code{\link[base]{environment}} used to warm-start and store the inner
#'   \code{gamma} solution. If \code{cache$gamma} is missing or the wrong length,
#'   it is initialized from \code{config$gamma} (if present) or zeros of length
#'   \code{ad_fun@sizes["gamma"]}. Both functions update \code{cache$gamma} after
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
#' ad_fun <- ad_fun(data, config)
#'
#' val <- outer_fn(x = x0, data = data, config = config, cache = cache, ad_fun = ad_fun)
#' gr <- outer_gr(x = x0, data = data, config = config, cache = cache, ad_fun = ad_fun)
#' }
#'
#' @name outer_optim_wrappers
#' @rdname outer_optim_wrappers
#' @export
outer_fn <- function(
  x, config, cache, ad_fun, control_inner = list(), ...
) {
  assign("last_par_fn", x, cache)
  num_gamma <- as.integer(ad_fun@sizes["gamma"])
  cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

  result <- adlaplace::inner_opt(
    parameters = x,
    gamma = cache$gamma,
    ad_fun = ad_fun,
    control = control_inner,
    deriv = FALSE,
    verbose = isTRUE(config$verbose)
  )

  assign("gamma", result$opt$solution, cache)
  result$neg_log_lik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr <- function(
  x, config, cache, ad_fun, control_inner = list(), ...
) {
  assign("last_par_gr", x, cache)
  num_gamma <- as.integer(ad_fun@sizes["gamma"])
  cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

  result <- adlaplace::log_lik_laplace(
    x = x, config = config,
    gamma = cache$gamma,
    control = control_inner,
    ad_fun = ad_fun,
    deriv = TRUE, ...
  )
  assign("gamma", result$opt$solution, cache)
  result$grad
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
#'   Length must equal \code{num_beta + num_theta} from \code{ad_fun@sizes}.
#' @param config A list of backend options (\code{verbose}, \code{package}, etc.).
#'   When \code{ad_fun} is missing, must also include \code{beta}, \code{gamma},
#'   and \code{theta} so \code{ad_fun()} can be built from \code{data}.
#' @param gamma Optional numeric vector of starting values for the inner
#'   parameter \code{gamma}. If missing, uses \code{config$gamma} when present,
#'   otherwise a zero vector of length \code{ad_fun@sizes["gamma"]}.
#' @param control List of control parameters passed *as-is* to the backend inner
#'   optimizer (e.g., \code{report.level}, \code{report.freq}). See backend
#'   documentation (e.g., \pkg{trustOptim}) for supported options.
#' @param ad_fun Optional \code{ad_fun} object from \code{\link{ad_fun}}.
#'   This is a single backend handle (no separate inner/outer handles). If
#'   missing, it will be constructed automatically using \code{data}.
#' @param data Optional data list used to build \code{ad_fun} when \code{ad_fun}
#'   is not supplied.
#' @param package Character scalar naming the backend package to use for
#'   \code{ad_fun()} and \code{inner_opt()}. Defaults to \code{"adlaplace"}.
#'   Other backends must export a compatible \code{ad_fun()} builder and \code{inner_opt()}.
#' @param deriv Logical scalar. If \code{TRUE}, include derivative quantities in
#'   the output (gradient, intermediate derivatives).
#'
#' @details
#' The default \pkg{adlaplace} backend uses a single AD handle. This function
#' passes that handle to \code{inner_opt(..., ad_fun = ad_fun)}. Parameter block
#' sizes are read from \code{ad_fun@sizes}; \code{config} is not required to
#' contain \code{beta}, \code{gamma}, or \code{theta} when \code{ad_fun} is
#' supplied.
#'
#' When \code{deriv=FALSE}, the return value is the \code{\link{inner_opt}()} list
#' (profile likelihood, nested \code{gradient} and \code{hessian}, and \code{opt}).
#' When \code{deriv=TRUE}, the same list is augmented with \code{grad},
#' \code{deriv}, and \code{extra} from the internal profile-derivative helper.
#'
#' @return A list. With \code{deriv=FALSE}, same as \code{\link{inner_opt}}. With
#' \code{deriv=TRUE}, additionally:
#' \describe{
#'   \item{grad}{Gradient of \code{neg_log_lik} w.r.t. outer \code{x}.}
#'   \item{deriv}{Data frame of profile derivative pieces.}
#'   \item{extra}{Intermediate objects (\code{dU}, \code{trace3}, Cholesky factors).}
#' }
#'
#' @note The Laplace approximation assumes the inner Hessian is positive definite
#'   at the optimum. If \code{\link{inner_opt}} fails to converge or returns a
#'   non-invertible Hessian, \code{log_lik_laplace()} issues a warning or error.
#'
#' @seealso
#' \code{\link{ad_fun}}, \code{\link{inner_opt}}
#'
#' @export
log_lik_laplace <- function(
  x,
  config = list(),
  gamma,
  control = list(report.level = 4, report.freq = 1),
  ad_fun,
  data,
  package = c(config[["package"]], "adlaplace")[1L],
  deriv = FALSE
) {
  if (missing(ad_fun)) {
    if (missing(data)) {
      stop("at least one of data and ad_fun must be supplied", call. = FALSE)
    }
    ad_fun <- adlaplace::ad_fun(data, config)
  }

  if (!is(ad_fun, "ad_fun")) {
    stop("ad_fun must be an ad_fun object", call. = FALSE)
  }
  sz <- ad_fun@sizes
  if (length(sz) < 3L || any(is.na(sz[c("beta", "gamma", "theta")]))) {
    stop(
      "ad_fun@sizes must contain finite beta, gamma, and theta",
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
    ad_fun = ad_fun,
    deriv = deriv,
    verbose = isTRUE(config[["verbose"]])
  )
  result <- result_inner[setdiff(names(result_inner), c("gradient", "hessian"))]
  result$extra <- result_inner[c("gradient", "hessian")]

  if (!deriv) {
    return(result)
  }

  the_deriv <- log_lik_deriv(
    full_parameters = result$full_parameters,
    hessian_pack = result$extra$hessian,
    grad = result$extra$gradient$outer,
    ad_fun = ad_fun,
    verbose = isTRUE(config[["verbose"]])
  )
  result <- c(result, the_deriv[setdiff(names(the_deriv), "extra")])
  result$extra <- c(result$extra, the_deriv$extra)
  result
}

log_lik_deriv <- function(
  full_parameters,
  hessian_pack,
  grad,
  ad_fun,
  verbose = FALSE
) {
  if (!is(ad_fun, "ad_fun")) {
    stop("ad_fun must be an ad_fun object", call. = FALSE)
  }
  sz <- ad_fun@sizes
  if (length(sz) < 3L || any(is.na(sz[c("beta", "gamma", "theta")]))) {
    stop(
      "ad_fun@sizes must contain finite beta, gamma, and theta",
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
  result$grad <- result$deriv$d_neg_log_lik
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
#'   \code{inner_opt(..., deriv = TRUE)} with the same \code{extra$hessian} layout).
#'
#' @return A sparse or dense matrix \eqn{H^{-1/2}} with \code{nrow = Ngamma}.
#'
#' @seealso \code{\link{log_lik_laplace}}, \code{\link{rmvnldl}}
#' @export
laplace_half_H_inv <- function(laplace) {
  if (!is.list(laplace)) {
    stop("laplace must be a list (output of log_lik_laplace)", call. = FALSE)
  }

  half_h <- laplace$extra$half_H_inv
  if (is.null(half_h) && is.list(laplace$extra)) {
    half_h <- laplace$extra$extra$half_H_inv
  }
  if (!is.null(half_h)) {
    return(half_h)
  }

  chol_inner <- NULL
  if (is.list(laplace$extra) && is.list(laplace$extra$hessian)) {
    chol_inner <- laplace$extra$hessian$chol_inner
  }
  if (is.null(chol_inner)) {
    stop(
      "laplace must contain extra$half_H_inv (from log_lik_laplace(..., deriv = TRUE)) ",
      "or extra$hessian$chol_inner",
      call. = FALSE
    )
  }

  ldl <- as_ldl_list(chol_inner)
  half_H_inv_from_ldl(ldl)
}
