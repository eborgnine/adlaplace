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
#' When \code{deriv=FALSE}, the return value is the \code{inner_opt()} list
#' (profile likelihood, nested \code{gradient} and \code{hessian}, and \code{opt}).
#' When \code{deriv=TRUE}, the same list is augmented with \code{grad},
#' \code{deriv}, and \code{extra} from \code{\link[=log_lik_deriv]{log_lik_deriv()}}.
#'
#' @return A list. With \code{deriv=FALSE}, same as \code{inner_opt()}. With
#' \code{deriv=TRUE}, additionally:
#' \describe{
#'   \item{grad}{Gradient of \code{neg_log_lik} w.r.t. outer \code{x}.}
#'   \item{deriv}{Data frame of profile derivative pieces (see \code{log_lik_deriv}).}
#'   \item{extra}{Intermediate objects (\code{dU}, \code{trace3}, Cholesky factors).}
#' }
#'
#' @note The Laplace approximation assumes the inner Hessian is positive definite
#'   at the optimum. If \code{inner_opt()} fails to converge or returns a
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

  if (length(ad_fun@chol_inner_list) > 0L) {
    cil <- ad_fun@chol_inner_list
    result$extra$hessian$perm <- cil$perm
    result$extra$hessian$perm_inv <- cil$perm_inv
    if (is.list(result$extra$hessian$chol_inner)) {
      result$extra$hessian$chol_inner$perm <- cil$perm
      result$extra$hessian$chol_inner$perm_inv <- cil$perm_inv
    }
  }

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
