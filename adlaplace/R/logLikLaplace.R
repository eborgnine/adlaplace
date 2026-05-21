#' Log-likelihood with inner Laplace optimization
#'
#' Evaluates the (profiled) log-likelihood for a hierarchical model by solving an
#' inner optimization problem for the latent vector \code{gamma} (e.g. random
#' effects) given outer parameters \code{x} (typically \code{beta} and \code{theta}).
#' The function delegates the inner optimization to \code{inner_opt()} from the
#' selected backend package (by default \pkg{adlaplace}).
#'
#' @param x Numeric vector of outer parameters, expected to have length
#'   \code{length(config$beta) + length(config$theta)}. Elements are split
#'   into \code{beta} (first \code{length(config$beta)} elements) and \code{theta}
#'   (remaining elements) before passing to the inner optimizer.
#' @param config A list containing model dimensions/starting values and backend
#'   options. Must include \code{beta}, \code{gamma}, and \code{theta}; may also
#'   include \code{package}, \code{verbose}, and \code{num_threads}.
#' @param gamma Optional numeric vector of starting values for the inner
#'   parameter \code{gamma}. If missing, defaults to \code{config$gamma}.
#' @param control List of control parameters passed *as-is* to the backend inner
#'   optimizer (e.g., \code{report.level}, \code{report.freq}). See backend
#'   documentation (e.g., \pkg{trustOptim}) for supported options.
#' @param ad_fun Optional AD object returned by the backend \code{get_ad_fun()}.
#'   This is a single backend handle (no separate inner/outer handles). If
#'   missing, it will be constructed automatically using \code{data}.
#' @param data Optional data list used to build \code{ad_fun} when \code{ad_fun}
#'   is not supplied.
#' @param package Character scalar naming the backend package to use for
#'   \code{get_ad_fun()} and \code{inner_opt()}. Defaults to \code{"adlaplace"}.
#'   Other backends must export \code{getAdFun_r()} and \code{inner_opt()}.
#' @param deriv Logical scalar. If \code{TRUE}, include derivative quantities in
#'   the output (gradient, intermediate derivatives).
#'
#' @details
#' The parameter vector \code{x} is split into \code{beta} and \code{theta} and
#' inserted into \code{config} (as \code{config_inner$beta} and
#' \code{config_inner$theta}) before calling the backend inner optimizer.
#'
#' The default \pkg{adlaplace} backend uses a single AD handle. This function
#' passes that handle to \code{inner_opt(..., ad_fun = ad_fun)}.
#'
#' The inner objective is treated as negative joint log density; this function
#' returns both the Laplace-approximated log-likelihood (\code{logLik}) and its
#' negation (\code{fval}) for optimizer convenience.
#'
#' @return A list with components including:
#' \describe{
#'   \item{logLik}{Laplace-approximated log-likelihood at the optimized inner
#'   solution.}
#'   \item{fval}{\code{-logLik}.}
#'   \item{parameters}{The input outer parameter vector \code{x}.}
#'   \item{full_parameters}{Concatenation of \code{beta}, optimized \code{gamma},
#'   and \code{theta}.}
#'   \item{hessian}{List with \code{H} (outer Hessian sparse matrix), \code{cholInner}
#'   (sparse LDL decomposition of inner Hessian), and \code{halfLogDet} (half the
#'   log-determinant of the inner Hessian).}
#'   \item{opt}{Inner optimizer outputs (excluding Hessian objects), including
#'   \code{solution}, \code{gradient}, \code{iterations}, \code{status}, and
#'   \code{trust.radius}.}
#'   \item{grad}{When \code{deriv=TRUE}: gradient of the Laplace-approximated
#'   log-likelihood w.r.t. \code{x}, with sign convention matching \code{-logLik}.}
#'   \item{deriv, extra}{When \code{deriv=TRUE}: intermediate derivative pieces
#'   returned by \code{logLikDeriv()}.}
#' }
#'
#' @note The Laplace approximation assumes the inner Hessian \code{Hinner} is
#'   positive definite at the optimum. If \code{inner_opt()} fails to converge
#'   or returns a non-invertible Hessian, \code{logLikLaplace()} issues a warning
#'   or error.
#'
#'
#' @seealso
#' \code{\link[adlaplace]{get_ad_fun}}, \code{\link[adlaplace]{inner_opt}}
#'
#' @export
logLikLaplace <- function(x, config,
                          gamma,
                          control = list(report.level = 4, report.freq = 1),
                          ad_fun, data,
                          package = c(config$package, "adlaplace")[1],
                          deriv = FALSE) {
  Nbeta <- length(config$beta)
  Ntheta <- length(config$theta)
  Ngamma <- length(config$gamma)

  config_inner <- config
  config_inner$beta <- x[seq.int(1, length.out = Nbeta)]
  config_inner$theta <- x[seq.int(Nbeta + 1, length.out = Ntheta)]
  config_inner$deriv <- deriv

  if (!missing(gamma)) {
    config_inner$gamma <- gamma
    if (length(config$gamma) != length(config_inner$gamma)) {
      warning("gamma is the wrong length; resetting to config$gamma")
      config_inner$gamma <- config$gamma
    }
  }


  if (missing(ad_fun)) {
    if (missing(data)) {
      stop("at least one of data and ad_fun must be supplied")
    }
    ad_fun <- adlaplace::get_ad_fun(data, config, package = package)
  }

  if (deriv) {
    Sgamma1 <- seq.int(Nbeta + 1, length.out = Ngamma)
  } else {
    Sgamma1 <- seq.int(1, length.out = Ngamma)
  }


  inner_res <- inner_opt(
    x,
    config_inner$gamma,
    config = config_inner,
    control = control,
    ad_fun = ad_fun,
    deriv = deriv
  )

  Hresult <- inner_res$hessian
  if (deriv) {
    Houter <- Hresult
    Hinner <- Hresult[Sgamma1, Sgamma1]
  } else {
    Houter <- NULL
    Hinner <- Hresult
  }

  result <- list(
    logLik = inner_res$log_lik,
    parameters = x,
    hessian = list(
      outer = Houter,
      inner = Hinner,
      chol_inner = inner_res$chol_inner
    ),
    opt = inner_res[
      c(
        "iterations", "trust.radius", "status", "method",
        "fval", "half_log_det",
        "solution", "gradient", "full_parameters"
      )
    ]
  )

  if (!deriv) {
    return(result)
  }

  theDeriv <- logLikDeriv(
    full_parameters = result$opt$full_parameters,
    hessian_pack = result$hessian,
    grad = inner_res$gradient,
    config, ad_fun
  )

  result$grad <- -theDeriv$deriv$dL
  result$deriv <- theDeriv$deriv
  result$extra <- theDeriv$extra

  return(result)
}
