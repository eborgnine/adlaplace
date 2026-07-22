#' @include log_lik.R model_data.R obs_groups.R
NULL

#' Fit a hierarchical model by Laplace-approximate maximum likelihood
#'
#' One-call interface to the \pkg{adlaplace} pipeline: parse the formula with
#' \code{\link{model_data}()}, partition observations into AD shards with
#' \code{\link{obs_groups}()}, build the AD handle with \code{\link{ad_pack}()},
#' maximize the Laplace-approximate marginal likelihood with
#' \code{\link[stats]{optim}()} using the exact profile gradient
#' (\code{\link{outer_fn}} / \code{\link{outer_gr}}), and return a fitted-model
#' object with standard accessor methods.
#'
#' @param formula Model formula built from model terms, e.g.
#'   \code{nbinom(y) ~ x + iid(g)}, a list of term objects from
#'   \code{\link{collect_terms}}, or a precomputed \code{\link{model_data}}
#'   result (in which case \code{data} is ignored).
#' @param data Data frame containing the variables referenced in \code{formula}.
#'   Not needed when \code{formula} is already a \code{model_data} result.
#' @param config List of backend options merged over the defaults
#'   \code{list(transform_theta = TRUE, verbose = verbose)}. When
#'   \code{config$obs_groups} is missing it is filled by
#'   \code{obs_groups(A, num_groups = num_groups)}, or with
#'   \code{elgm_matrix} when the model data includes one (e.g. case-crossover).
#' @param num_threads Positive integer, OpenMP thread count for the AD handle.
#' @param num_groups Target number of observation shards for
#'   \code{\link{obs_groups}()}.
#' @param control List of \code{\link[stats]{optim}} control options, merged
#'   over the defaults \code{list(maxit = 200L, parscale = <from terms>)}.
#' @param control_inner List of control options for the inner (random effects)
#'   trust-region optimizer.
#' @param method \code{\link[stats]{optim}} method; bounds are used for
#'   \code{"L-BFGS-B"} and \code{"Brent"} only.
#' @param hessian Logical. When \code{TRUE} (default), \code{optim} numerically
#'   differentiates the exact AD profile gradient at the optimum, giving the
#'   observed information used by \code{\link{vcov.adlaplace_fit}},
#'   \code{\link{summary.adlaplace_fit}}, and
#'   \code{\link{confint.adlaplace_fit}}.
#' @param verbose Logical, print progress.
#' @param na_omit Passed to \code{\link{model_data}()}.
#'
#' @return An object of class \code{"adlaplace_fit"}: a list with components
#' \describe{
#'   \item{call, formula}{The matched call and model formula.}
#'   \item{par_info}{Parameter table with \code{mle}, \code{se},
#'     \code{mle_internal}, \code{se_internal}, bounds, and \code{log}.}
#'   \item{gamma}{Named random-effect modes at the optimum.}
#'   \item{nobs}{Number of observations.}
#'   \item{details}{Laplace evaluation at the optimum, including
#'     \code{outer_opt}, \code{inner_opt}, \code{vcov}, \code{log_lik},
#'     \code{gradient}, and \code{hessian}.}
#'   \item{model_data, ad_pack, config, cache}{Objects needed by the accessor
#'     methods and for warm restarts.}
#' }
#'
#' @seealso \code{\link{summary.adlaplace_fit}}, \code{\link{predict.adlaplace_fit}},
#'   \code{\link{log_lik_laplace}}, \code{\link{model_data}}
#'
#' @examples
#' \dontrun{
#' fit <- adlaplace(
#'   nbinom(y, lower = 1e-9, init = 0.15) ~ x + iid(g, init = 0.1),
#'   data = dat, num_threads = 2L
#' )
#' summary(fit)
#' confint(fit)
#' }
#' @export
adlaplace <- function(
  formula,
  data,
  config = list(),
  num_threads = 1L,
  num_groups = 100L,
  control = list(),
  control_inner = list(maxit = 100L, report.level = 0, report.freq = 0),
  method = "L-BFGS-B",
  hessian = TRUE,
  verbose = FALSE,
  na_omit = TRUE
) {
  cl <- match.call()
  if (is.list(formula) &&
    all(c("terms", "term_data", "observations") %in% names(formula))) {
    md <- formula
    formula <- NULL
  } else {
    md <- model_data(
      formula = formula, data = data,
      verbose = verbose, na_omit = na_omit
    )
  }

  config <- utils::modifyList(
    list(transform_theta = TRUE, verbose = verbose),
    config
  )
  if (is.null(config$obs_groups)) {
    A <- md$term_data$A
    if (!is.null(A) && ncol(A) > 0L) {
      elgm <- md$term_data$elgm_matrix
      if (!is.null(elgm) && methods::is(elgm, "Matrix") && ncol(elgm) > 0L) {
        config$obs_groups <- obs_groups(
          A,
          elgm_matrix = elgm,
          num_groups = num_groups,
          min_groups = min(
            as.integer(num_groups),
            as.integer(num_threads) * 4L
          )
        )
      } else {
        config$obs_groups <- obs_groups(A, num_groups = num_groups)
      }
    }
  }
  if (is.null(config$num_threads)) {
    config$num_threads <- as.integer(num_threads)[1L]
  }

  af <- ad_pack(md, config, num_threads = config$num_threads)

  par_meta <- md$term_data$info$parameters
  labels <- par_meta$label

  control_use <- utils::modifyList(
    list(maxit = 200L, parscale = par_meta$parscale),
    control
  )

  cache <- new.env(parent = emptyenv())
  cache$gamma <- rep(0, nrow(md$term_data$info$gamma))

  optim_args <- list(
    par = par_meta$init,
    fn = outer_fn,
    gr = outer_gr,
    method = method,
    control = control_use,
    hessian = hessian,
    config = config,
    ad_pack = af,
    cache = cache,
    control_inner = control_inner
  )
  if (method %in% c("L-BFGS-B", "Brent")) {
    optim_args$lower <- par_meta$lower
    optim_args$upper <- par_meta$upper
  }
  opt <- do.call(stats::optim, optim_args)

  details <- log_lik_laplace(
    x = opt$par,
    gamma = cache$gamma,
    ad_pack = af,
    config = config,
    control = control_inner,
    deriv = TRUE
  )
  cache$gamma <- details$inner_opt$solution

  vc <- NULL
  if (hessian && !is.null(opt$hessian) && all(is.finite(opt$hessian))) {
    vc <- tryCatch(
      solve(opt$hessian),
      error = function(e) {
        warning(
          "outer Hessian is not invertible; vcov unavailable: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(vc)) {
      dimnames(vc) <- list(labels, labels)
    }
  }

  details$outer_opt <- opt
  details$vcov <- vc

  mle_internal <- stats::setNames(opt$par, labels)
  se_internal <- if (!is.null(vc)) {
    sqrt(pmax(0, diag(vc)))
  } else {
    rep(NA_real_, length(labels))
  }
  names(se_internal) <- labels

  log_idx <- par_meta$log %in% TRUE
  mle <- mle_internal
  mle[log_idx] <- exp(mle_internal[log_idx])
  se <- se_internal
  se[log_idx] <- NA_real_

  par_info <- data.frame(
    label = par_meta$label,
    mle = unname(mle),
    se = unname(se),
    mle_internal = unname(mle_internal),
    se_internal = unname(se_internal),
    init = par_meta$init,
    lower = par_meta$lower,
    upper = par_meta$upper,
    parscale = par_meta$parscale,
    log = par_meta$log,
    stringsAsFactors = FALSE
  )

  out <- list(
    call = cl,
    formula = formula,
    model_data = md,
    ad_pack = af,
    config = config,
    details = details,
    par_info = par_info,
    gamma = stats::setNames(
      details$inner_opt$solution,
      md$term_data$info$gamma$gamma_label
    ),
    nobs = length(md$term_data$y),
    cache = cache
  )
  class(out) <- "adlaplace_fit"
  out
}
