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
#'   \code{list(transform_theta = TRUE, verbose = verbose, num_threads = 1L,
#'   num_shards = 100L)}. When \code{config$obs_groups} is missing it is filled
#'   by \code{obs_groups(A, num_shards = config$num_shards)}, or with
#'   \code{elgm_matrix} when the model data includes one (e.g. case-crossover).
#'   \code{config$num_threads} is the OpenMP thread count for the AD handle;
#'   \code{config$num_shards} is the target number of observation shards for
#'   \code{\link{obs_groups}()}.
#' @param control List of \code{\link[stats]{optim}} control options, merged
#'   over the defaults \code{list(maxit = 200L, parscale = <from terms>)}.
#' @param control_inner List of control options for the inner (random effects)
#'   trust-region optimizer.
#' @param method \code{\link[stats]{optim}} method; bounds are used for
#'   \code{"L-BFGS-B"} and \code{"Brent"} only.
#' @param hessian Logical. When \code{NULL} (default), set to \code{TRUE} if
#'   there are at most five outer parameters and \code{FALSE} otherwise.
#'   When \code{TRUE}, \code{\link[stats]{optimHess}} numerically differentiates
#'   the exact AD profile gradient at the optimum, giving the observed
#'   information used by \code{\link{vcov.adlaplace_fit}},
#'   \code{\link{summary.adlaplace_fit}}, and \code{\link{confint.adlaplace_fit}}.
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
#' If \code{optim} throws (for example L-BFGS-B when \code{fn} is not finite),
#' the same class is returned with \code{par_info} estimates set to \code{NA},
#' \code{details$outer_opt = NULL}, and \code{details$error} set to the
#' condition message. \code{ad_pack} and \code{cache} are still present so the
#' Laplace objective can be inspected or re-run without rebuilding tapes:
#' \preformatted{
#' stats::optim(
#'   fit$cache$last_par_fn,
#'   fn = adlaplace::outer_fn, gr = adlaplace::outer_gr,
#'   method = fit$method, control = fit$control,
#'   lower = fit$par_info$lower, upper = fit$par_info$upper,
#'   config = fit$config, ad_pack = fit$ad_pack,
#'   cache = fit$cache, control_inner = fit$control_inner
#' )
#' }
#'
#' @seealso \code{\link{summary.adlaplace_fit}}, \code{\link{predict.adlaplace_fit}},
#'   \code{\link{log_lik_laplace}}, \code{\link{model_data}}
#'
#' @examples
#' \dontrun{
#' fit <- adlaplace(
#'   nbinom(y, lower = 1e-9, init = 0.15) ~ x + iid(g, init = 0.1),
#'   data = dat, config = list(num_threads = 2L)
#' )
#' summary(fit)
#' confint(fit)
#' }
#' @export
adlaplace <- function(
  formula,
  data,
  config = list(),
  control = list(),
  control_inner = list(maxit = 100L, report.level = 0, report.freq = 0),
  method = "L-BFGS-B",
  hessian = NULL,
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
    list(
      transform_theta = TRUE,
      verbose = verbose,
      num_threads = 1L,
      num_shards = 100L
    ),
    config
  )
  config$num_threads <- as.integer(config$num_threads)[1L]
  config$num_shards <- as.integer(config$num_shards)[1L]
  config <- ensure_config_obs_groups(
    config,
    A = md$term_data$A,
    elgm_matrix = md$term_data$elgm_matrix,
    num_shards = config$num_shards,
    num_threads = config$num_threads
  )

  af <- ad_pack(md, config, num_threads = config$num_threads)

  par_meta <- md$term_data$info$parameters
  labels <- par_meta$label
  if (is.null(hessian)) {
    hessian <- length(par_meta$init) <= 5L
  }
  hessian <- isTRUE(hessian)

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
    hessian = FALSE,
    config = config,
    ad_pack = af,
    cache = cache,
    control_inner = control_inner
  )
  if (method %in% c("L-BFGS-B", "Brent")) {
    optim_args$lower <- par_meta$lower
    optim_args$upper <- par_meta$upper
  }
  opt <- tryCatch(
    do.call(stats::optim, optim_args),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    return(.adlaplace_fit_optim_failed(
      err = opt,
      cl = cl,
      formula = formula,
      md = md,
      af = af,
      config = config,
      control_use = control_use,
      control_inner = control_inner,
      method = method,
      par_meta = par_meta,
      cache = cache
    ))
  }

  # Always run a full Laplace evaluation at the MLE for fit$details
  # (hessians, H_inv, dU, etc.). Slim outer_fn/outer_gr cache must not be reused.
  details <- log_lik_laplace(
    x = opt$par,
    gamma = cache$gamma,
    ad_pack = af,
    config = config,
    control = control_inner,
    deriv = TRUE,
    return_hessians = TRUE
  )
  cache$gamma <- details$inner_opt$solution
  cache$fg_x <- as.numeric(opt$par)
  cache$neg_log_lik <- details$neg_log_lik
  cache$d_neg_log_lik <- as.numeric(details$deriv$d_neg_log_lik)
  if (exists("fg_result", envir = cache, inherits = FALSE)) {
    rm("fg_result", envir = cache)
  }

  if (hessian) {
    details_at_opt <- details
    gamma_at_opt <- details$inner_opt$solution
    hess <- tryCatch(
      stats::optimHess(
        par = opt$par,
        fn = outer_fn,
        gr = outer_gr,
        control = control_use,
        config = config,
        ad_pack = af,
        cache = cache,
        control_inner = control_inner
      ),
      error = function(e) e
    )
    if (!inherits(hess, "error") && all(is.finite(hess))) {
      opt$hessian <- hess
    } else if (inherits(hess, "error")) {
      warning(
        "outer Hessian failed; vcov unavailable: ",
        conditionMessage(hess),
        call. = FALSE
      )
    }
    # Restore the optimum evaluation after finite-difference probes.
    details <- details_at_opt
    cache$gamma <- gamma_at_opt
    cache$fg_x <- as.numeric(opt$par)
    cache$neg_log_lik <- details_at_opt$neg_log_lik
    cache$d_neg_log_lik <- as.numeric(details_at_opt$deriv$d_neg_log_lik)
  }

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
    control = control,
    control_inner = control_inner,
    method = method,
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

#' Build an adlaplace_fit after outer optim throws
#'
#' @noRd
.adlaplace_fit_optim_failed <- function(
    err,
    cl,
    formula,
    md,
    af,
    config,
    control_use,
    control_inner,
    method,
    par_meta,
    cache) {
  err_msg <- conditionMessage(err)
  labels <- par_meta$label
  last_x <- if (!is.null(cache$last_par_fn)) {
    as.numeric(cache$last_par_fn)
  } else {
    as.numeric(par_meta$init)
  }
  names(last_x) <- labels
  log_idx <- par_meta$log %in% TRUE
  last_nat <- last_x
  last_nat[log_idx] <- exp(last_x[log_idx])
  g <- cache$gamma

  message("outer optim failed: ", err_msg)
  message(
    "last outer parameters (internal):\n",
    paste(sprintf("  %s = %s", labels, format(last_x, digits = 6L)), collapse = "\n")
  )
  message(
    "last outer parameters (natural scale for log slots):\n",
    paste(sprintf("  %s = %s", labels, format(last_nat, digits = 6L)), collapse = "\n")
  )
  if (!is.null(cache$last_par_gr)) {
    gr <- as.numeric(cache$last_par_gr)
    names(gr) <- labels
    message(
      "last outer gradient:\n",
      paste(sprintf("  %s = %s", labels, format(gr, digits = 6L)), collapse = "\n")
    )
  }
  message(
    "gamma length: ", length(g),
    "; any non-finite gamma: ", any(!is.finite(g)),
    "; any non-finite last x: ", any(!is.finite(last_x))
  )

  n_par <- length(labels)
  par_info <- data.frame(
    label = par_meta$label,
    mle = rep(NA_real_, n_par),
    se = rep(NA_real_, n_par),
    mle_internal = rep(NA_real_, n_par),
    se_internal = rep(NA_real_, n_par),
    init = par_meta$init,
    lower = par_meta$lower,
    upper = par_meta$upper,
    parscale = par_meta$parscale,
    log = par_meta$log,
    stringsAsFactors = FALSE
  )
  gamma_lab <- md$term_data$info$gamma$gamma_label
  if (is.null(g)) {
    g <- numeric()
    gamma_lab <- character()
  }
  out <- list(
    call = cl,
    formula = formula,
    model_data = md,
    ad_pack = af,
    config = config,
    control = control_use,
    control_inner = control_inner,
    method = method,
    details = list(outer_opt = NULL, vcov = NULL, error = err_msg),
    par_info = par_info,
    gamma = stats::setNames(as.numeric(g), gamma_lab),
    nobs = length(md$term_data$y),
    cache = cache
  )
  class(out) <- "adlaplace_fit"
  out
}
