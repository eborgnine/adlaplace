#' Fit Hierarchical Non-Linear Models
#'
#' @description
#' Fits a hierarchical case-crossover model with a \code{dirichlet_multinom} response
#' term, using \pkg{adlaplace} automatic differentiation and Laplace approximation.
#'
#' @param formula Model formula with a \code{dirichlet_multinom(...)} response on the LHS.
#' @param data Data frame containing variables referenced in \code{formula}.
#' @param config Configuration list passed to \pkg{adlaplace}. Common entries:
#'   \describe{
#'     \item{\code{transform_theta}}{If \code{TRUE} (default), optimize
#'       log-scale hyperparameters.}
#'     \item{\code{num_threads}}{OpenMP threads for inner optimization and
#'       derivative evaluation (default \code{1L}).}
#'     \item{\code{num_shards}}{Target number of observation shards for
#'       parallel evaluation (default \code{1000L}).}
#'     \item{\code{num_sim}}{Number of conditional simulation draws from
#'       \code{adlaplaceHgp::cond_sim_iwp} (default \code{500}).}
#'     \item{\code{verbose}}{Verbosity level; values above \code{1} enable
#'       extra \pkg{adlaplace} logging.}
#'   }
#' @param control Control list passed to outer \code{optim}.
#' @param control_inner Control list passed to inner optimization.
#' @param for_dev If \code{TRUE}, return intermediate objects for development.
#' @param ... Unused.
#'
#' @return
#' When \code{for_dev = TRUE}, a development bundle of class \code{c("hnlm_dev", "hnlm")}.
#' Otherwise a fitted object of class \code{hnlm} with primary slots
#' \code{coefficients}, \code{log_lik}, \code{optim}, \code{converged},
#' \code{extra} (Laplace output), \code{hessian}, \code{info}, optional
#' \code{sample}, and \code{call} (fit metadata including \code{config},
#' \code{terms}, and the original \code{match.call}).
#' @export
hnlm <- function(
  formula,
  data,
  config = list(transform_theta = TRUE),
  control = list(
    maxit = 1000,
    trace = 3,
    REPORT = 1
  ),
  control_inner = list(report.level = 0),
  for_dev = FALSE,
  ...
) {
  call <- match.call()
  config <- hnlm_merge_config(config)

  model_terms <- if (methods::is(formula, "formula")) {
    adlaplace::collect_terms(
      formula = stats::update.formula(formula, . ~ . - 1),
      verbose = config$verbose
    )
  } else {
    formula
  }

  resp_idx <- which(vapply(
    model_terms,
    function(xx) methods::is(xx, "dirichlet_multinom"),
    logical(1)
  ))
  if (length(resp_idx) != 1L) {
    stop("formula must include exactly one dirichlet_multinom(...) response term",
      call. = FALSE
    )
  }
  resp <- model_terms[[resp_idx]]
  if (!length(resp@by)) {
    stop("dirichlet_multinom(..., by = ...) is required", call. = FALSE)
  }

  model_data <- adlaplace::model_data(
    formula = model_terms,
    data = data,
    verbose = config$verbose,
    na_omit = TRUE
  )

  if (config$verbose) {
    cat("number per strata\n")
    print(table(diff(model_data$data$elgm_matrix@p)))
    cat("\ncollecting terms\n")
  }

  verbose_orig <- config$verbose
  config$verbose <- config$verbose > 1L

  config$opt <- as.list(
    model_data$data$info$parameters[c("init", "lower", "upper", "parscale")]
  )
  control$parscale <- config$opt$parscale

  cache <- new.env(parent = emptyenv())
  cache$gamma <- rep(0, nrow(model_data$data$info$gamma))

  if (verbose_orig) {
    cat("getting shards...")
  }
  config$shards <- hnlm_build_shards(model_data, config)

  if (for_dev) {
    return(hnlm_dev_bundle(
      model_data = model_data,
      config = config,
      formula = formula,
      control = control,
      control_inner = control_inner,
      cache = cache,
      call = call,
      verbose_orig = verbose_orig
    ))
  }

  if (!length(cache$gamma)) {
    if (verbose_orig) {
      cat("no gammas, only one layer of optimization\n")
    }
    return(hnlm_fit_flat(
      model_data = model_data,
      config = config,
      control = control,
      control_inner = control_inner,
      cache = cache,
      call = call
    ))
  }

  hnlm_fit_laplace(
    model_data = model_data,
    config = config,
    control = control,
    control_inner = control_inner,
    cache = cache,
    call = call,
    verbose_orig = verbose_orig
  )
}

#' @keywords internal
hnlm_merge_config <- function(config) {
  defaults <- list(
    verbose = FALSE,
    transform_theta = TRUE,
    num_threads = 1L,
    num_shards = 1000L,
    num_sim = 500L
  )
  defaults <- defaults[setdiff(names(defaults), names(config))]
  c(config, defaults)
}

#' @keywords internal
hnlm_build_shards <- function(model_data, config) {
  adlaplace::ad_shards(
    A = model_data$data$A,
    elgm_matrix = model_data$data$elgm_matrix,
    num_shards = config$num_shards,
    min_groups = min(config$num_shards, config$num_threads * 4L)
  )
}

#' @keywords internal
hnlm_dev_bundle <- function(
  model_data,
  config,
  formula,
  control,
  control_inner,
  cache,
  call,
  verbose_orig
) {
  if (verbose_orig) {
    cat("done.\n")
    cat(
      "getting AD fun, ",
      paste(dim(config$shards), collapse = ","), "shards..."
    )
  }
  ad_fun <- adlaplace::ad_fun(
    x = model_data,
    config,
    num_threads = config$num_threads
  )
  if (verbose_orig) {
    cat("done.\n")
  }
  structure(
    list(
      model_data = model_data,
      config = config,
      formula = formula,
      data = model_data$data$data,
      terms = model_data$terms,
      control = control,
      control_inner = control_inner,
      cache = cache,
      ad_fun = ad_fun,
      call = call
    ),
    class = c("hnlm_dev", "hnlm", "list")
  )
}

#' @keywords internal
hnlm_fit_flat <- function(
  model_data,
  config,
  control,
  control_inner,
  cache,
  call
) {
  ad_fun <- adlaplace::ad_fun(
    x = model_data,
    config,
    num_threads = config$num_threads
  )
  optim_result <- stats::optim(
    par = config$opt$init,
    fn = function(x) adlaplace::joint_log_dens(ad_fun, x),
    gr = function(x) adlaplace::grad(ad_fun, x),
    method = "L-BFGS-B",
    lower = config$opt$lower,
    upper = config$opt$upper,
    control = control
  )
  hessian_outer <- adlaplace::hessian(ad_fun, optim_result$par)
  log_lik <- -optim_result$value

  coefficients <- try(adlaplace::format_parameters(
    info = ad_fun@info,
    gamma = cache$gamma,
    parameters = optim_result$par
  ))
  coefficients <- hnlm_attach_parameter_se(coefficients, hessian_outer)

  hnlm_assemble(
    coefficients = coefficients,
    log_lik = log_lik,
    optim_result = optim_result,
    laplace = list(log_lik = log_lik),
    hessian = list(outer = hessian_outer, inner = NULL, var_iid = NULL),
    ad_fun = ad_fun,
    sample = NULL,
    call = call,
    config = config,
    model_data = model_data,
    control = control,
    control_inner = control_inner,
    cache = cache
  )
}

#' @keywords internal
hnlm_fit_laplace <- function(
  model_data,
  config,
  control,
  control_inner,
  cache,
  call,
  verbose_orig
) {
  if (verbose_orig) {
    cat("optimizing initial, lower, upper\n")
    to_print <- do.call(cbind, config$opt)
    rownames(to_print) <- model_data$data$info$parameters$label
    print(to_print)
    cat("threads: ", config$num_threads, "\n")
  }

  fit <- adlaplace::adlaplace(
    formula = model_data,
    config = config,
    num_threads = config$num_threads,
    num_shards = config$num_shards,
    control = control,
    control_inner = control_inner,
    method = "L-BFGS-B",
    hessian = TRUE,
    verbose = verbose_orig
  )

  ad_fun <- fit$ad_fun
  optim_result <- fit$optim
  cache$gamma <- fit$cache$gamma
  config$gamma <- fit$cache$gamma
  laplace <- fit$details

  coefficients <- try(adlaplace::format_parameters(
    info = ad_fun@info,
    gamma = fit$cache$gamma,
    parameters = optim_result$par
  ))

  if (verbose_orig) {
    cat("hessian of parameters\n")
    cat("conditional samples\n")
  }

  sample <- hnlm_build_simulation(laplace, model_data, config)
  var_out <- hnlm_build_variance(
    coefficients = coefficients,
    laplace = laplace,
    model_data = model_data,
    hessian_outer = optim_result$hessian,
    vcov = fit$vcov
  )
  optim_result$hessian <- NULL

  hnlm_assemble(
    coefficients = var_out$coefficients,
    log_lik = laplace$log_lik,
    optim_result = optim_result,
    laplace = laplace,
    hessian = var_out$hessian,
    ad_fun = ad_fun,
    sample = sample,
    call = call,
    config = config,
    model_data = model_data,
    control = control,
    control_inner = control_inner,
    cache = cache
  )
}

#' @keywords internal
hnlm_build_simulation <- function(laplace, model_data, config) {
  try(adlaplaceHgp::cond_sim_iwp(
    fit = laplace,
    model_data = model_data,
    n = config$num_sim
  ))
}

#' @keywords internal
hnlm_build_variance <- function(
  coefficients,
  laplace,
  model_data,
  hessian_outer,
  vcov = NULL
) {
  hessian <- list(outer = hessian_outer, inner = NULL, var_iid = NULL)
  coefficients <- hnlm_attach_parameter_se(
    coefficients, hessian$outer,
    vcov = vcov
  )

  n_beta <- nrow(model_data$data$info$beta)
  if (!inherits(coefficients, "try-error") &&
    !is.null(coefficients$gamma) &&
    nrow(coefficients$gamma) > 0L) {
    h_pack <- laplace$extra$hessian
    if (is.list(h_pack) && !is.null(h_pack$inner)) {
      hessian$inner <- h_pack$inner
    } else if (is.list(h_pack) && !is.null(h_pack$H)) {
      n_gamma <- nrow(coefficients$gamma)
      if (nrow(h_pack$H) == n_gamma) {
        hessian$inner <- h_pack$H
      } else {
        seq_inner <- seq(
          from = max(c(0L, n_beta)) + 1L,
          length.out = n_gamma
        )
        hessian$inner <- h_pack$H[seq_inner, seq_inner, drop = FALSE]
      }
    }
    which_is_iid <- grepl("iid", coefficients$gamma$model)
    if (any(which_is_iid) &&
      !is.null(hessian$inner) &&
      requireNamespace("WoodburyMatrix", quietly = TRUE)) {
      H_inner <- hessian$inner
      Dinv <- H_inner[which_is_iid, which_is_iid]
      for_var_years <- WoodburyMatrix::WoodburyMatrix(
        A = Matrix::solve(Dinv),
        B = H_inner[!which_is_iid, !which_is_iid],
        X = H_inner[which_is_iid, !which_is_iid],
        symmetric = TRUE
      )
      hessian$var_iid <- WoodburyMatrix::solve(for_var_years)
    }
  }

  list(hessian = hessian, coefficients = coefficients)
}

#' @keywords internal
hnlm_assemble <- function(
  coefficients,
  log_lik,
  optim_result,
  laplace,
  hessian,
  ad_fun,
  sample,
  call,
  config,
  model_data,
  control,
  control_inner,
  cache
) {
  info <- if (methods::is(ad_fun, "ad_fun")) {
    ad_fun@info
  } else if (!is.null(model_data$data$info)) {
    model_data$data$info
  } else {
    NULL
  }
  structure(
    list(
      coefficients = coefficients,
      log_lik = log_lik,
      optim = optim_result,
      converged = isTRUE(optim_result$convergence == 0L),
      extra = laplace,
      hessian = hessian,
      info = info,
      sample = sample,
      call = list(
        config = config,
        terms = model_data$terms,
        control = control,
        control_inner = control_inner,
        cache = cache,
        call = call
      )
    ),
    class = c("hnlm", "list")
  )
}
