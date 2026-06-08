#' Fit Hierarchical Non-Linear Models
#'
#' @description
#' Fits a hierarchical case-crossover model with a \code{dirichlet_multinom} response
#' term, using \pkg{adlaplace} automatic differentiation and Laplace approximation.
#'
#' @param formula Model formula with a \code{dirichlet_multinom(...)} response on the LHS.
#' @param data Data frame containing variables referenced in \code{formula}.
#' @param config Configuration list. Common entries: \code{transform_theta},
#'   \code{num_threads}, \code{num_shards}, \code{verbose}.
#' @param control Control list passed to outer \code{optim}.
#' @param control_inner Control list passed to inner optimization.
#' @param for_dev If \code{TRUE}, return intermediate objects for development.
#' @param ... Unused.
#'
#' @return Fitted model list, or development objects when \code{for_dev = TRUE}.
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
  config_defaults <- list(
    verbose = FALSE,
    transform_theta = TRUE,
    num_threads = 1L,
    num_shards = 1000L
  )
  config_defaults <- config_defaults[
    setdiff(names(config_defaults), names(config))
  ]
  config <- c(config, config_defaults)

  if (methods::is(formula, "formula")) {
    model_terms <- adlaplace::collect_terms(
      formula = stats::update.formula(formula, . ~ . - 1),
      package = c("hpolcc", "adlaplaceHgp"),
      verbose = config$verbose
    )
  } else {
    model_terms <- formula
  }

  resp_idx <- which(vapply(
    model_terms,
    function(xx) methods::is(xx, "dirichlet_multinom"),
    logical(1)
  ))
  if (length(resp_idx) != 1L) {
    stop("formula must include exactly one dirichlet_multinom(...) response term",
      call. = FALSE)
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
  cc_matrix <- model_data$data$elgm_matrix

  if (config$verbose) {
    cat("number per strata\n")
    print(table(diff(cc_matrix@p)))
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

  config$shards <- adlaplace::ad_shards(
    A = model_data$data$A,
    elgm_matrix = cc_matrix,
    num_shards = config$num_shards,
    min_groups = min(config$num_shards, config$num_threads * 4L)
  )

  if (verbose_orig) {
    cat("done.\n")
    cat(
      "getting AD fun, ",
      paste(dim(config$shards), collapse = ","), "shards\n"
    )
  }

  ad_fun <- adlaplace::ad_fun(
    model_data,
    config,
    num_threads = config$num_threads
  )

  if (for_dev) {
    return(list(
      model_data = model_data,
      config = config,
      formula = formula,
      data = model_data$data$data,
      control = control,
      control_inner = control_inner,
      cache = cache,
      ad_fun = ad_fun
    ))
  }

  if (!length(cache$gamma)) {
    if (verbose_orig) {
      cat("no gammas, only one layer of optimization\n")
    }

    mle <- stats::optim(
      par = config$opt$init,
      fn = function(x) adlaplace::joint_log_dens(ad_fun, x),
      gr = function(x) adlaplace::grad(ad_fun, x),
      method = "L-BFGS-B",
      lower = config$opt$lower,
      upper = config$opt$upper,
      control = control
    )
    mle$hessian <- adlaplace::hessian(ad_fun, mle$par)
    return(mle)
  }

  if (verbose_orig) {
    cat("optimizing initial, lower, upper\n")
    to_print <- do.call(cbind, config$opt)
    rownames(to_print) <- model_data$data$info$parameters$label
    print(to_print)
    cat("threads: ", config$num_threads, "\n")
  }

  mle <- try(stats::optim(
    par = config$opt$init,
    fn = adlaplace::outer_fn,
    gr = adlaplace::outer_gr,
    method = "L-BFGS-B",
    control = control,
    lower = config$opt$lower,
    upper = config$opt$upper,
    config = config,
    ad_fun = ad_fun,
    cache = cache,
    control_inner = control_inner
  ))

  config$gamma <- get("gamma", cache)

  result <- list(
    opt = mle,
    objects = list(
      config = config,
      formula = formula,
      terms = model_data$terms,
      parameters_info = model_data$data$info,
      random_info = model_data$data$info$gamma,
      control_inner = control_inner,
      control = control,
      cache = cache,
      ad_fun = ad_fun,
      model_data = model_data
    )
  )
  if (verbose_orig) {
    cat("one last evaluation of likelihood\n")
  }

  par_name <- grep("solution|par", names(result$opt), value = TRUE)[1]
  result$extra <- try(adlaplace::log_lik_laplace(
    x = result$opt[[par_name]],
    gamma = result$objects$cache$gamma,
    config = result$objects$config,
    control = result$objects$control_inner,
    ad_fun = ad_fun,
    deriv = TRUE
  ))
  result$extra$parameters_orig <- result$parameters
  result$parameters <- try(format_parameters(result))

  if (verbose_orig) {
    cat("hessian of parameters\n")
  }
  result$hessian_parameters <- try(
    Matrix::forceSymmetric(numDeriv::jacobian(
      func = adlaplace::outer_gr,
      x = result$opt[[par_name]],
      config = result$objects$config,
      control_inner = result$objects$control_inner,
      ad_fun = ad_fun,
      cache = result$objects$cache
    ), "U")
  )
  if (verbose_orig) {
    cat("conditional samples\n")
  }
  result$sample <- try(adlaplaceHgp::cond_sim_iwp(
    fit = result,
    n = c(result$objects$config$num_sim, 500)[1]
  ))

  seq_inner <- seq(
    from = max(c(0, nrow(result$parameters$beta))) + 1,
    length.out = nrow(result$parameters$gamma)
  )

  H_inner <- result$extra$hessian$H_inner <-
    result$extra$hessian$H[seq_inner, seq_inner]

  which_is_iid <- grepl("iid", result$parameters$gamma$model)
  if (any(which_is_iid) && requireNamespace("WoodburyMatrix", quietly = TRUE)) {
    Dinv <- H_inner[which_is_iid, which_is_iid]

    for_var_years <- WoodburyMatrix::WoodburyMatrix(
      A = Matrix::solve(Dinv),
      B = H_inner[!which_is_iid, !which_is_iid],
      X = H_inner[which_is_iid, !which_is_iid],
      symmetric = TRUE
    )
    result$extra$hessian$var_iid <- WoodburyMatrix::solve(for_var_years)
  }

  result
}
