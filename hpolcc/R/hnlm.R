#' @import data.table
#' @importFrom adlaplace getAdFun jointLogDens grad hess
#' @importFrom adlaplace traceHinvT outer_fn outer_gr
#' @importFrom adlaplace logLikLaplace inner_opt
#' @export
hnlm <- function(
  formula,
  data,
  cc_design = c(),
  config = list(
    dirichlet = TRUE,
    boundary_is_random = FALSE,
    transform_theta = TRUE
  ),
  control = list(
    maxit = 1000,
    trace = 3,
    REPORT = 1
  ),
  control_inner = list(report.level = 0),
  for_dev = FALSE,
  ...
) {
  # Check inputs
  if (!is(formula, "formula")) {
    stop("formula must be a formula.")
  }

  config_defaults <- list(
    verbose = FALSE,
    transform_theta = TRUE,
    num_threads = 1,
    num_groups = 1000L,
    dirichlet_init = 1e-3,
    dirichlet_lower = 0,
    dirichlet_upper = Inf,
    package = "hpolcc"
  )

  config_defaults <- config_defaults[
    setdiff(names(config_defaults), names(config))
  ]
  config <- c(config, config_defaults)

  # Order the rows of data appropriately.
  if (is.character(cc_design)) {
    cc_design <- ccDesign(strat_vars = cc_design)
  }
  if (is.null(cc_design$strat_vars) &&
    is.null(cc_design$time_var)) {
    stop("Provide a valid stratification (or time) variable.")
  }

  model_terms <- adlaplace::collect_terms(
    formula,
    package = "hpolcc", verbose = config$verbose
  )

  if (!any(sapply(model_terms, methods::slot, "type") == "family")) {
    model_terms <- c(
      model_terms,
      adlaplace::overdispersion()
    )
  }

  covariates <- unique(unlist(lapply(model_terms, methods::slot, "term")))
  outcome_var <- all.vars(formula)[1]
  random_slope_terms = unique(unlist(sapply(model_terms[
    grep("^rs", sapply(model_terms, class))
  ], slot, "mult")))

  strat_time_vars <- unique(c(cc_design$strat_vars, cc_design$time_var))

  data.table::setDT(data)

  strat_time_vars <- strat_time_vars[
    order(sapply(
      data[, strat_time_vars, with = FALSE],
      function(xx) length(unique(xx))
    ), decreasing = FALSE)
  ]


  required_vars <- unique(c(covariates, strat_time_vars, random_slope_terms))
  if (config$verbose) {
    cat("variables:\n")
    print(required_vars)
  }
  # Remove rows with NA values in required variables
  data <- data[complete.cases(data[, ..required_vars])]

  if (anyNA(data[[outcome_var]])) {
    warning("missing values in outcome, treating as zeros")
    data[is.na(get(outcome_var)), (outcome_var) := 0]
  }
  data.table::setorderv(data, strat_time_vars)

  n_per_strata <- data[
    , list(
      outcome_sum = sum(get(outcome_var)),
      n_rows = .N
    ),
    by = strat_time_vars
  ]
  n_per_strata <- n_per_strata[n_per_strata$outcome_sum > 0, ]

  data_sub <- n_per_strata[data, on = strat_time_vars, nomatch = 0]

  if(!nrow(data_sub)) {
    warning("no data left after removing missings")
  }

  # setup the data for case-crossover
  if (config$verbose) {
    cat("setting strata\n")
  }

  cc_matrix <- setStrata(
    cc_design = cc_design,
    data = data_sub,
    outcome = outcome_var
  )

  if (config$verbose) {
    cat("numer per strata\n")
    print(table(diff(cc_matrix@p)))
    cat("\ncollecting terms\n")
  }

  # Use adlaplace::model_setup for standardized model preparation
  model_stuff <- adlaplace::model_setup(
    formula = formula,
    data = data_sub,
    verbose = config$verbose
  )

  # Extract components from model_setup result
  tmb_data <- model_stuff$data
  parameters_info <- model_stuff$info

  # Add case-crossover specific components
  tmb_data$elgm_matrix <- cc_matrix
  tmb_data <- formatHpolData(tmb_data)

  verbose_orig <- config$verbose
  config$verbose <- config$verbose > 1

  config$beta <- beta_setup$init
  config$theta <- theta_setup$init
  config$gamma <- rep(0, nrow(gamma_setup))

  if (is.null(beta_setup)) {
    beta_theta_names <- names(theta_setup)
    config$beta <- numeric(0)
  } else {
    beta_theta_names <- intersect(colnames(beta_setup), colnames(theta_setup))
  }

  beta_theta_names <- setdiff(beta_theta_names, "order")

  parameters_info <- list(
    beta = beta_setup,
    gamma = gamma_setup,
    theta = theta_setup,
    parameters = rbind(
      beta_setup[, beta_theta_names],
      theta_setup[, beta_theta_names]
    )
  )

  if (verbose_orig) {
    cat("getting groups...")
  }

  config$groups <- adlaplace::adFun_groups(
    ATp = tmb_data$ATp,
    elgm_matrix = tmb_data$elgm_matrix,
    Ngroups = config$num_groups
  )

  if (verbose_orig) {
    cat("done.")
  }


  cache <- new.env()
  assign("gamma", config$gamma, cache)

  config$opt <- as.list(
    parameters_info$parameters[c("init", "lower", "upper", "parscale")]
  )

  if (for_dev) {
    return(list(
      tmb_data = tmb_data,
      config = config,
      formula = formula,
      terms = model_terms,
      data = data_sub,
      control = control,
      control_inner = control_inner,
      cache = cache,
      parameters_info = parameters_info
    ))
  }

  if (verbose_orig) {
    cat(
      paste(
        "getting AD fun, ",
        paste(dim(config$groups), collapse = ","), "groups\n"
      )
    )
  }

  ad_fun <- adlaplace::getAdFun(
    tmb_data,
    config,
    package = config$package
  )

  if (!length(cache$gamma)) {
    if (verbose_orig) {
      cat("no gammas, only one layer of optimizatino")
    }
    # no gammas, no inner opt
    # Add parscale to the control if it exists in theta_info
    control_list <- list(
      fnscale = -1,
      trace = 3,
      REPORT = 1,
      maxit = 1000
    )
    control_list$parscale <- config$opt$parscale

    mle <- stats::optim(
      par = config$opt$init,
      fn = adlaplace::jointLogDens,
      gr = adlaplace::grad,
      method = "L-BFGS-B",
      lower = config$opt$lower,
      upper = config$opt$upper,
      backendContext = ad_fun,
      control = control_list
    )
    mle$hessian <- adlaplace::hess(
      mle$par,
      backendContext = ad_fun
    )
    return(mle)
  }

  cache <- new.env(parent = emptyenv())
  cache$gamma <- config$gamma
  if (verbose_orig) {
    cat("optimizing")
    cat(" initial, lower , upper\n")
    to_print <- do.call(cbind, config$opt)
    rownames(to_print) <- parameters_info$parameters$label
    print(to_print)
  }

  control$parscale <- config$theta_info$parscale

  mle <- try(stats::optim(
    par = config$opt$init,
    fn = adlaplace::outer_fn,
    gr = adlaplace::outer_gr,
    method = "L-BFGS-B",
    control = control,
    lower = config$opt$lower,
    upper = config$opt$upper,
    config = config,
    adFun = ad_fun,
    cache = cache,
    control_inner = control_inner
  ))

  config$gamma <- get("gamma", cache)

  result <- list(
    opt = mle,
    objects = list(
      tmb_data = tmb_data,
      config = config,
      formula = formula,
      terms = model_terms,
      parameters_info = parameters_info,
      random_info = random_info,
      control_inner = control$inner,
      control = control,
      cache = cache,
      data=data_sub
    )
  )
  if (verbose_orig) {
    cat("done")
  }

  result$extra <- try(adlaplace::logLikLaplace(
    result$opt[[grep("solution|par", names(result$opt), value = TRUE)[1]]],
    gamma = result$objects$config$gamma,
    data = result$objects$tmb_data,
    config = result$objects$config,
    control = control_inner,
    adFun = ad_fun,
    deriv = 1
  ))
  result$extra$parameters_orig <- result$parameters
  result$parameters <- try(format_parameters(result))

    result$hessian_parameters <- try(
      numDeriv::jacobian(
        adlaplace::outer_gr,
        x = mle$solution,
        package = "hpolcc",
        data = tmb_data,
        config = config,
        control_inner = control_inner,
        adFun = ad_fun,
        cache = cache
      )
    )

  result$sample <- try(cond_sim_iwp(
    fit = result,
    n = c(result$objects$config$num_sim, 500)[1]
  ))

  result
}
