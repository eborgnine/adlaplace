
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
    fnscale = -1,
    trace = 3,
    REPORT = 1,
    maxit = 1000,
    parscale = NULL
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
    Ngroups = 1e4,
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

  model_terms <- collect_terms(formula)

  if(!any(sapply(model_terms, slot, "type") == "familly")) {
    model_terms = c(model_terms,
      overdispersion()
    )
  }

  covariates <- unique(unlist(lapply(model_terms, slot, "term")))
  outcome_var <- all.vars(formula)[1]

  strat_time_vars <- unique(c(cc_design$strat_vars, cc_design$time_var))

  data.table::setDT(data)

  strat_time_vars <- strat_time_vars[
    order(sapply(
      data[, strat_time_vars, with = FALSE],
      function(xx) length(unique(xx))
    ), decreasing = FALSE)
  ]

  required_vars <- unique(c(covariates, strat_time_vars))
  data <- data[complete.cases(data[, ..required_vars])]
  if (anyNA(data[[outcome_var]])) {
    warning("missing values in outcome, treating as zeros")
    data[is.na(get(outcome_var)), (outcome_var) := 0]
  }
  data.table::setorderv(data, strat_time_vars)

  n_per_strata <- data[
    , .(
      outcome_sum = sum(get(outcome_var)),
      n_rows = .N
    ),
    by = strat_time_vars
  ]
  n_per_strata <- n_per_strata[n_per_strata$outcome_sum > 0, ]

  data_sub <- n_per_strata[data, on = strat_time_vars, nomatch = 0]


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

  # add group info to model
  # to do, also for iid
  # create a function terms_extra(term_list, data)
  var_seq <- sapply(model_terms, slot, "by")
  has_group_var <- which(unlist(lapply(var_seq, length))>0)
  var_seq = unlist(var_seq[has_group_var])

  unique_values <- lapply(unique(var_seq), function(xx) unique(data_sub[[xx]]))
  names(unique_values) <- unique(var_seq)

  unique_values_string <- lapply(unique_values, function(xx) {
    formatC(xx, width = ceiling(log10(max(xx))), flag = "0")
  })

  for (D in has_group_var) {
    model_terms[[D]]@by_levels <- unique_values[[model_terms[[D]]@by]]
    model_terms[[D]]@by_labels <- unique_values_string[[model_terms[[D]]@by]]
  }

  precision_list <- lapply(model_terms, precision)

  theta_info_list <- lapply(model_terms, theta_info)

  # to do: class for overdisp.  can add f(model="overdisp", options), if no f(model="overdisp") in formula it's added by default
  # to do: add fixed=value to classes. (so can do overdisp with fixed=0)
  theta_info_list <- c(
    theta_info_list,
    list(overdisp = data.frame(
      var = NA,
      model = "overdisp",
      name = "overdisp",
      global = TRUE,
      order = NA,
      init = config$dirichlet_init,
      lower = config$dirichlet_lower,
      upper = config$dirichlet_upper,
      parscale = 1
    ))
  )

  theta_setup <- do.call(rbind, theta_info_list)
  theta_setup$id <- seq.int(0, len = nrow(theta_setup))

  beta_setup <- do.call(rbind, lapply(model_terms, beta_info, data = data_sub))

  random_info_list <- lapply(model_terms, random_info, data = data_sub)
  gamma_setup <- do.call(rbind, random_info_list)
  gamma_setup$id <- seq.int(0, len = nrow(gamma_setup))
  gamma_setup$theta_id <- theta_setup[match(
    gamma_setup$name, theta_setup$name
  ), "id"]

  terms_with_gamma <- sapply(model_terms, slot, "name") %in%
    unique(gamma_setup$name)

  design_list_x <- lapply(model_terms[!terms_with_gamma], design,
    data = data_sub
  )
  design_list_a <- parallel::mclapply(model_terms[terms_with_gamma], design,
    data = data_sub,
    mc.cores = config$num_threads
  )

  gamma_dims <- data.frame(
    design=sapply(design_list_a, ncol),
    gamma=unlist(sapply(random_info_list, nrow))
  )
  if(!all(gamma_dims$design == gamma_dims$gamma)) {
    warning("gamma sizes wrong")
    print(gamma_dims)
  }
  a_matrix <- do.call(cbind, design_list_a)
  x_matrix <- do.call(cbind, design_list_x)
  
  beta_reorder = match(colnames(x_matrix), beta_setup$beta_name)
  beta_setup = beta_setup[beta_reorder, ]
  gamma_reorder = match(colnames(a_matrix), gamma_setup$gamma_label)
  gamma_setup = gamma_setup[gamma_reorder, ]

  gamma_setup_sub = gamma_setup[!is.na(gamma_setup$theta_id), ]
  gamma_theta_map = Matrix::sparseMatrix(
    i = gamma_setup_sub$id, j=gamma_setup_sub$theta_id,
    x = 1.0,
    index1=FALSE
  )

  # Create block-diagonal precision matrix, excluding NULL elements
  valid_precision <- precision_list[!sapply(precision_list, is.null)]
  if (length(valid_precision) > 0) {
    precision_matrix <- do.call(Matrix::bdiag, valid_precision)
  } else {
    precision_matrix <- Matrix::Diagonal(nrow = 0, ncol = 0)
  }

  if (any(duplicated(colnames(a_matrix)))) {
    warning("duplicated column names of random effects design matrix, perhaps same term in the model twice?")
  }
  if(!(all(colnames(a_matrix) == gamma_setup$gamma_name)) ) {
    warning("some names of A don't match up")
    print(table(colnames(a_matrix) %in% gamma_setup$gamma_label))
    print(str(setdiff(colnames(a_matrix), gamma_setup$gamma_label)))
    print(str(setdiff(gamma_setup$gamma_label, colnames(a_matrix))))
  }

  if (any(duplicated(colnames(x_matrix)))) {
    warning("duplicated column names of fixed effects design matrix, perhaps same term in the model twice?")
  }
  if (!(all(colnames(x_matrix) == beta_setup$beta_name))) {
    warning("some names of X don't match beta_setup")
    print(table(colnames(x_matrix) %in% beta_setup$beta_name))
    print(str(setdiff(colnames(x_matrix), beta_setup$beta_name)))
    print(str(setdiff(beta_setup$beta_name, colnames(x_matrix))))
  }

  if(config$transform_theta) {
    for(d_par in c("init","lower","upper")) {
      theta_setup[[d_par]] = log(theta_setup[[d_par]])
      if(any(is.na(theta_setup[[d_par]]))) {
        warning("theta ", d_par, "values out of range")
      }
    }
  }

  tmb_data <- list(
    X = x_matrix,
    A =a_matrix,
    y = data_sub[[outcome_var]],
    Q = precision_matrix,
    map = gamma_theta_map,
    elgm_matrix = cc_matrix
  )
  tmb_data <- formatHpolData(tmb_data)

  verbose_orig <- config$verbose
  config$verbose <- config$verbose > 1

  config$beta <- beta_setup$init
  config$theta <- theta_setup$init
  config$gamma <- rep(0, nrow(gamma_setup))

  beta_theta_names = intersect(colnames(beta_setup), colnames(theta_setup))
  beta_theta_names = setdiff(beta_theta_names, "order")

  parameters_info <- list(
    beta = beta_setup,
    gamma = gamma_setup,
    theta = theta_setup,
    parameters = rbind(
      beta_setup[,beta_theta_names],
      theta_setup[,beta_theta_names]
    )
  )

  if (verbose_orig) {
    cat("getting groups...")
  }

  for_groups_x <- tmb_data$XTp
  for_groups_A <- tmb_data$ATp
  for_groups_x@x <- rep(1, length(for_groups_x@x))
  for_groups_A@x <- rep(1, length(for_groups_A@x))

  for_groups <- rbind(for_groups_x, for_groups_A) %*% tmb_data$elgm_matrix
  for_groups@x <- rep(1, length(for_groups@x))
  config$groups <- adlaplace::adFun_groups(
    ATp = for_groups,
    Ngroups = config$Ngroups
  )
  if (verbose_orig) {
    cat("done.")
  }


  cache <- new.env()
  assign("gamma", config$gamma, cache)

  config$opt <- as.list(
    parameters_info$parameters[c("init","lower","upper","parscale")]
  )

  if (for_dev) {
    return(list(
      tmb_data = tmb_data,
      config = config,
      formula = formula,
      terms = terms,
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
    if(verbose_orig) {
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
    mle$hessian <- adaplace::hess(
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
    print(do.call(rbind, config$opt))
  }

  # Add parscale to control if it exists in theta_info
  if (!is.null(config$theta_info) && "parscale" %in% names(config$theta_info)) {
    control$parscale <- config$theta_info$parscale
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
      terms = terms,
      parameters_info = parameters_info,
      random_info = random_info,
      control_inner = control$inner,
      control = control
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

  result$parameters <- try(format_parameters(
    x = result$extra$full_parameters,
    result$objects$parameters_info
  ))

  if (FALSE) {
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
  }

  result$sample <- try(cond_sim_iwp(
    fit = result,
    n = c(result$objects$config$num_sim, 500)[1]
  ))

  result
}
