#' Fit Hierarchical Non-Linear Models
#'
#' @description
#' This function fits a hierarchical model using the specified formula and data,
#' incorporating case-crossover designs, fixed and random effects, and precision
#' matrices for random effects.
#'
#' @param formula A formula object specifying the model to be fitted.
#' @param data A data frame containing the variables specified in the formula.
#' @param cc_design An object specifying the case-crossover design.
#' @param config A list of configuration options including:
#'        \itemize{
#'          \item dirichlet: Logical; whether to use Dirichlet distribution (default: TRUE)
#'          \item boundary_is_random: Logical; whether boundary should be treated as random (default: FALSE)
#'          \item transform_theta: Logical; whether to transform theta parameters (default: TRUE)
#'        }
#' @param control A list of control options for optimization including:
#'        \itemize{
#'          \item maxit: Maximum number of iterations (default: 1000)
#'          \item trace: Level of tracing output (default: 3)
#'          \item REPORT: Reporting frequency (default: 1)
#'        }
#' @param control_inner A list of control options for inner optimization including:
#'        \itemize{
#'          \item report.level: Reporting level for inner optimization (default: 0)
#'        }
#' @param for_dev Logical; if TRUE, returns intermediate objects for development (default: FALSE).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The function handles fixed effects, random effects, and their associated
#' precision matrices. It also optimizes the model using TMB with options for
#' additional preprocessing and handling specific random effect structures.
#'
#' @return
#' A list containing the fitted model object and related information.
#'
#' @seealso
#' \code{\link{f}} for specifying model terms
#'
#' @examples
#' # Example usage
#' # data <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
#' # result <- hnlm(y ~ x1 + x2, data = data)
#'
#' @import data.table
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

  if (methods::is(formula, "formula")) {
    model_terms <- adlaplace::collect_terms(
      stats::update.formula(formula, . ~ . - 1), # no intercept
      package = "hpolcc", verbose = config$verbose
    )
  } else {
    model_terms <- formula
  }

  if (!any(sapply(model_terms, class) == "overdispersion")) {
    model_terms <- c(
      model_terms,
      adlaplace::overdispersion()
    )
  }

  covariates <- unique(
    unlist(
      lapply(model_terms, methods::slot, "term")
    ),
    value = TRUE, invert = TRUE
  )
  covariates <- unique(unlist(strsplit(covariates, ":")))
  random_slope_terms <- unique(unlist(sapply(model_terms[
    grep("^rs", sapply(model_terms, class))
  ], methods::slot, "mult")))

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
    cat("excluding data colnames:\n")
    print(setdiff(required_vars, colnames(data)))
  }
  # Remove rows with NA values in required variables

  data <- data[stats::complete.cases(data[, required_vars, with = FALSE])]


  if (config$verbose) {
    cat("data has ", nrow(data), " rows\n")
  }

  the_response <- which(
    unlist(lapply(model_terms, function(xx) any(class(xx) == "response")))
  )
  if (length(the_response) != 1) {
    warning("cant find response variable")
  }
  outcome_var <- model_terms[[the_response[1]]]@term

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

  if (!nrow(data_sub)) {
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
    formula = model_terms,
    data = data_sub,
    verbose = config$verbose
  )

  # log transform
  if (config$transform_theta) {
    transform_rows <- model_stuff$info$theta$type == "random"
    transform_cols <- c("init", "lower", "upper")
    model_stuff$info$theta[transform_rows, transform_cols] <- log(
      model_stuff$info$theta[transform_rows, transform_cols]
    )
    all_par_cols <- colnames(model_stuff$info$parameters)
    model_stuff$info$parameters <- rbind(
      model_stuff$info$beta[, all_par_cols],
      model_stuff$info$theta[, all_par_cols]
    )
  }

  # Add case-crossover specific components
  model_stuff$data$elgm_matrix <- cc_matrix
  model_stuff$data$y <- data_sub[[outcome_var]]


  verbose_orig <- config$verbose
  config$verbose <- config$verbose > 1

  config$beta <- model_stuff$info$beta$init
  config$theta <- model_stuff$info$theta$init
  config$gamma <- rep(0, nrow(model_stuff$info$gamma))

  for (D in c("beta", "theta", "gamma")) {
    if (is.null(config[[D]])) {
      config[[D]] <- numeric(0)
    }
  }

  if (verbose_orig) {
    cat("getting groups...")
  }

  config$groups <- adlaplace::adFun_groups(
    ATp = rbind(model_stuff$data$XTp, model_stuff$data$ATp),
    elgm_matrix = model_stuff$data$elgm_matrix,
    Ngroups = config$num_groups,
    min_groups = config$num_threads * 4
  )

  if (verbose_orig) {
    cat("done.")
  }


  cache <- new.env()
  assign("gamma", config$gamma, cache)

  config$opt <- as.list(
    model_stuff$info$parameters[c("init", "lower", "upper", "parscale")]
  )

  if (for_dev) {
    return(list(
      model = model_stuff,
      config = config,
      formula = formula,
      terms = model_terms,
      data = data_sub,
      control = control,
      control_inner = control_inner,
      cache = cache
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
    model_stuff$data,
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
    rownames(to_print) <- model_stuff$info$parameters$label
    print(to_print)
    cat("threads: ", config$num_threads, "\n")
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
      #      tmb_data = model_stuff$data,
      config = config,
      formula = formula,
      terms = model_terms,
      parameters_info = model_stuff$info,
      random_info = random_info,
      control_inner = control$inner,
      control = control,
      cache = cache,
      #      data = data_sub,
      ad_fun = ad_fun
    )
  )
  if (verbose_orig) {
    cat("one last evaulation of likelihhood\n")
  }

  result$extra <- try(adlaplace::logLikLaplace(
    x = result$opt[[grep("solution|par", names(result$opt), value = TRUE)[1]]],
    gamma = result$objects$cache$gamma,
    data = model_stuff$data, # result$objects$tmb_data,
    config = result$objects$config,
    control = result$objects$control_inner,
    adFun = ad_fun,
    deriv = 1
  ))
  result$extra$parameters_orig <- result$parameters
  result$parameters <- try(format_parameters(result))

  if (verbose_orig) {
    cat("hessian of parameters\n")
  }
  result$hessian_parameters <- try(
    Matrix::forceSymmetric(numDeriv::jacobian(
      func = adlaplace::outer_gr,
      x = result$opt[[grep("solution|par", names(result$opt), value = TRUE)[1]]],
      package = "hpolcc",
      data = model_stuff$data, # result$objects$tmb_data,
      config = result$objects$config,
      control_inner = result$objects$control_inner,
      adFun = ad_fun,
      cache = result$objects$cache
    ), "U")
  )
  if (verbose_orig) {
    cat("conditional samples\n")
  }
  result$sample <- try(cond_sim_iwp(
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
  if (any(which_is_iid) & requireNamespace("WoodburyMatrix", quietly = TRUE)) {
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
