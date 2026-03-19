#' Fit (some) hierarchical non-linear models (that we like)
#'
#' @description
#' This function fits a hierarchical model using the specified formula
#' and data, incorporating case-crossover designs, fixed and random
#' effects, and precision matrices for random effects. It allows for
#' flexible inclusion of stratification variables, time variables, and
#' complex random effect structures.
#'
#' @param formula A formula object specifying the model to be fitted.
#' @param data A data frame containing the variables specified in the
#'   formula and any additional variables required for the model.
#' @param cc_design An object specifying the case-crossover design,
#'   including stratification and time variables. Defaults to the output
#'   of `ccDesign()`.
#' @param dirichlet If `TRUE`, fit a dirichlet-multinomial model with a
#'   time-within-strata gamma random effect; otherwise fit a multinomial
#'   model.
#' @param for_dev Logical; if `TRUE`, the function returns intermediate
#'   objects for development purposes. Defaults to `FALSE`.
#' @param verbose logical, if `TRUE` occasional information is printed
#'
#' @return A list containing the fitted TMB object, the formula, terms
#'   used in the model, the case-crossover design, and information about
#'   the gamma and theta parameters.
#' @details The function handles fixed effects, random effects, and
#'   their associated precision matrices. It also optimizes the model
#'   using TMB with options for additional preprocessing and handling
#'   specific random effect structures.
#'
#' @import Matrix
#'
#' @examples
#' # See vignette for basic usage
#'
#' @useDynLib hpolcc
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
    maxit = 2000,
    start.trust.radius = 0.1,
    report.level = 4,
    report.freq = 1,
    report.header.freq = 10,
    report.precision = 5
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
    dirichlet = TRUE,
    dirichlet_init = 0.1,
    Ngroups = 1e4,
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

  terms1 <- collect_terms(formula)
  covariates = unlist(lapply(terms1, "[[", "var"))
  outcome_var <- all.vars(formula)[1]

  strat_time_vars <- unique(c(cc_design$strat_vars, cc_design$time_var))

  strat_time_vars <- strat_time_vars[
    order(sapply(
      data[, strat_time_vars, with = FALSE],
      function(xx) length(unique(xx))
    ), decreasing = FALSE)
  ]
  
  data.table::setDT(data)
  
  required_vars <- unique(c(covariates, strat_time_vars))
  data <- data[complete.cases(data[, ..required_vars])]
  if (anyNA(data[[outcome_var]])) {
    warning("missing values in outcome, treating as zeros")
    data[is.na(get(outcome_var)), (outcome_var) := 0]
  }
  data.table::setorderv(data, strat_time_vars)

  # setup the data for case-crossover
  if (config$verbose) {
    cat("setting strata\n")
  }

  cc_matrix <- setStrata(
    cc_design = cc_design,
    data = data,
    outcome = outcome_var
  )

  if (config$verbose) {
    cat("numer per strata\n")
    print(table(apply(cc_matrix, 2, sum)))
    cat("\ncollecting terms\n")
  }
  # setup of the design matrices and other parameters
  # terms carries all the information throughout
  terms <- lapply(
    terms1,
    get_extra,
    data = data,
    cc_matrix = cc_matrix
  )
  for (term_idx in seq_along(terms)) {
    terms[[term_idx]]$id <- term_idx
  }

  is_fpoly <- which(unlist(lapply(terms, "[[", "model")) == "fpoly")
  is_hrpoly <- which(unlist(lapply(terms, "[[", "model")) == "hrpoly")
  is_asis <- which(unlist(lapply(terms, "[[", "run_as_is")))
  is_random <- setdiff(seq_along(terms), c(is_fpoly, is_hrpoly, is_asis))

  x_as_is <- lapply(terms[is_asis], function(xx, data) {
    res <- Matrix::sparse.model.matrix(xx$f, data, drop.unused.levels = FALSE)
    if (is.factor(data[[xx$var]])) {
      res <- res[, -1, drop = FALSE]
    }
    res
  }, data = data)
  names(x_as_is) <- unlist(lapply(terms[is_asis], "[[", "var"))

  x_fpoly <- lapply(terms[is_fpoly], function(xx, data) {
    x_sub <- as(
      poly(
        data[[xx$var]] - xx$ref_value,
        degree = xx$p,
        raw = TRUE,
        simple = TRUE
      ),
      "TsparseMatrix"
    )
    colnames(x_sub) <- paste0(
      xx$var,
      "_fpoly_",
      seq(from = 1, len = ncol(x_sub))
    )
    x_sub
  }, data = data)
  names(x_fpoly) <- unlist(lapply(terms[is_fpoly], "[[", "var"))

  a_random <- lapply( #parallel::mc
    terms[c(is_hrpoly, is_random)],
    get_design,
    data = data#, mc.cores = config$num_threads
  )

  qs <- lapply(terms[c(is_hrpoly, is_random)], get_precision)
  names(qs) <- names(a_random) <- paste(
    unlist(lapply(terms[is_random], "[[", "var")),
    unlist(lapply(terms[is_random], "[[", "model")),
    sep = "_"
  )

  fpoly_random <- unlist(lapply(terms[is_fpoly], "[[", "boundary_is_random"))
  is_boundary <- is_fpoly[fpoly_random]
  x_fpoly_random <- x_fpoly[fpoly_random]
  x_fpoly_fixed <- x_fpoly[!fpoly_random]

  if (any(fpoly_random) && is.null(config$prec_boundary)) {
    config$prec_boundary <- 0
  }

  # random boundary X's go in A
  q_fpoly <- lapply(
    x_fpoly_random,
    function(xx, prec) Matrix::Diagonal(ncol(xx), x = prec),
    prec = config$prec_boundary
  )
  is_boundary <- is_fpoly
  x_list <- x_as_is

  a_list <- c(x_fpoly_random, a_random)
  q_all <- c(q_fpoly, qs)
  q <- Matrix::bdiag(q_all[!sapply(q_all, is.null)])

  x_list <- c(x_as_is, x_fpoly_fixed)

  # A, Gamma
  if (length(a_list)) {
    a_matrix <- do.call(cbind, a_list) |> as("TsparseMatrix")
  } else {
    a_matrix <- matrix(nrow = 0, ncol = 0) |> as("TsparseMatrix")
  }

  gamma_setup <- lapply(
    terms[c(is_boundary, is_hrpoly, is_random)],
    get_gamma_setup
  )
  gamma_info <- do.call(rbind, gamma_setup)
  gamma_info$global <- as.logical(
    pmin(gamma_info$group == "GLOBAL", TRUE, na.rm = TRUE)
  )

  missing_gamma_info <- setdiff(colnames(a_matrix), gamma_info$name)
  missing_gamma_a <- setdiff(gamma_info$name, colnames(a_matrix))
  if (length(missing_gamma_info)) {
    warning(" missing gamma info ", paste(missing_gamma_info, collapse = ","))
  }
  if (length(missing_gamma_a)) {
    warning(" missing gamma A ", paste(missing_gamma_a, collapse = ","))
  }

  gamma_info <- gamma_info[order(match(gamma_info$name, colnames(a_matrix))), ]
  gamma_info$gamma_id <- seq(0L, len = nrow(gamma_info))
  if (any(is.na(gamma_info$var))) {
    warning("some columns of design matrix not found in gamma")
  }

  if (length(x_list)) {
    x_matrix <- do.call(cbind, x_list)
    beta_info <- data.frame(
      var = rep(names(x_list), unlist(lapply(x_list, ncol))),
      name = colnames(x_matrix)
    )
  } else {
    x_matrix <- matrix(nrow = nrow(data), ncol = 0)
    beta_info <- data.frame()
  }

  # theta setup
  theta_setup <- lapply(
    terms[c(is_hrpoly, is_random)],
    get_theta_setup,
    theta_info = list()
  )

  theta_setup <- c(
    theta_setup,
    list(data.frame(
      var = "overdisp",
      model = "overdisp",
      global = TRUE,
      order = NA,
      init = config$dirichlet_init,
      name = "overdisp"
    ))
  )

  theta_info <- do.call(rbind, theta_setup)
  if (config$transform_theta) {
    theta_info$init <- pmax(-15, log(theta_info$init))
    theta_info$log <- TRUE
  } else {
    theta_info$log <- FALSE
  }
  theta_info$theta_id <- seq(0L, len = nrow(theta_info))

  gamma_theta <- merge(
    gamma_info,
    theta_info,
    by = c("var", "model", "order", "global"),
    all.x = TRUE,
    all.y = TRUE,
    suffixes = c("_gamma", "_theta")
  )

  if (sum(is.na(gamma_theta$theta_id)) > 10) {
    warning("problem matching gamma and theta")
  }

  any_na <- is.na(gamma_theta$theta_id) | is.na(gamma_theta$gamma_id)

  gamma_theta_both <- gamma_theta[!any_na, ]
  # map matrix column theta, row gamma
  gamma_theta_map <- Matrix::sparseMatrix(
    i = gamma_theta_both$gamma_id,
    j = gamma_theta_both$theta_id,
    x = rep(1L, nrow(gamma_theta_both)),
    index1 = FALSE,
    dims = c(nrow(gamma_info), nrow(theta_info))
  )

  tmb_data <- list(
    X = x_matrix,
    A = a_matrix,
    y = data[[all.vars(formula)[1]]],
    Q = q,
    map = gamma_theta_map,
    elgm_matrix = cc_matrix
  )
  tmb_data <- hpolcc:::formatHpolData(tmb_data)
  gamma_info$matchA <- match(gamma_info$name, rownames(tmb_data$ATp))

  verbose_orig <- config$verbose
  config$verbose <- config$verbose > 1

  if (!length(config$beta)) {
    config$beta <- rep(0, nrow(tmb_data$XTp))
  } else {
    config$beta <- rep_len(config$beta, nrow(tmb_data$XTp))
  }
  if (!length(config$theta)) {
    config$theta <- theta_info$init
  } else {
    config$theta <- rep_len(config$theta, length(theta_info$init))
  }
  if (!length(config$gamma)) {
    config$gamma <- rep(0, nrow(tmb_data$ATp))
  } else {
    config$gamma <- rep_len(config$gamma, nrow(tmb_data$ATp))
  }

  parameters_info <- list(
    beta = beta_info,
    gamma = gamma_info,
    theta = theta_info
  )

  if (verbose_orig) {
    cat("getting groups...")
  }

  for_groups <- tmb_data$ATp %*% tmb_data$elgm_matrix
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
  # some checks
  if (!all(parameters_info$gamma$name == rownames(tmb_data$ATp))) {
    warning("names of gamma don't match up")
    setdiff(parameters_info$gamma$name, rownames(tmb_data$ATp))
  }

  if (verbose_orig) {
    cat(for_dev)
  }
  if (for_dev) {
    return(list(
      tmb_data = tmb_data,
      config = config,
      formula = formula,
      terms = terms,
      theta_info = theta_info,
      gamma_info = gamma_info,
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
        paste(dim(config$groups), collapse = ","), "groups"
      )
    )
  }

  ad_fun <- adlaplace::getAdFun(tmb_data, config, package = config$package)

  cache <- new.env(parent = emptyenv())
  cache$gamma <- config$gamma

  x0 <- c(config$beta, config$theta)

  if (verbose_orig) {
    cat("optimizing")
  }
  mle <- trustOptim::trust.optim(
    x = x0,
    fn = adlaplace::outer_fn,
    gr = adlaplace::outer_gr,
    method = "SR1",
    config = config,
    adFun = ad_fun,
    cache = cache,
    control = control,
    control_inner = control_inner
  )

  config$gamma <- get("gamma", cache)

  result <- list(
    opt = mle,
    objects = list(
      tmb_data = tmb_data,
      config = config,
      formula = formula,
      terms = terms,
      parameters_info = parameters_info,
      gamma_info = gamma_info,
      control_inner = control$inner
    )
  )
  if (verbose_orig) {
    cat("done")
  }

  if (FALSE) {
    ad_fun <- adlaplace::getAdFun(
      result$objects$tmb_data,
      result$objects$config,
      package = "hpolcc"
    )
    control_inner <- list(report.level = 0)
  }

  result$extra <- try(adlaplace::logLikLaplace(
    result$opt$solution,
    gamma = result$objects$config$gamma,
    data = result$objects$tmb_data,
    config = result$objects$config,
    control = control_inner,
    adFun = ad_fun,
    deriv = 1
  ))

  result$parameters <- try(hpolcc:::format_parameters(
    x = result$extra$fullParameters,
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
