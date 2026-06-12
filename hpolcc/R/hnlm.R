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

  config$shards <- adlaplace::ad_shards(
    A = model_data$data$A,
    elgm_matrix = model_data$data$elgm_matrix,
    num_shards = config$num_shards,
    min_groups = min(config$num_shards, config$num_threads * 4L)
  )

  if (verbose_orig) {
    cat("done.\n")
    cat(
      "getting AD fun, ",
      paste(dim(config$shards), collapse = ","), "shards..."
    )
  }

  ad_fun <- adlaplace::ad_fun(
    model_data,
    config,
    num_threads = config$num_threads
  )
  if (verbose_orig) {
    cat("done.\n")
  }

  if (for_dev) {
    return(structure(
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

  optim_result <- try(stats::optim(
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
    control_inner = control_inner,
    hessian = TRUE
  ))

  config$gamma <- get("gamma", cache)

  if (verbose_orig) {
    cat("one last evaluation of likelihood\n")
  }

  par_name <- grep("solution|par", names(optim_result), value = TRUE)[1]
  par_opt <- optim_result[[par_name]]

  laplace <- try(adlaplace::log_lik_laplace(
    x = par_opt,
    gamma = cache$gamma,
    config = config,
    control = control_inner,
    ad_fun = ad_fun,
    deriv = TRUE
  ))

  coefficients <- try(adlaplace::format_parameters(
    info = ad_fun@info,
    gamma = cache$gamma, parameters = par_opt
  ))

  if (verbose_orig) {
    cat("hessian of parameters\n")
  }


  if (verbose_orig) {
    cat("conditional samples\n")
  }
  sample <- try(adlaplaceHgp::cond_sim_iwp(
    fit = laplace,
    model_data = model_data,
    n = c(config$num_sim, 500)[1]
  ))

  hessian <- list(outer = optim_result$hessian, inner = NULL, var_iid = NULL)
  optim_result$hessian <- NULL
  coefficients <- hnlm_attach_parameter_se(coefficients, hessian$outer)


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

  structure(
    list(
      coefficients = coefficients,
      log_lik = laplace$log_lik,
      optim = optim_result,
      converged = isTRUE(optim_result$convergence == 0L),
      extra = laplace,
      hessian = hessian,
      info = ad_fun@info,
      # model_data = model_data,
      #      ad_fun = ad_fun,
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
