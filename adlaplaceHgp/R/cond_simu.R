get_terms_pred <- function(terms, length.out = 100) {
  smodel <- unlist(lapply(terms, class))
  is_iwp <- which(smodel %in% c("iwp", "hiwp", "rsiwp"))
  is_rsiwp <- which(smodel[is_iwp] %in% c("rsiwp"))

  svar <- unlist(lapply(terms[is_iwp], methods::slot, "term"))
  sknots <- lapply(terms[is_iwp], methods::slot, "knots")

  smin <- unlist(lapply(sknots, min))
  smax <- unlist(lapply(sknots, max))

  pred_seq <- mapply(
    function(var, from, to, length.out) {
      result <- data.frame(seq(from = from, to = to, length.out = length.out))
      names(result) <- var
      result
    },
    from = smin, to = smax, var = svar,
    MoreArgs = list(length.out = length.out), SIMPLIFY = FALSE
  )
  names(pred_seq) <- svar

  for (D in is_rsiwp) {
    pred_seq[[D]][[
      terms[is_iwp][[D]]@mult
    ]] <- 1
  }

  pred_seq
}

get_group_effect <- function(
  a_here,
  sim_global_here,
  random_info_here,
  sim_gamma,
  d_var,
  d_group,
  probs = c(0.025, 0.5, 0.975)
) {
  gamma_here <- random_info_here[
    which(random_info_here$by == d_group), ,
    drop = FALSE
  ]
  gamma_here$name_for_a <- gsub("_g[[:digit:]]+$", "", gamma_here$gamma_label)
  gamma_here$name_for_a <- gsub("_hiwp_", "_iwp_", gamma_here$name_for_a)
  gamma_here$name_for_a <- gsub("_hrpoly_", "_rpoly_", gamma_here$name_for_a)

  sim_h_here <- a_here[, gamma_here$name_for_a] %*% sim_gamma[gamma_here$gamma_label, ]
  sim_here <- sim_h_here + sim_global_here

  as.matrix(exp(sim_here))
}

get_one_envelope <- function(x, probs) {
  if (requireNamespace("GET", quietly = TRUE) && !is.null(x)) {
    result <- GET::central_region(
      GET::create_curve_set(
        list(obs = x)
      ),
      probs = probs
    )
    result <- result[, c("lo", "central", "hi")]
  } else {
    result <- NULL
  }
  result
}

normalize_group_weights <- function(weights_here, d_groups) {
  if (is.null(weights_here)) {
    stats::setNames(rep(1 / length(d_groups), length(d_groups)), d_groups)
  } else {
    if (!all(d_groups %in% names(weights_here))) {
      warning(
        "Not all groups have weights, or missing names, giving weight 0."
      )
      weights_here[setdiff(d_groups, names(weights_here))] <- 0
    }
    weights_here
  }
}

compute_sim_by_group <- function(
  a_here,
  sim_global_here,
  random_info_here,
  sim_gamma,
  d_var,
  d_groups,
  probs
) {
  sim_by_group <- mapply(
    get_group_effect,
    d_group = d_groups,
    MoreArgs = list(
      a_here = a_here,
      sim_global_here = sim_global_here,
      random_info_here = random_info_here,
      sim_gamma = sim_gamma,
      d_var = d_var,
      probs = probs
    ),
    SIMPLIFY = FALSE
  )
  stats::setNames(sim_by_group, d_groups)
}

get_one_variable_quantiles <- function(
  d_var,
  sim_f,
  new_xa,
  sim_gamma,
  random_info,
  weights_here,
  probs,
  probs_envelope
) {
  a_here <- new_xa[[d_var]]$random
  sim_global_here <- sim_f[[d_var]]
  random_info_here <- random_info[random_info$term == d_var, , drop = FALSE]
  d_groups <- setdiff(unique(random_info_here$by), NA)

  weights_here <- normalize_group_weights(weights_here, d_groups)

  if (!ncol(a_here)) {
    warning("empty A")
  }

  sim_by_group <- compute_sim_by_group(
    a_here = a_here,
    sim_global_here = sim_global_here,
    random_info_here = random_info_here,
    sim_gamma = sim_gamma,
    d_var = d_var,
    d_groups = d_groups,
    probs = probs
  )

  group_quantiles <- lapply(sim_by_group, function(sim_here) {
    t(apply(sim_here, 1, stats::quantile, probs = probs))
  })

  group_envelope <- lapply(
    sim_by_group,
    get_one_envelope,
    probs = probs_envelope
  )

  var_weights <- weights_here[as.character(d_groups)]
  weighted_average <- Reduce(
    `+`,
    Map(function(sim_here, wt) sim_here * wt, sim_by_group, var_weights)
  )

  weighted_envelope <- NULL
  weighted_quantiles <- NULL
  if (!is.null(weighted_average)) {
    weighted_envelope <- try(get_one_envelope(
      weighted_average,
      probs = probs_envelope
    ))
    weighted_quantiles <- try(t(
      apply(weighted_average, 1, stats::quantile, probs = probs)
    ))
  }

  list(
    group_quantiles = group_quantiles,
    group_envelope = group_envelope,
    weighted_quantiles = weighted_quantiles,
    weighted_envelope = weighted_envelope
  )
}

get_group_quantiles <- function(
  sim_f,
  new_xa,
  sim_gamma,
  random_info,
  weights = NULL,
  probs = c(0.025, 0.5, 0.975),
  probs_envelope = c(0.1, 0.9)
) {
  weighted_quantiles <- vector("list", length(sim_f))
  names(weighted_quantiles) <- names(sim_f)

  group_envelope <- weighted_envelope <-
    group_quantiles <-
    weighted_average <- weighted_quantiles

  for (d_var in names(sim_f)) {
    if (is.list(weights)) {
      weights_here <- weights[[d_var]]
    } else {
      weights_here <- weights
    }

    one_var <- get_one_variable_quantiles(
      d_var = d_var,
      sim_f = sim_f,
      new_xa = new_xa,
      sim_gamma = sim_gamma,
      random_info = random_info,
      weights_here = weights_here,
      probs = probs,
      probs_envelope = probs_envelope
    )

    group_quantiles[[d_var]] <- one_var$group_quantiles
    group_envelope[[d_var]] <- one_var$group_envelope
    weighted_quantiles[[d_var]] <- one_var$weighted_quantiles
    weighted_envelope[[d_var]] <- one_var$weighted_envelope
  }

  list(
    quantiles = list(
      group = group_quantiles,
      weighted = weighted_quantiles
    ),
    envelope = list(
      group = group_envelope,
      weighted = weighted_envelope
    )
  )
}

#' @keywords internal
.as_named_param_vector <- function(x, label_col, value_col) {
  if (is.data.frame(x)) {
    if (!all(c(label_col, value_col) %in% names(x))) {
      stop(
        "expected columns ", label_col, " and ", value_col,
        call. = FALSE
      )
    }
    stats::setNames(x[[value_col]], x[[label_col]])
  } else {
    out <- as.numeric(x)
    if (!is.null(names(x))) {
      names(out) <- names(x)
    }
    out
  }
}

#' @keywords internal
.validate_random_info <- function(random_info) {
  if (!is.data.frame(random_info)) {
    stop("random_info must be a data.frame", call. = FALSE)
  }
  required <- c("gamma_label", "term", "model", "by")
  missing <- setdiff(required, names(random_info))
  if (length(missing) > 0L) {
    stop(
      "random_info must contain columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(random_info)
}

cond_sim_gamma <- function(half_H_inv, gamma_mode, gamma_label, n) {
  if (is.null(half_H_inv)) {
    stop("half_H_inv is required", call. = FALSE)
  }
  ngamma <- nrow(half_H_inv)
  gamma_hat <- .as_named_param_vector(
    gamma_mode,
    label_col = "gamma_label",
    value_col = "mode"
  )
  if (length(gamma_hat) != ngamma) {
    stop(
      "length(gamma_mode) (", length(gamma_hat),
      ") must equal nrow(half_H_inv) (", ngamma, ")",
      call. = FALSE
    )
  }
  if (missing(gamma_label)) {
    gamma_label <- names(gamma_hat)
  }
  if (length(gamma_label) != ngamma) {
    stop("length(gamma_label) must equal nrow(half_H_inv)", call. = FALSE)
  }

  sim_ind <- matrix(stats::rnorm(n * ngamma), ngamma, n)
  sim_gamma_1 <- as.matrix(half_H_inv %*% sim_ind)
  sim_gamma <- sim_gamma_1 + matrix(
    gamma_hat,
    length(gamma_hat),
    ncol(sim_gamma_1)
  )
  rownames(sim_gamma) <- gamma_label
  sim_gamma
}

select_iwp_terms <- function(terms) {
  terms_vars <- lapply(terms, methods::slot, "term")
  terms_no_vars <- unlist(lapply(terms_vars, length)) == 0
  terms_have_vars <- terms[!terms_no_vars]
  terms_vars <- terms_vars[!terms_no_vars]
  terms_type <- unlist(lapply(terms_have_vars, methods::slot, "type"))

  terms_has_by <- vapply(terms_have_vars, function(x) {
    if (methods::.hasSlot(x, "by")) {
      length(methods::slot(x, "by")) > 0
    } else {
      FALSE
    }
  }, logical(1))

  terms_classes <- unlist(lapply(terms_have_vars, class))
  is_iwp <- which(terms_classes %in% c("iwp", "hiwp", "rsiwp"))

  list(
    terms_have_vars = terms_have_vars,
    terms_vars = terms_vars,
    terms_type = terms_type,
    terms_has_by = terms_has_by,
    is_iwp = is_iwp
  )
}

build_iwp_simulations <- function(
  terms_have_vars,
  terms_vars,
  terms_type,
  terms_has_by,
  newx,
  beta_hat,
  sim_gamma
) {
  design_list <- sim_global <- fixed_pred <- sim_f <- list()

  for (D in names(newx)) {
    newx_here <- newx[[D]]
    which_is_D <- terms_vars == D
    which_here <- which(which_is_D & !terms_has_by)
    design_list_here <- mapply(adlaplace::design,
      term = terms_have_vars[which_here],
      MoreArgs = list(data = newx_here), SIMPLIFY = FALSE
    )
    is_beta_here <- which(terms_type[which_here] == "fixed")
    is_gamma_here <- which(terms_type[which_here] == "random")
    design_list[[D]] <- list(
      fixed = do.call(cbind, design_list_here[is_beta_here]),
      random = do.call(cbind, design_list_here[is_gamma_here])
    )
    sim_global[[D]] <-
      design_list[[D]]$random %*% sim_gamma[colnames(design_list[[D]]$random), ]
    if (!is.null(design_list[[D]]$fixed) && ncol(design_list[[D]]$fixed) > 0L) {
      fixed_pred[[D]] <-
        design_list[[D]]$fixed %*% beta_hat[colnames(design_list[[D]]$fixed)]
    } else {
      fixed_pred[[D]] <- rep(0, nrow(newx_here))
    }
    sim_f[[D]] <- sim_global[[D]] + drop(fixed_pred[[D]])
  }

  list(
    design_list = design_list,
    sim_f = sim_f
  )
}

append_common_sim_summaries <- function(result, sim_f, probs, probs_envelope) {
  result$x <- lapply(result$x, "[[", 1)
  result$sim <- lapply(sim_f, exp)
  result$quantiles$common <- lapply(result$sim, function(sim_here) {
    t(apply(sim_here, 1, stats::quantile, probs = probs))
  })
  result$envelope$common <- lapply(
    result$sim,
    get_one_envelope,
    probs = probs_envelope
  )
  result
}

#' Build inputs for conditional IWP simulation
#'
#' Extracts flat arguments for \code{\link{cond_sim_iwp_at}} from the output of
#' \code{\link[adlaplace]{log_lik_laplace}} and \code{\link[adlaplace]{model_data}}.
#'
#' @param laplace Output of \code{log_lik_laplace(..., deriv = TRUE)}.
#' @param model_data Output of \code{model_data()}.
#'
#' @return A list with \code{terms}, \code{random_info}, \code{beta},
#'   \code{gamma_mode}, and \code{half_H_inv}.
#' @export
cond_sim_iwp_inputs <- function(laplace, model_data) {
  if (!is.list(laplace) || is.null(laplace$full_parameters)) {
    stop("laplace must be output of log_lik_laplace(...)", call. = FALSE)
  }
  if (!is.list(model_data) || is.null(model_data$data$info)) {
    stop("model_data must be output of model_data()", call. = FALSE)
  }

  info <- model_data$data$info
  n_beta <- nrow(info$beta)
  n_gamma <- nrow(info$gamma)
  n_theta <- nrow(info$theta)
  full <- as.numeric(laplace$full_parameters)
  n_expected <- n_beta + n_gamma + n_theta
  if (length(full) != n_expected) {
    stop(
      "length(laplace$full_parameters) (", length(full),
      ") must equal n_beta + n_gamma + n_theta (", n_expected, ")",
      call. = FALSE
    )
  }

  beta <- stats::setNames(
    full[seq_len(n_beta)],
    info$beta$beta_label
  )
  gamma_mode <- stats::setNames(
    full[seq(n_beta + 1L, length.out = n_gamma)],
    info$gamma$gamma_label
  )

  list(
    terms = model_data$terms,
    random_info = info$gamma,
    beta = beta,
    gamma_mode = gamma_mode,
    half_H_inv = adlaplace::laplace_half_H_inv(laplace)
  )
}

#' Conditional simulation for IWP/HIWP terms (flat arguments)
#'
#' @description
#' Draw conditional simulations of IWP/HIWP/RSIWP model components from explicit
#' inputs. Prefer \code{\link{cond_sim_iwp}} when you have a Laplace \code{fit}
#' and \code{model_data} objects; use this function for custom pipelines.
#'
#' @param terms Named list of model term objects from \code{model_data()$terms}.
#' @param random_info Data frame from \code{model_data()$data$info$gamma} with
#'   columns \code{gamma_label}, \code{term}, \code{model}, and \code{by}.
#' @param beta Named numeric vector of fixed-effect MLEs, or a data frame with
#'   columns \code{beta_label} and \code{mle}.
#' @param gamma_mode Named numeric vector of random-effect modes at the Laplace
#'   inner optimum, or a data frame with columns \code{gamma_label} and \code{mode}.
#' @param half_H_inv Matrix \eqn{H^{-1/2}} for inner random effects; from
#'   \code{\link[adlaplace]{laplace_half_H_inv}(laplace)}.
#' @param newx Optional list of prediction data frames, one per variable.
#'   When \code{NULL}, a default grid is built from term knot ranges.
#' @param n Number of conditional draws.
#' @param weights Optional group weights (named vector or list by variable).
#' @param probs Quantile levels for summaries.
#' @param probs_envelope Envelope probability levels (requires \pkg{GET}).
#'
#' @return List with simulated curves, quantiles, and envelopes.
#' @export
cond_sim_iwp_at <- function(
  terms,
  random_info,
  beta,
  gamma_mode,
  half_H_inv,
  newx = NULL,
  n = 500,
  weights = NULL,
  probs = c(0.025, 0.5, 0.975),
  probs_envelope = c(0.1, 0.9)
) {
  .validate_random_info(random_info)

  term_info <- select_iwp_terms(terms)
  beta_hat <- .as_named_param_vector(beta, label_col = "beta_label", value_col = "mle")
  gamma_label <- random_info$gamma_label
  sim_gamma <- cond_sim_gamma(half_H_inv, gamma_mode, gamma_label, n)

  if (is.null(newx)) {
    newx <- get_terms_pred(term_info$terms_have_vars[term_info$is_iwp])
  }

  sims <- build_iwp_simulations(
    terms_have_vars = term_info$terms_have_vars,
    terms_vars = term_info$terms_vars,
    terms_type = term_info$terms_type,
    terms_has_by = term_info$terms_has_by,
    newx = newx,
    beta_hat = beta_hat,
    sim_gamma = sim_gamma
  )

  result <- get_group_quantiles(
    sim_f = sims$sim_f,
    new_xa = sims$design_list,
    sim_gamma = sim_gamma,
    random_info = random_info,
    weights = weights,
    probs = probs,
    probs_envelope = probs_envelope
  )

  result$x <- newx
  result$gamma <- sim_gamma
  append_common_sim_summaries(result, sims$sim_f, probs, probs_envelope)
}

#' Conditional simulation for IWP/HIWP terms
#'
#' @description
#' Draw conditional simulations of the Gaussian process model components,
#' then summarize the resulting linear predictors and group-level effect curves.
#'
#' @param fit Output of \code{\link[adlaplace]{log_lik_laplace}(..., deriv = TRUE)}.
#' @param model_data Output of \code{\link[adlaplace]{model_data}()}.
#' @param newx Optional list of prediction data frames, one per variable.
#'   When \code{NULL}, a default prediction grid is built from term knot ranges.
#' @param n Number of conditional draws to simulate.
#' @param weights Optional weights used to average group-level effects.
#' @param probs Numeric vector of probabilities used when computing quantiles.
#' @param probs_envelope Numeric vector of probabilities for envelopes.
#'
#' @return A list with components \code{sim}, \code{quantiles}, and \code{envelope}.
#' @seealso \code{\link{cond_sim_iwp_at}}, \code{\link{cond_sim_iwp_inputs}}
#' @export
cond_sim_iwp <- function(
  fit,
  model_data,
  newx = NULL,
  n = 500,
  weights = NULL,
  probs = c(0.025, 0.5, 0.975),
  probs_envelope = c(0.1, 0.9)
) {
  inputs <- cond_sim_iwp_inputs(fit, model_data)
  cond_sim_iwp_at(
    terms = inputs$terms,
    random_info = inputs$random_info,
    beta = inputs$beta,
    gamma_mode = inputs$gamma_mode,
    half_H_inv = inputs$half_H_inv,
    newx = newx,
    n = n,
    weights = weights,
    probs = probs,
    probs_envelope = probs_envelope
  )
}
