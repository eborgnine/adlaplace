get_terms_pred <- function(terms, length.out=100) {

  smodel <- unlist(lapply(terms, class))
  is_iwp <- which(smodel %in% c("iwp","rsiwp"))
  is_rsiwp <- which(smodel[is_iwp] %in% c("rsiwp"))

  svar <- unlist(lapply(terms[is_iwp], slot, "term"))
  sknots <- lapply(terms[is_iwp], slot, "knots")
  sref <- unlist(lapply(terms[is_iwp], slot, "ref_value"))
  smin = unlist(lapply(sknots, min)) + sref
  smax = unlist(lapply(sknots, max)) + sref

  pred_seq <- mapply(function(var, from, to, length.out) {
    result = data.frame(seq(from=from, to=to, length.out=length.out))
    names(result) = var
    result
  },
    from=smin, to=smax, var = svar, 
    MoreArgs=list(length.out = length.out), SIMPLIFY=FALSE
    )
  names(pred_seq) = svar

  for(D in is_rsiwp) {
    pred_seq[[D]][[
      terms[is_iwp][[D]]@mult
    ]] = 1
  }

  return(pred_seq)
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
  gamma_here$name_for_a = gsub("_g[[:digit:]]+$", "", gamma_here$gamma_label)
  gamma_here$name_for_a = gsub("_hiwp_", "_iwp_", gamma_here$name_for_a)
  gamma_here$name_for_a = gsub("_hrpoly_", "_rpoly_", gamma_here$name_for_a)
  # to do:  won't be rpoly if boundary is fixed. 

  sim_h_here = a_here[,gamma_here$name_for_a] %*% sim_gamma[gamma_here$gamma_label, ]
  sim_here = sim_h_here + sim_global_here

  as.matrix(exp(sim_here))
}
get_one_envelope <- function(x, probs) {
  if (requireNamespace("GET", quietly = TRUE) & !is.null(x)) {
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
    a_here <- new_xa[[d_var]]$random
    sim_global_here <- sim_f[[d_var]]
    random_info_here <- random_info[random_info$term == d_var, , drop = FALSE]
    d_groups <- setdiff(unique(random_info_here$by), NA)

    if (is.list(weights)) {
      weights_here <- weights[[d_var]]
    } else {
      weights_here <- weights
    }
    if (is.null(weights_here)) {
      weights_here <- setNames(
        rep(1 / length(d_groups), length(d_groups)),
        d_groups
      )
    } else {
      if (!all(d_groups %in% names(weights_here))) {
        warning(
          "Not all groups have weights, or missing names, giving weight 0."
        )
        weights_here[setdiff(d_groups, names(weights_here))] <- 0
      }
    }

    if (!ncol(a_here)) {
      warning("empty A")
    }

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
    names(sim_by_group) <- d_groups

    group_quantiles[[d_var]] <- lapply(sim_by_group, function(sim_here) {
      t(apply(sim_here, 1, stats::quantile, probs = probs))
    })

    group_envelope[[d_var]] <- lapply(sim_by_group,
      get_one_envelope,
      probs = probs_envelope
    )

    var_weights <- weights_here[as.character(d_groups)]
    weighted_average[[d_var]] <- Reduce(
      `+`,
      Map(function(sim_here, wt) sim_here * wt, sim_by_group, var_weights)
    )

    if (!is.null(weighted_average[[d_var]])) {
      weighted_envelope[[d_var]] <- try(get_one_envelope(
        weighted_average[[d_var]],
        probs = probs_envelope
      ))
      weighted_quantiles[[d_var]] <- try(t(
        apply(weighted_average[[d_var]], 1, stats::quantile, probs = probs)
      ))
    }
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

cond_sim_gamma <- function(fit, n) {
  half_h <- fit$extra$extra$halfHinv
  # note tcrossprod(half_h) = Hinv
  ngamma <- nrow(half_h)

  gamma_hat <- fit$parameters$gamma$mode
  sim_ind <- matrix(stats::rnorm(n * ngamma), ngamma, n)
  sim_gamma_1 <- as.matrix(half_h %*% sim_ind)

  sim_gamma <- sim_gamma_1 + matrix(
    gamma_hat,
    length(gamma_hat),
    ncol(sim_gamma_1)
  )
  rownames(sim_gamma) <- fit$parameters$gamma$gamma_label

  sim_gamma
}

get_gamma_sim <- function(fit, term, n) {
  gamma_sim <- cond_sim_gamma(fit, n)
  gamma_here <- grep(term, fit$random_info$term)
  gamma_sim[gamma_here, , drop = FALSE]
}

#' @export
cond_sim <- function(fit, term, newx, n = 500) {
  terms_here <- grep(term, unlist(lapply(fit$terms, function(xx) xx$var)))
  model_here <- unlist(lapply(fit$terms[terms_here], function(xx) xx$model))

  if (any(model_here %in% c("hiwp", "iwp"))) {
    return(cond_sim_iwp(fit, newx = newx, n = n))
  }
  if (any(model_here %in% c("iid"))) {
    return(cond_sim_iid(fit, term, n))
  }

  stop("No supported model found for term.")
}

#' Conditional simulation for IWP/HIWP terms
#'
#' @description
#' Draw conditional simulations of the Gaussian process model components,
#' then summarize the
#' resulting linear predictors and group-level effect curves.
#'
#' @param fit A fitted `hnlm()` result object.
#' @param newx Optional list of prediction data frames, one per variable.
#'   When omitted, a default prediction grid is built from the term ranges.
#' @param n Number of conditional draws to simulate.
#' @param weights Optional weights used to average group-level effects.
#'   Supply either a named numeric vector keyed by group or a named list of
#'   such vectors keyed by variable. If `NULL`, equal weights are used.
#' @param probs Numeric vector of probabilities used when computing
#'   quantiles.
#'
#' @return A list with components:
#' \describe{
#'   \item{`sim`}{Simulated linear predictors for each variable.}
#'   \item{`quantiles`}{Pointwise quantiles of `sim` using `probs`.}
#'   \item{`group_quantiles`}{Pointwise quantiles of transformed
#'   group-level effects for each variable and group.}
#'   \item{`weighted_quantiles`}{Pointwise quantiles of the weighted
#'   average group-level effects (on the exponential scale) for each variable.}
#' }
#'
#' @export
cond_sim_iwp <- function(
  fit,
  newx,
  n = 500,
  weights = NULL,
  probs = c(0.025, 0.5, 0.975),
  probs_envelope = c(0.1, 0.9)
) {
  terms <- fit$objects$terms
  terms_vars = lapply(terms, slot, "term")
  terms_no_vars = unlist(lapply(terms_vars, length))==0
  terms_have_vars = terms[!terms_no_vars]
  terms_vars = terms_vars[!terms_no_vars]
  terms_type = unlist(lapply(terms_have_vars, slot, "type"))

  terms_has_by = unlist(lapply(lapply(terms_have_vars, slot, "by"), length))>0
  terms_classes = unlist(lapply(terms_have_vars, class))
  is_iwp = which(terms_classes %in% c("rsiwp","iwp"))

  vars_to_sim <- unique(unlist(terms_vars[is_iwp]))

  sim_gamma <- cond_sim_gamma(fit, n)
  beta_hat = fit$parameters$beta$mle
  names(beta_hat) = fit$parameters$beta$beta_label

  if(missing(newx)) {
    newx = get_terms_pred(terms_have_vars[is_iwp])
  }
  design_list = sim_global = fixed_pred = sim_f = list()
  for(D in names(newx)) {
    newx_here = newx[[D]]
    which_is_D = terms_vars == D
    which_here = which(which_is_D & !terms_has_by)
    design_list_here = mapply(adlaplace::design, 
      term = terms_have_vars[which_here],
      MoreArgs = list(data = newx_here), SIMPLIFY=FALSE)
    is_beta_here = which(terms_type[which_here]  == "fixed")
    is_gamma_here = which(terms_type[which_here]  == "random")
    design_list[[D]] = list(
      fixed = do.call(cbind, design_list_here[is_beta_here]),
      random = do.call(cbind,design_list_here[is_gamma_here])
    )
    sim_global[[D]] = design_list[[D]]$random %*% sim_gamma[colnames(design_list[[D]]$random ), ]
    if(!is.null(design_list[[D]]$fixed )) {
      fixed_pred[[D]] = design_list[[D]]$fixed %*% beta_hat[colnames(design_list[[D]]$fixed)]
    } else {
      fixed_pred[[D]] = rep(0, nrow(newx_here))
    }
    sim_f[[D]] = sim_global[[D]] + drop(fixed_pred[[D]])
  }

  #D=1;matplot(unlist(newx[[D]]),sim_f[[D]], type="l" )

  result <- get_group_quantiles(
    sim_f = sim_f,
    new_xa = design_list,
    sim_gamma = sim_gamma,
    random_info = fit$objects$parameters_info$gamma,
    weights = weights,
    probs = probs,
    probs_envelope = probs_envelope
  )


  result$x <- lapply(newx, "[[", 1)
  result$sim <- lapply(sim_f, exp)
  result$quantiles$common <- lapply(result$sim, function(sim_here) {
    t(apply(sim_here, 1, stats::quantile, probs = probs))
  })
  result$envelope$common <- lapply(result$sim,
    get_one_envelope,
    probs = probs_envelope
  )

  result
}

cond_sim_iid <- function(fit, term, n) {
  terms_here <- grep(term, unlist(lapply(fit$terms, function(xx) xx$var)))
  model_here <- unlist(lapply(fit$terms[terms_here], function(xx) xx$model))

  if (!all(model_here == "iid")) {
    warning("model should be iid to use cond_sim_iid")
  }
  gamma_here <- grep(term, fit$random_info$term)
  sx1 <- colnames(fit$obj$env$data$A)[gamma_here]
  sx <- gsub(paste0("(factor[(])?", term, "[)]?"), "", sx1)

  gamma_sim <- get_gamma_sim(fit, term, n)
  list(
    x = sx,
    y = gamma_sim,
    fixed_mean = 0
  )
}
