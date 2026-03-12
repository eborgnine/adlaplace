get_terms_pred <- function(terms) {
  svar <- unlist(lapply(terms, "[[", "var"))
  smodel <- unlist(lapply(terms, "[[", "model"))

  sref1 <- lapply(terms, "[[", "ref_value")
  sref <- rep(NA_real_, length(svar))
  sref[unlist(lapply(sref1, length)) > 0] <- unlist(sref1)

  is_hiwp <- which(smodel %in% c("iwp", "hiwp"))
  sref <- sref[is_hiwp]
  svar <- svar[is_hiwp]
  srange <- lapply(terms[is_hiwp], "[[", "range")
  pred_seq <- lapply(srange, function(xx) seq(min(xx), max(xx), len = 100))

  sgroup <- lapply(terms[is_hiwp], "[[", "group_var")

  sref_index <- rep(NA_integer_, length(pred_seq))
  names(pred_seq) <- names(sref_index) <- names(sgroup) <- names(sref) <- svar
  for (d in svar) {
    sref_index[d] <- which.min(abs(sref[d] - pred_seq[[d]]))
  }

  list(
    pred_seq = pred_seq,
    sgroup = sgroup,
    sref = sref,
    sref_index = sref_index
  )
}

build_pred_df <- function(terms_pred) {
  pred_df <- vector("list", length(terms_pred$pred_seq))
  names(pred_df) <- names(terms_pred$pred_seq)

  for (d_var in names(terms_pred$pred_seq)) {
    pred_df_here <- data.frame(x = terms_pred$pred_seq[[d_var]])
    names(pred_df_here) <- d_var

    group_var <- terms_pred$sgroup[[d_var]]
    if (!is.null(group_var) && nzchar(group_var)) {
      pred_df_here[[group_var]] <- "GLOBAL"
    }

    pred_df[[d_var]] <- pred_df_here
  }

  pred_df
}


get_group_effect <- function(
  a_here,
  sim_global_here,
  gamma_info_here,
  sim_gamma,
  d_var,
  d_group,
  probs = c(0.025, 0.5, 0.975)
) {
  gamma_here <- gamma_info_here[
    which(gamma_info_here$group == d_group), ,
    drop = FALSE
  ]
  a_here_new_names <- gsub(
    "_GLOBAL_",
    paste0("_", d_group, "_"),
    colnames(a_here)
  )
  a_here_new_names <- gsub(
    "_fpoly_",
    paste0("_hrpoly_", d_group, "_"),
    a_here_new_names
  )
  matched_cols <- match(gamma_here$name, a_here_new_names)
  matched_rows <- gamma_here$gamma_id + 1L

  if (any(is.na(matched_cols))) {
    warning(
      "missing group-level design columns for ",
      d_var,
      " / ",
      d_group
    )
    return(NULL)
  }

  sim_here_r <- a_here[, matched_cols, drop = FALSE] %*%
    sim_gamma[matched_rows, , drop = FALSE]
  as.matrix(exp(sim_here_r + sim_global_here))
}
get_one_envelope <- function(x, probs) {
  if (requireNamespace("GET", quietly = TRUE)) {
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
  gamma_info,
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
    a_here <- new_xa[[d_var]]$A
    sim_global_here <- sim_f[[d_var]]
    gamma_info_here <- gamma_info[gamma_info$var == d_var, , drop = FALSE]
    d_groups <- setdiff(unique(gamma_info_here$group), c(NA, "GLOBAL"))

    if(is.list(weights)) {
      weights_here = weights[[d_var]]
    } else {
      weights_here = weights
    }
    if(is.null(weights_here)) {
      weights_here = setNames(rep(1/length(d_groups), length(d_groups)),
        d_groups
      )
    } else {
      if(!all(d_groups %in% names(weights_here))) {
        warning(
          "Not all groups have weights, or missing names, giving weight 0."
        )
        weights_here[setdiff(d_groups, names(weights_here))] <- 0
      }
    }

    if (!length(d_groups) || !ncol(a_here)) {
      warning("size mismatch groups and A")
    }

    sim_by_group <- mapply(
      get_group_effect,
      d_group = d_groups,
      MoreArgs = list(
        a_here = a_here,
        sim_global_here = sim_global_here,
        gamma_info_here = gamma_info_here,
        sim_gamma = sim_gamma,
        d_var = d_var,
        probs = probs
      ),
      SIMPLIFY = FALSE
    )
    names(sim_by_group) <- d_groups

    group_quantiles[[d_var]] <- lapply(sim_by_group, function(sim_here) {
      t(apply(sim_here, 1, quantile, probs = probs))
    })

    group_envelope[[d_var]] <- lapply(sim_by_group,
      get_one_envelope, probs = probs_envelope
    )

    var_weights <- weights_here[d_groups]
    weighted_average[[d_var]] <- Reduce(
      `+`,
      Map(function(sim_here, wt) sim_here * wt, sim_by_group, var_weights)
    )

    weighted_envelope[[d_var]] <- get_one_envelope(
      weighted_average[[d_var]],
      prob = probs_envelope
    )
    weighted_quantiles[[d_var]] <- t(
      apply(weighted_average[[d_var]], 1, quantile, probs = probs)
    )
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
  half_h <- fit$extra$deriv$extra$halfH
  # note tcrossprod(half_h) = Hinv
  ngamma <- nrow(half_h)

  gamma_hat <- fit$extra$inner$solution
  sim_ind <- matrix(rnorm(n * ngamma), ngamma, n)
  sim_gamma_1 <- as.matrix(half_h %*% sim_ind)

  sim_gamma <- sim_gamma_1 + matrix(
    gamma_hat,
    length(gamma_hat),
    ncol(sim_gamma_1)
  )
  rownames(sim_gamma) <- names(gamma_hat)

  sim_gamma
}

get_gamma_sim <- function(fit, term, n) {
  gamma_sim <- cond_sim_gamma(fit, n)
  gamma_here <- grep(term, fit$gamma_info$var)
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
  parameters_info <- fit$objects$parameters_info
  gamma_info <- fit$objects$gamma_info
  
  beta <- fit$extra$fullParameters[seq(1, len = nrow(parameters_info$beta))]

  sim_gamma <- hpolcc:::cond_sim_gamma(fit, n)
  rownames(sim_gamma) <- parameters_info$gamma$name

  terms_pred <- hpolcc:::get_terms_pred(terms)

    if (!missing(newx)) {
    terms_pred$pred_df <- newx
  } else {
    terms_pred$pred_df <- hpolcc:::build_pred_df(terms_pred)
  }

    new_xa <- mapply(
    hpolcc:::get_new_xa,
    df = terms_pred$pred_df,
    MoreArgs = list(terms = terms),
    SIMPLIFY = FALSE
  )

      sim_f <- mapply(
    function(xa, gamma, beta) {
      names_both_beta <- intersect(names(beta), colnames(xa$X))
      names_both_gamma <- intersect(rownames(gamma), colnames(xa$A))
      if(length(names_both_beta)) {
      fixed_part <-
        xa$X[, names_both_beta, drop = FALSE] %*%
        beta[names_both_beta]
      } else {
        fixed_part = matrix(0, nrow(xa$X), 1)
      }

      random_part <-
        xa$A[, names_both_gamma, drop = FALSE] %*%
        gamma[names_both_gamma, , drop = FALSE]

      as.matrix(
        random_part + fixed_part[, rep(1, ncol(random_part)),
          drop = FALSE
        ]
      )
    },
    xa = new_xa,
    MoreArgs = list(beta = beta, gamma = sim_gamma),
    SIMPLIFY = FALSE
  )

  result <- hpolcc:::get_group_quantiles(
    sim_f = sim_f,
    new_xa = new_xa,
    sim_gamma = sim_gamma,
    gamma_info = gamma_info,
    weights = weights,
    probs = probs,
    probs_envelope = probs_envelope
  )

  
  result$x <- lapply(terms_pred$pred_df, "[[", 1)
  result$sim <- lapply(sim_f, exp)
  result$quantiles$common <- lapply(result$sim, function(sim_here) {
    t(apply(sim_here, 1, quantile, probs = probs))
  })
  result$envelope$common <- lapply(result$sim,
    get_one_envelope, probs = probs_envelope
  )

  result
}

cond_sim_iid <- function(fit, term, n) {
  terms_here <- grep(term, unlist(lapply(fit$terms, function(xx) xx$var)))
  model_here <- unlist(lapply(fit$terms[terms_here], function(xx) xx$model))

  if (!all(model_here == "iid")) {
    warning("model should be iid to use cond_sim_iid")
  }
  gamma_here <- grep(term, fit$gamma_info$var)
  sx1 <- colnames(fit$obj$env$data$A)[gamma_here]
  sx <- gsub(paste0("(factor[(])?", term, "[)]?"), "", sx1)

  gamma_sim <- get_gamma_sim(fit, term, n)
  list(
    x = sx,
    y = gamma_sim,
    fixed_mean = 0
  )
}
