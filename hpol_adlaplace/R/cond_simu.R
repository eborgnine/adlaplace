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

get_var_weights <- function(weights, d_var, d_groups) {
  if (is.null(weights)) {
    var_weights <- rep(1 / length(d_groups), length(d_groups))
    names(var_weights) <- d_groups
    return(var_weights)
  }

  if (is.list(weights)) {
    var_weights <- weights[[d_var]]
  } else {
    var_weights <- weights
  }

  if (is.null(var_weights)) {
    stop("Missing weights for variable ", d_var)
  }

  if (is.null(names(var_weights))) {
    if (length(var_weights) != length(d_groups)) {
      stop("Unnamed weights for ", d_var, " must match the number of groups")
    }
    names(var_weights) <- d_groups
  }

  var_weights <- var_weights[d_groups]
  if (any(is.na(var_weights))) {
    stop("Weights for ", d_var, " must be named for all groups")
  }

  var_weights / sum(var_weights)
}

get_group_quantiles <- function(
  sim_f,
  new_xa,
  sim_gamma,
  gamma_info,
  weights = NULL
) {
  group_quantiles <- vector("list", length(sim_f))
  weighted_average <- vector("list", length(sim_f))
  weighted_quantiles <- vector("list", length(sim_f))
  names(group_quantiles) <- names(sim_f)
  names(weighted_average) <- names(sim_f)
  names(weighted_quantiles) <- names(sim_f)

  for (d_var in names(sim_f)) {
    a_here <- new_xa[[d_var]]$A
    sim_global_here <- sim_f[[d_var]]
    gamma_info_here <- gamma_info[gamma_info$var == d_var, , drop = FALSE]
    d_groups <- setdiff(unique(gamma_info_here$group), "GLOBAL")

    if (!length(d_groups) || !ncol(a_here)) {
      next
    }

    sim_by_group <- lapply(d_groups, function(d_group) {
      gamma_here <- gamma_info_here[
        gamma_info_here$group == d_group,
        ,
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
      sim_here_r + sim_global_here
    })
    names(sim_by_group) <- d_groups
    sim_by_group <- Filter(Negate(is.null), sim_by_group)
    if (!length(sim_by_group)) {
      next
    }

    d_groups_present <- names(sim_by_group)
    sim_effects <- lapply(sim_by_group, function(sim_here) {
      100 * (exp(sim_here) - 1)
    })

    group_quantiles[[d_var]] <- lapply(sim_effects, function(sim_here) {
      t(apply(sim_here, 1, quantile, probs = c(0.025, 0.5, 0.975)))
    })

    var_weights <- get_var_weights(weights, d_var, d_groups_present)
    weighted_average[[d_var]] <- Reduce(
      `+`,
      Map(function(sim_here, wt) sim_here * wt, sim_effects, var_weights)
    )
    weighted_quantiles[[d_var]] <- t(apply(
      weighted_average[[d_var]],
      1,
      quantile,
      probs = c(0.025, 0.5, 0.975)
    ))
  }

  list(
    group_quantiles = group_quantiles,
    weighted_average = weighted_average,
    weighted_quantiles = weighted_quantiles
  )
}

cond_sim_gamma <- function(fit, n) {
  half_h <- fit$extra$halfHinv
  # note tcrossprod(half_h) = Hinv
  ngamma <- nrow(half_h)

  gamma_hat <- fit$opt$solution
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

#' @export
cond_sim_iwp <- function(
  fit,
  newx,
  n,
  weights = NULL
) {
  terms <- fit$objects$terms
  parameters_info <- fit$objects$parameters_info
  gamma_info <- fit$objects$gamma_info
  fit <- fit$extra

  beta <- fit$fullParameters[seq(1, len = nrow(parameters_info$beta))]

  sim_gamma <- cond_sim_gamma(fit, n)
  rownames(sim_gamma) <- parameters_info$gamma$name

  terms_pred <- get_terms_pred(terms)

  if (!missing(newx)) {
    terms_pred$pred_df <- newx
  } else {
    terms_pred$pred_df <- build_pred_df(terms_pred)
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
      fixed_part <-
        xa$X[, names_both_beta, drop = FALSE] %*%
        beta[names_both_beta]

      random_part <-
        xa$A[, names_both_gamma, drop = FALSE] %*%
        gamma[names_both_gamma, , drop = FALSE]

      random_part + fixed_part[, rep(1, ncol(random_part)), drop = FALSE]
    },
    xa = new_xa,
    MoreArgs = list(beta = beta, gamma = sim_gamma)
  )

  group_summary <- get_group_quantiles(
    sim_f = sim_f,
    new_xa = new_xa,
    sim_gamma = sim_gamma,
    gamma_info = gamma_info,
    weights = weights
  )

  list(
    sim = sim_f,
    x = terms_pred,
    gamma = sim_gamma,
    XA = new_xa,
    group_quantiles = group_summary$group_quantiles,
    weighted_average = group_summary$weighted_average,
    weighted_quantiles = group_summary$weighted_quantiles
  )
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
