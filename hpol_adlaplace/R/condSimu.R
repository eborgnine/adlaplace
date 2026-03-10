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

cond_sim_gamma <- function(fit, n) {
  half_h <- fit$extra$halfHinv
  # note tcrossprod(half_h) = Hinv
  ngamma <- nrow(half_h)

  gamma_hat <- fit$opt$solution
  sim_ind <- matrix(rnorm(n * ngamma), ngamma, n)
  sim_gamma_1 <- as.matrix(half_h %*% sim_ind)

  sim_gamma <- sim_gamma_1 + matrix(gamma_hat, length(gamma_hat), ncol(sim_gamma_1))
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
    return(cond_sim_iwp(fit, fit$terms, fit$parameters_info, n, newx))
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
  n
) {
  terms <- fit$objects$terms
  parameters_info <- fit$objects$parameters_info
  fit <- fit$extra

  beta <- fit$fullParameters[seq(1, len = nrow(parameters_info$beta))]

  sim_gamma <- cond_sim_gamma(fit, n)
  rownames(sim_gamma) <- parameters_info$gamma$name

  terms_pred <- get_terms_pred(terms)

  if (!missing(newx)) {
    terms_pred$predDf <- newx
  }

  new_xa <- mapply(
    getNewXA,
    df = terms_pred$predDf,
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

  list(sim = sim_f, x = terms_pred, gamma = sim_gamma, XA = new_xa)
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
